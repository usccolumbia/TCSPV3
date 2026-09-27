"""BERTOS evidence without collapsing mixed predictions to the last atom."""
from functools import lru_cache
import itertools
import math
import numpy as np
from pymatgen.core import Composition, Element


def neutral(formula, states):
    d = Composition(formula).reduced_composition.get_el_amt_dict()
    return bool(states) and set(d) == set(states) and abs(sum(n * states[e] for e, n in d.items())) < 1e-6


def evidence_from_probabilities(formula, probabilities, max_assignments=4):
    d = Composition(formula).reduced_composition.get_el_amt_dict()
    if all(Element(e).is_metal for e in d):
        return {'assignments': [], 'neutral_mass': 0., 'reason': 'all_metal'}
    elements = sorted(d)
    if set(elements) != set(probabilities):
        return {'assignments': [], 'neutral_mass': 0., 'reason': 'missing_species'}
    options = [sorted(probabilities[e].items(), key=lambda x: (-x[1], x[0]))[:3] for e in elements]
    found = []
    for choice in itertools.product(*options):
        ox = {e: int(v[0]) for e, v in zip(elements, choice)}
        if neutral(formula, ox):
            # Independent species probabilities are heuristic evidence, not calibrated confidence.
            mass = math.prod(float(v[1]) for v in choice)
            found.append({'states': ox, 'mass': mass})
    found.sort(key=lambda x: (-x['mass'], tuple(x['states'].values())))
    mass = sum(x['mass'] for x in found)
    kept = found[:max_assignments]
    total = sum(x['mass'] for x in kept)
    for x in kept:
        x['weight'] = x['mass'] / total if total else 0.
    return {'assignments': kept, 'neutral_mass': min(1., mass),
            'reason': 'neutral_evidence' if kept else 'no_neutral_top3_assignment'}


def legacy_evidence(formula, states):
    return {'assignments': [{'states': states, 'weight': 1., 'mass': 1.}], 'neutral_mass': 1.} if neutral(formula, states) else {'assignments': [], 'neutral_mass': 0.}


def predict(formulas, database, device='cpu'):
    import torch
    tok, model = load_bertos(str(database), device)
    output = {}
    formulas = sorted(set(formulas))
    for start in range(0, len(formulas), 64):
        fs = formulas[start:start+64]
        sequences = [[e for e, n in Composition(f).reduced_composition.get_el_amt_dict().items() for _ in range(int(n))] for f in fs]
        texts = [' '.join(s) for s in sequences]
        inputs = tok(texts, padding=True, return_tensors='pt')
        if inputs['input_ids'].shape[1] > model.config.max_position_embeddings:
            raise ValueError('Formula exceeds BERTOS context length')
        with torch.inference_mode():
            prob = model(**{k: v.to(device) for k, v in inputs.items()}).logits.softmax(-1).cpu().numpy()
        for f, seq, ids, p in zip(fs, sequences, inputs['input_ids'], prob):
            if int((ids != tok.pad_token_id).sum()) != len(seq)+2 or tok.unk_token_id in ids.tolist():
                output[f] = {'assignments': [], 'neutral_mass': 0., 'reason': 'token_coverage_error'}
                continue
            species = {e: np.mean([p[i+1] for i, s in enumerate(seq) if s == e], axis=0) for e in set(seq)}
            probabilities = {e: {int(i)-5: float(v) for i, v in enumerate(pp)} for e, pp in species.items()}
            output[f] = evidence_from_probabilities(f, probabilities)
            output[f]['species_probabilities'] = probabilities
            output[f]['legacy_last_atom'] = {e: int(p[i+1].argmax())-5 for i, e in enumerate(seq)}
        if start % 1024 == 0:
            print('BERTOS', start, '/', len(formulas), flush=True)
    return output


@lru_cache(maxsize=1)
def load_bertos(database, device):
    import torch
    from transformers import AutoModelForTokenClassification, BertTokenizerFast
    from pathlib import Path
    root = Path(database) / 'BERTOS'
    tok = BertTokenizerFast.from_pretrained(str(root/'tokenizer'), do_lower_case=False)
    model = AutoModelForTokenClassification.from_pretrained(str(root/'trained_models/ICSD_CN')).eval().to(device)
    return tok, model
