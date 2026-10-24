import csv
import argparse
import importlib.util
import json
import os
import re
import sys
import warnings
from pathlib import Path


def quiet():
    """Hide library warnings and Warp's startup messages (use --show-warnings to keep them)."""
    warnings.filterwarnings('ignore')
    os.environ.setdefault('PYTHONWARNINGS', 'ignore')  # also for worker processes
    if importlib.util.find_spec('warp') is not None:
        import warp
        warp.config.quiet = True                      # older Warp
        if hasattr(warp, 'LOG_WARNING'):
            warp.config.log_level = warp.LOG_WARNING  # newer Warp


if '--show-warnings' not in sys.argv:
    quiet()  # before the imports below, which raise some of the warnings
from pymatgen.core import Composition, Structure
from .core import Predictor, dump, sha256, symmetry_vote
from .oxidation import predict
from . import orb, __version__


def row_folder_name(row, formula):
    safe = re.sub(r'[^A-Za-z0-9.+-]+', '_', formula).strip('._-')[:80]
    return f'row_{row}_{safe or "empty"}'


def atom_count(formula):
    amounts = Composition(formula).get_el_amt_dict().values()
    total = sum(amounts)
    if total <= 0 or any(v <= 0 or abs(v-round(v)) > 1e-7 for v in amounts) or abs(total-round(total)) > 1e-7:
        raise ValueError('Full formula requires positive integer atom counts')
    return int(round(total))


def read_formula_csv(path, column=None):
    with Path(path).open(encoding='utf-8-sig', newline='') as handle:
        reader = csv.DictReader(handle)
        headers = reader.fieldnames or []
        if column is not None:
            if column not in headers:
                raise ValueError('CSV column not found: '+column)
        else:
            composition = [h for h in headers if h.strip().lower() == 'composition']
            formulas = [h for h in headers if 'formula' in h.lower()]
            candidates = composition or formulas
            if len(candidates) != 1:
                raise ValueError('Expected one composition/formula column; use --column. Found: '+repr(candidates))
            column = candidates[0]
        records = [(row, (record.get(column) or '').strip()) for row, record in enumerate(reader, 2)]
        if not records:
            raise ValueError('CSV has no data rows')
        return column, records


def clear_previous(root):
    """Remove the files a previous single-formula run wrote, and nothing else."""
    for f in [root/'prediction.json', root/'oxidation.json', *root.glob('tcspv3_*.cif'), *root.glob('native/tcspv3_*.cif')]:
        f.unlink()
    if (root/'native').is_dir() and not any((root/'native').iterdir()):
        (root/'native').rmdir()


BANNER = f"""\
 ==================================================================
   TCSP 3.0   Template-based Crystal Structure Prediction
   Version {__version__}

   Machine Learning and Evolution Laboratory
   University of South Carolina, 2026

   License  MIT
   Code     https://github.com/usccolumbia/TCSPV3
   Data     https://doi.org/10.5281/zenodo.23627915
 =================================================================="""


class HelpFormatter(argparse.ArgumentDefaultsHelpFormatter):
    """Show defaults in --help, except empty ones."""
    def _get_help_string(self, action):
        if action.default in (None, False):
            return action.help
        return super()._get_help_string(action)


def build_parser():
    p = argparse.ArgumentParser(description='Standalone TCSPV3; no CSP dispatcher integration',
                                formatter_class=HelpFormatter)
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--formula', help='One composition, e.g. SrTiO3')
    inputs.add_argument('--csv', help='CSV file of compositions/formulas')
    p.add_argument('--column', help='CSV column override')
    p.add_argument('--database', default=str(Path(__file__).resolve().parent / 'database'), help='Data folder')
    p.add_argument('--index', default=str(Path(__file__).resolve().parent / '731K_index_minus_MP20_test.pkl'), help='Template index')
    p.add_argument('--output', '--out', dest='out', default='./cif_output', help='Output directory')
    p.add_argument('--cif-store', help='Read-only CIF mmap .bin (default) or trusted legacy .pkl')
    p.add_argument('--exclude-formulas', help='JSON list of additional excluded compositions')
    p.add_argument('--mode', choices=['baseline', 'pn', 'union', 'v3'], default='v3', help='Retrieval mode')
    p.add_argument('--formula-type', choices=['full', 'reduced'], default='full',
                   help='full: match input atom count; reduced: allow multiple template-cell Z values')
    p.add_argument('--pool-size', type=int, default=100, help='Valid candidates kept before ranking')
    p.add_argument('--topk', type=int, default=20, help='Number of structures to output')
    p.add_argument('--device', default='cpu', choices=['cpu', 'cuda'], help='Device for BERTOS')
    p.add_argument('--symmetry-vote', action='store_true', help='Experimental')
    p.add_argument('--ranker', choices=['orb-v3', 'native'], default='orb-v3',
                   help='Final stage. orb-v3: relax the candidates with ORB-v3 and rank by ORB-v3 energy per atom; '
                        'native: keep the chemical-score order, no relaxation')
    p.add_argument('--orb-device', default='auto',
                   help='Device for the ORB-v3 stage: auto (least-busy GPU, else CPU), cpu, cuda or cuda:N')
    p.add_argument('--orb-weights', help='Local ORB-v3 weights file (default: downloaded by orb-models)')
    p.add_argument('--orb-steps', type=int, default=orb.STEPS, help='Maximum ORB-v3 relaxation steps')
    p.add_argument('--orb-fmax', type=float, default=orb.FMAX, help='ORB-v3 force threshold in eV/A')
    p.add_argument('--show-warnings', action='store_true', help='Show library warnings')
    p.add_argument('--no-banner', action='store_true', help='Do not print the software banner')
    p.add_argument('--version', action='version', version=f'TCSP 3.0 (TCSPV3 {__version__}), MIT License')
    return p


def main():
    p = build_parser()
    args = p.parse_args()
    if not args.no_banner:
        print(BANNER, file=sys.stderr, flush=True)
    if args.pool_size < 1 or not 0 < args.topk <= args.pool_size:
        p.error('Require pool-size >= topk > 0')
    if not re.fullmatch(r'auto|cpu|cuda(:\d+)?', args.orb_device):
        p.error('--orb-device must be auto, cpu, cuda or cuda:N')
    for required in [Path(args.index), Path(args.database)/'data/matscholar-embedding.json',
                     Path(args.database)/'BERTOS/tokenizer/vocab.txt',
                     Path(args.database)/'BERTOS/trained_models/ICSD_CN/config.json']:
        if not required.exists():
            p.error('Missing local dependency: ' + str(required))
    root = Path(args.out)
    if args.column and not args.csv:
        p.error('--column requires --csv')
    if args.csv:
        try:
            column, records = read_formula_csv(args.csv, args.column)
        except (ValueError, OSError) as exc:
            p.error(str(exc))
        if root.exists() and any(root.iterdir()):
            p.error(f'{root} is not empty; use --output to choose a new folder for this batch')
    else:
        column, records = None, [(None, args.formula)]
        if (root/'prediction.json').exists():
            print(f'Warning: {root} already has a prediction; it will be overwritten. Use --output to write elsewhere.', file=sys.stderr)
            clear_previous(root)
    excluded = json.loads(Path(args.exclude_formulas).read_text()) if args.exclude_formulas else []
    stage = None
    native_k = args.topk
    if args.ranker == 'orb-v3':
        try:
            stage = orb.OrbStage(args.orb_device, weights_path=args.orb_weights, steps=args.orb_steps, fmax=args.orb_fmax)
        except RuntimeError as exc:
            p.error(str(exc))
        native_k = min(max(args.topk, 20), args.pool_size)   # the ORB-v3 stage relaxes at least the 20 native candidates
    engine = Predictor(args.database, args.index, excluded, cif_store=args.cif_store)
    index_hash = sha256(args.index)
    evidence_cache = {}
    summary = []
    for row, formula in records:
        destination = root / row_folder_name(row, formula) if args.csv else root
        try:
            if not formula:
                raise ValueError('Empty formula')
            target_atoms = atom_count(formula) if args.formula_type == 'full' else None
            plan = engine.retrieve(formula)
            needed = set([formula]+plan['embedding']+plan['pn']) - evidence_cache.keys()
            if needed:
                evidence_cache.update(predict(list(needed), args.database, args.device))
            evidence = {f:evidence_cache[f] for f in set([formula]+plan['embedding']+plan['pn'])}
            legacy = {f: v.get('legacy_last_atom', {}) for f, v in evidence.items()}
            result = engine.generate(plan, args.mode, evidence, legacy, args.pool_size, native_k, target_atoms=target_atoms)
            if args.symmetry_vote:
                result['selected'] = symmetry_vote(result['selected'])
            native = result['selected']
            if stage is not None:
                ranked, orb_failures = orb.rank(native, stage.relax, args.topk)
                if not ranked:
                    raise RuntimeError('ORB-v3 stage produced no valid candidate: ' + json.dumps(orb_failures)[:500])
                result['native_selected_ids'] = [c['candidate_id'] for c in native]
                result['selected'] = ranked
                result['orb_v3_failures'] = orb_failures
            result['provenance'] = {'index_sha256': index_hash, 'excluded_formulas': excluded,
                'arguments': vars(args), 'csv_row': row, 'input_formula': formula,
                'csv_column': column, 'formula_type': args.formula_type,
                'target_atoms': target_atoms, 'reference_structure_used': False, 'ranker': args.ranker,
                'final_stage': stage.provenance() if stage is not None else 'native chemical-score order, no relaxation'}
            dump(destination/'prediction.json', result)
            dump(destination/'oxidation.json', evidence)
            for i, c in enumerate(result['selected'], 1):
                Structure.from_dict(c['structure']).to(filename=str(destination/('tcspv3_%03d.cif' % i)))
            if stage is not None:   # keep the unrelaxed native candidates for inspection
                (destination/'native').mkdir(parents=True, exist_ok=True)
                for i, c in enumerate(native, 1):
                    Structure.from_dict(c['structure']).to(filename=str(destination/'native'/('tcspv3_%03d.cif' % i)))
            status = {'row':row, 'formula':formula, 'status':'ok', 'output':str(destination.resolve()),
                      'candidates':len(result['selected']), 'failures':len(result['failures'])}
        except Exception as exc:
            if not args.csv:
                raise
            status = {'row':row, 'formula':formula, 'status':'error', 'error':str(exc), 'candidates':0}
            dump(destination/'error.json', status)
        summary.append(status)
        print(json.dumps(status), flush=True)
        if args.csv:
            dump(root/'batch_summary.json', {'column':column, 'input_csv':str(Path(args.csv).resolve()), 'rows':summary})
    if args.csv:
        with (root/'batch_summary.csv').open('w', newline='') as handle:
            writer = csv.DictWriter(handle, fieldnames=['row','formula','status','candidates','failures','output','error'])
            writer.writeheader();writer.writerows(summary)


if __name__ == '__main__':
    main()
