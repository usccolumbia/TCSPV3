import csv
import argparse
import json
import re
from pathlib import Path
from pymatgen.core import Composition, Structure
from .core import Predictor, dump, sha256, symmetry_vote
from .oxidation import predict


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


def build_parser():
    p = argparse.ArgumentParser(description='Standalone TCSPV3; no CSP dispatcher integration')
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--formula')
    inputs.add_argument('--csv', help='CSV file of compositions/formulas')
    p.add_argument('--column', help='CSV column override')
    p.add_argument('--database', default=str(Path(__file__).resolve().parent / 'database'))
    p.add_argument('--index', default=str(Path(__file__).resolve().parent / '731K_index_minus_MP20_test.pkl'))
    p.add_argument('--output', '--out', dest='out', default='./cif_output',
                   help='Output directory (default: ./cif_output)')
    p.add_argument('--cif-store', help='Read-only CIF mmap .bin (default) or trusted legacy .pkl')
    p.add_argument('--exclude-formulas', help='JSON list of additional excluded compositions')
    p.add_argument('--mode', choices=['baseline', 'pn', 'union', 'v3'], default='v3')
    p.add_argument('--formula-type', choices=['full', 'reduced'], default='full',
                   help='full: match input atom count (default); reduced: allow multiple template-cell Z values')
    p.add_argument('--pool-size', type=int, default=100)
    p.add_argument('--topk', type=int, default=20)
    p.add_argument('--device', default='cpu', choices=['cpu', 'cuda'])
    p.add_argument('--symmetry-vote', action='store_true', help='Experimental, off by default')
    return p


def main():
    p = build_parser()
    args = p.parse_args()
    if args.pool_size < 1 or not 0 < args.topk <= args.pool_size:
        p.error('Require pool-size >= topk > 0')
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
            p.error('Batch output directory must be empty or new')
    else:
        column, records = None, [(None, args.formula)]
        if (root/'prediction.json').exists():
            p.error('Output already contains a prediction; choose a new directory')
    excluded = json.loads(Path(args.exclude_formulas).read_text()) if args.exclude_formulas else []
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
            result = engine.generate(plan, args.mode, evidence, legacy, args.pool_size, args.topk, target_atoms=target_atoms)
            if args.symmetry_vote:
                result['selected'] = symmetry_vote(result['selected'])
            result['provenance'] = {'index_sha256': index_hash, 'excluded_formulas': excluded,
                'arguments': vars(args), 'csv_row': row, 'input_formula': formula,
                'csv_column': column, 'formula_type': args.formula_type,
                'target_atoms': target_atoms, 'reference_structure_used': False}
            dump(destination/'prediction.json', result)
            dump(destination/'oxidation.json', evidence)
            for i, c in enumerate(result['selected'], 1):
                Structure.from_dict(c['structure']).to(filename=str(destination/('tcspv3_%03d.cif' % i)))
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
