"""Parallel CSV TCSPV3 prediction with persistent per-process model/data caches."""
import csv
import json
import math
import multiprocessing as mp
import os
import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from pathlib import Path

if __package__ in (None, ''):
    sys.path.insert(0, str(Path(__file__).resolve().parent))
from tcspv3.__main__ import atom_count, build_parser, read_formula_csv, row_folder_name
from tcspv3.core import Predictor, dump, sha256, symmetry_vote
from tcspv3.oxidation import predict
from pymatgen.core import Structure

_CONTEXT = None


def default_workers():
    import psutil
    idle = sum(1 - x / 100 for x in psutil.cpu_percent(interval=1, percpu=True))
    cpu_limit = max(1, math.floor(idle * .6))
    # The mmap archive is shared; allow RAM for offsets, index, Python, and BERTOS per worker.
    memory_limit = max(1, math.floor((psutil.virtual_memory().available - 2e9) / 0.8e9))
    return min(cpu_limit, memory_limit), {'idle_cpu_equivalents': round(idle, 2),
        'cpu_limit': cpu_limit, 'memory_limit': memory_limit}


def init_worker(options):
    global _CONTEXT
    os.environ['OMP_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'
    import torch
    torch.set_num_threads(1)
    engine = Predictor(options['database'], options['index'], options['excluded'],
                       cif_store=options['cif_store'])
    _CONTEXT = {'options': options, 'engine': engine, 'evidence_cache': {},
                'index_hash': options['index_hash']}


def run_record(record):
    row, formula = record
    ctx = _CONTEXT
    o = ctx['options']
    destination = Path(o['out']) / row_folder_name(row, formula)
    try:
        if not formula:
            raise ValueError('Empty formula')
        target_atoms = atom_count(formula) if o['formula_type'] == 'full' else None
        plan = ctx['engine'].retrieve(formula)
        needed = set([formula] + plan['embedding'] + plan['pn'])
        missing = needed - ctx['evidence_cache'].keys()
        if missing:
            ctx['evidence_cache'].update(predict(list(missing), o['database'], o['device']))
        evidence = {f: ctx['evidence_cache'][f] for f in needed}
        legacy = {f: v.get('legacy_last_atom', {}) for f, v in evidence.items()}
        result = ctx['engine'].generate(plan, o['mode'], evidence, legacy,
             o['pool_size'], o['topk'], target_atoms=target_atoms)
        if o['symmetry_vote']:
            result['selected'] = symmetry_vote(result['selected'])
        result['provenance'] = {'index_sha256': ctx['index_hash'],
            'excluded_formulas': o['excluded'], 'arguments': o['arguments'],
            'csv_row': row, 'input_formula': formula, 'csv_column': o['column'],
            'formula_type': o['formula_type'], 'target_atoms': target_atoms,
            'reference_structure_used': False}
        dump(destination/'prediction.json', result)
        dump(destination/'oxidation.json', evidence)
        for i, c in enumerate(result['selected'], 1):
            Structure.from_dict(c['structure']).to(filename=str(destination/('tcspv3_%03d.cif' % i)))
        return {'row': row, 'formula': formula, 'status': 'ok',
                'output': str(destination.resolve()), 'candidates': len(result['selected']),
                'failures': len(result['failures'])}
    except Exception as exc:
        status = {'row': row, 'formula': formula, 'status': 'error',
                  'error': str(exc), 'candidates': 0}
        dump(destination/'error.json', status)
        return status


def main():
    p = build_parser()
    p.description = 'Parallel CSV TCSPV3 prediction; one model/index load per worker'
    p.add_argument('--workers', type=int, help='Worker processes; default 60%% of idle CPU capacity, capped by available RAM')
    args = p.parse_args()
    if args.formula or not args.csv:
        p.error('tcspx.py requires --csv; use tcsp.py for one formula')
    if args.device != 'cpu':
        p.error('Parallel mode currently supports --device cpu only')
    if args.workers is not None and args.workers < 1:
        p.error('--workers must be at least 1')
    if args.pool_size < 1 or not 0 < args.topk <= args.pool_size:
        p.error('Require pool-size >= topk > 0')
    for required in [Path(args.index), Path(args.database)/'data/matscholar-embedding.json',
                     Path(args.database)/'BERTOS/tokenizer/vocab.txt',
                     Path(args.database)/'BERTOS/trained_models/ICSD_CN/config.json']:
        if not required.exists():
            p.error('Missing local dependency: ' + str(required))
    try:
        column, records = read_formula_csv(args.csv, args.column)
    except (ValueError, OSError) as exc:
        p.error(str(exc))
    root = Path(args.out)
    if root.exists() and any(root.iterdir()):
        p.error('Batch output directory must be empty or new')
    root.mkdir(parents=True, exist_ok=True)
    default, resource = default_workers()
    workers = min(args.workers or default, len(records))
    if args.workers and args.workers > resource['memory_limit']:
        print('Warning: requested workers exceed estimated RAM-safe limit', file=sys.stderr)
    excluded = json.loads(Path(args.exclude_formulas).read_text()) if args.exclude_formulas else []
    options = vars(args).copy()
    options.update(excluded=excluded, index_hash=sha256(args.index), column=column)
    options['arguments'] = vars(args).copy()
    dump(root/'batch_config.json', {'workers': workers, 'worker_selection': resource,
        'input_csv': str(Path(args.csv).resolve()), 'column': column, 'options': vars(args)})
    summary = {}
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context('spawn'),
                             initializer=init_worker, initargs=(options,)) as executor:
        futures = {executor.submit(run_record, record): record[0] for record in records}
        for future in as_completed(futures):
            result = future.result()
            summary[result['row']] = result
            print(json.dumps(result), flush=True)
            ordered = [summary[n] for n in sorted(summary)]
            dump(root/'batch_summary.json', {'column': column, 'input_csv': str(Path(args.csv).resolve()),
                                              'workers': workers, 'rows': ordered})
    with (root/'batch_summary.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=['row','formula','status','candidates','failures','output','error'])
        writer.writeheader();writer.writerows([summary[n] for n in sorted(summary)])


if __name__ == '__main__':
    main()
