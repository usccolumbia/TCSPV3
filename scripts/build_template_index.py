"""Build a parallel, primitive-cell TCSPV3 index and mmap CIF store from local CIFs."""
import argparse
from collections import defaultdict
import hashlib
import importlib.metadata
import json
import multiprocessing as mp
import os
from pathlib import Path
import pickle
import re
import sys
import warnings

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import numpy as np
from pymatgen.core import Composition
from pymatgen.io.cif import CifParser, CifWriter
from tcspv3.core import props, sha256


def index_key(formula):
    counts = props(formula)[1]  # Match the loader, including pymatgen peroxide/element conventions.
    if len(counts) > 26:
        raise ValueError('Native index supports at most 26 element labels')
    key = ''.join(chr(65 + i) + (str(int(n)) if n != 1 else '') for i, n in enumerate(counts))
    decoded = tuple(sorted(float(n or 1) for n in re.findall(r'[A-Z](\d*)', key)))
    if decoded != counts:
        raise ValueError('Index key does not reproduce native multiplicity signature')
    return key


def cif_paths(root):
    for folder, directories, files in os.walk(root):
        directories.sort()
        for name in sorted(files):
            if name.lower().endswith('.cif') and not name.startswith('._'):
                yield Path(folder) / name


def parse_template(task):
    path, relative, tolerance, excluded_ids, excluded_formulas = task
    metadata = {'source_relative_path': relative, 'source_id': Path(relative).stem}
    try:
        if metadata['source_id'] in excluded_ids:
            return 'excluded', dict(metadata, reason='Excluded source ID'), None
        raw = Path(path).read_bytes()
        metadata['source_sha256'] = hashlib.sha256(raw).hexdigest()
        with warnings.catch_warnings():
            warnings.simplefilter('ignore')
            structures = CifParser.from_str(raw.decode('utf-8')).parse_structures(primitive=False)
        if len(structures) != 1:
            raise ValueError('Require exactly one structure per CIF; split multi-block CIFs first')
        structure = structures[0]
        if not structure.is_ordered:
            raise ValueError('Disordered/partial-occupancy structure; no arbitrary ordering chosen')
        structure.remove_oxidation_states()
        if not structure or structure.volume <= 0 or not np.isfinite(structure.lattice.matrix).all() or not np.isfinite(structure.frac_coords).all():
            raise ValueError('Empty or invalid structure geometry')
        if props(structure.composition.formula)[2] in excluded_formulas:
            return 'excluded', dict(metadata, reason='Excluded reduced composition'), None
        original_atoms = len(structure)
        primitive = structure.get_primitive_structure(tolerance=tolerance)
        formula = primitive.composition.formula.replace(' ', '')
        key = index_key(formula)
        if abs(Composition(formula).num_atoms - len(primitive)) > 1e-7:
            raise ValueError('Primitive formula/site-count mismatch')
        # Preserve actual primitive geometry: no idealized symmetry or energy filtering.
        text = str(CifWriter(primitive, significant_figures=12)).encode('utf-8')
        metadata.update(template_id='template_' + hashlib.sha256(relative.encode()).hexdigest()[:24],
                        original_atoms=original_atoms, primitive_atoms=len(primitive),
                        primitive_cell_formula=formula, index_key=key,
                        exported_sha256=hashlib.sha256(text).hexdigest())
        return 'accepted', metadata, text
    except Exception as exc:
        return 'rejected', dict(metadata, error=type(exc).__name__ + ': ' + str(exc)), None


def build(cif_dir, outdir, workers=1, tolerance=.01, exclude_ids=(), exclude_formulas=(), source_revision=None):
    root, out = Path(cif_dir).resolve(), Path(outdir).resolve()
    if not root.is_dir():
        raise ValueError('CIF directory does not exist: ' + str(root))
    if workers < 1 or tolerance <= 0:
        raise ValueError('Require workers >= 1 and primitive tolerance > 0')
    allocation = os.environ.get('SLURM_CPUS_PER_TASK')
    if allocation and workers > int(allocation):
        raise ValueError('Workers exceed SLURM_CPUS_PER_TASK')
    if out == root or root in out.parents:
        raise ValueError('Output must be outside the input CIF tree')
    out.mkdir(parents=True, exist_ok=False)  # Never replace an existing template library.
    (out / 'data').mkdir()
    exclusions = set(exclude_ids), {props(f)[2] for f in exclude_formulas}
    tasks = ((str(p), p.relative_to(root).as_posix(), tolerance, *exclusions) for p in cif_paths(root))
    counts = dict(accepted=0, excluded=0, rejected=0, scanned=0)
    index, offsets = defaultdict(list), {}
    archive = out / 'data/source_cifs.bin'
    pool = mp.get_context('spawn').Pool(workers) if workers > 1 else None
    results = pool.imap(parse_template, tasks, chunksize=8) if pool else map(parse_template, tasks)
    try:
        with archive.open('xb') as store, (out / 'template_metadata.jsonl').open('x') as metadata_file, (out / 'rejections.jsonl').open('x') as errors:
            for status, metadata, raw in results:
                counts['scanned'] += 1
                counts[status] += 1
                if status == 'accepted':
                    mid = metadata['template_id']
                    if mid in offsets:
                        raise ValueError('Duplicate/colliding template ID: ' + mid)
                    offsets[mid] = (store.tell(), len(raw))
                    store.write(raw)
                    index[metadata['index_key']].append((mid, metadata['primitive_cell_formula']))
                    metadata_file.write(json.dumps(metadata) + '\n')
                else:
                    errors.write(json.dumps(dict(status=status, **metadata)) + '\n')
                if counts['scanned'] % 1000 == 0:
                    print(json.dumps(counts), flush=True)
    finally:
        if pool:
            pool.close()
            pool.join()
    if not counts['accepted']:
        raise ValueError('No accepted templates; inspect rejections.jsonl; build is incomplete')
    index = {key: sorted(rows) for key, rows in sorted(index.items())}
    index_path = out / 'template_index.pkl'
    with index_path.open('xb') as handle:
        pickle.dump(index, handle, protocol=4)
    offsets_path = archive.with_suffix('.offsets.pkl')
    with offsets_path.open('xb') as handle:
        pickle.dump({'format': 'tcspv3-cif-mmap-v1', 'bytes': archive.stat().st_size,
                     'offsets': offsets}, handle, protocol=4)
    manifest = dict(format='tcspv3-local-template-library-v1', **counts,
                    cif_dir=str(root), source_revision=source_revision,
                    primitive_tolerance=tolerance, workers=workers,
                    excluded_ids=sorted(exclusions[0]), excluded_formulas=sorted(exclusions[1]),
                    pymatgen_version=importlib.metadata.version('pymatgen'),
                    structural_deduplication=False, energy_filter=False,
                    index_schema='native anonymous coefficient key -> (template_id, primitive_cell_formula)',
                    files={p.relative_to(out).as_posix(): {'bytes': p.stat().st_size, 'sha256': sha256(p)}
                           for p in [index_path, archive, offsets_path, out / 'template_metadata.jsonl', out / 'rejections.jsonl']})
    assert counts['scanned'] == sum(counts[k] for k in ('accepted', 'excluded', 'rejected'))
    (out / 'BUILD_COMPLETE.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--cif-dir', required=True, type=Path)
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--primitive-tolerance', type=float, default=.01)
    parser.add_argument('--exclude-ids', type=Path, help='JSON list of source filename stems')
    parser.add_argument('--exclude-formulas', type=Path, help='JSON list of excluded compositions')
    parser.add_argument('--source-revision', help='Recorded dataset revision/subset or local export provenance')
    args = parser.parse_args()
    print(json.dumps(build(args.cif_dir, args.outdir, args.workers, args.primitive_tolerance,
                           json.loads(args.exclude_ids.read_text()) if args.exclude_ids else (),
                           json.loads(args.exclude_formulas.read_text()) if args.exclude_formulas else (),
                           args.source_revision), indent=2))


if __name__ == '__main__':
    main()
