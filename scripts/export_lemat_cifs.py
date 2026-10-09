"""Stream a revision-pinned official LeMat subset into individual ordered bulk CIFs."""
import argparse
import hashlib
import json
from pathlib import Path
import re

import numpy as np
from pymatgen.core import Element, Structure
from pymatgen.io.cif import CifWriter


def row_to_structure(row):
    if row['nperiodic_dimensions'] != 3 or list(row['dimension_types']) != [1, 1, 1]:
        raise ValueError('Not a fully periodic bulk structure')
    definitions = json.loads(row['species']) if isinstance(row['species'], str) else row['species']
    species = {}
    for definition in definitions:
        symbols, concentrations = definition['chemical_symbols'], definition['concentration']
        if len(symbols) != 1 or len(concentrations) != 1 or concentrations[0] != 1:
            raise ValueError('Mixed/partial occupancy is not supported')
        Element(symbols[0])  # Reject vacancy/unknown element labels.
        if definition['name'] in species:
            raise ValueError('Duplicate species definition')
        species[definition['name']] = symbols[0]
    labels = [species[label] for label in row['species_at_sites']]
    structure = Structure(row['lattice_vectors'], labels, row['cartesian_site_positions'],
                          coords_are_cartesian=True)
    if len(structure) != row['nsites'] or not structure or not structure.is_ordered:
        raise ValueError('Site count/occupancy mismatch')
    if structure.volume <= 0 or not np.isfinite(structure.lattice.matrix).all() or not np.isfinite(structure.frac_coords).all():
        raise ValueError('Invalid geometry')
    return structure


def export_rows(rows, outdir, dataset, subset, revision, max_records=1000, exclude_ids=()):
    if max_records < 0:
        raise ValueError('max_records must be nonnegative; zero means all records')
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=False)
    (out / 'cifs').mkdir()
    excluded_ids = set(exclude_ids)
    counts = dict(scanned=0, accepted=0, excluded=0, rejected=0)
    with (out / 'export_metadata.jsonl').open('x') as metadata, (out / 'rejections.jsonl').open('x') as errors:
        for number, row in enumerate(rows):
            if max_records and number >= max_records:
                break
            counts['scanned'] += 1
            record = dict(record_number=number, immutable_id=row.get('immutable_id'),
                          functional=row.get('functional'), last_modified=row.get('last_modified'))
            if row.get('immutable_id') in excluded_ids:
                counts['excluded'] += 1
                errors.write(json.dumps(dict(status='excluded', reason='Excluded upstream ID', **record)) + '\n')
                continue
            try:
                structure = row_to_structure(row)
                safe = re.sub(r'[^A-Za-z0-9_-]+', '_', str(row.get('immutable_id') or 'unknown'))[:80]
                filename = f'lemat_{number:09d}_{safe}.cif'
                text = str(CifWriter(structure, significant_figures=12)).encode()
                (out / 'cifs' / filename).write_bytes(text)
                counts['accepted'] += 1
                metadata.write(json.dumps(dict(file='cifs/' + filename, sha256=hashlib.sha256(text).hexdigest(),
                                               actual_cell_formula=structure.composition.formula,
                                               nsites=len(structure), **record)) + '\n')
            except Exception as exc:
                counts['rejected'] += 1
                errors.write(json.dumps(dict(status='rejected', error=type(exc).__name__ + ': ' + str(exc), **record)) + '\n')
            if counts['scanned'] % 1000 == 0:
                print(json.dumps(counts), flush=True)
    assert counts['scanned'] == counts['accepted'] + counts['excluded'] + counts['rejected']
    if not counts['accepted']:
        raise ValueError('No accepted records; inspect rejections.jsonl')
    manifest = dict(dataset=dataset, subset=subset, revision=revision, split='train',
                    max_records=max_records, **counts, excluded_ids=sorted(excluded_ids),
                    note='HF split name is not a leakage-controlled benchmark split; no energy filtering')
    (out / 'EXPORT_COMPLETE.json').write_text(json.dumps(manifest, indent=2) + '\n')
    return manifest


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dataset', choices=['LeMaterial/LeMat-Bulk', 'LeMaterial/LeMat-BulkUnique'],
                        default='LeMaterial/LeMat-BulkUnique')
    parser.add_argument('--subset', default='unique_pbe')
    parser.add_argument('--revision', required=True, help='Frozen 40-character Hugging Face commit SHA')
    parser.add_argument('--outdir', required=True, type=Path)
    parser.add_argument('--max-records', type=int, default=1000, help='Small pilot by default; explicitly use 0 for all')
    parser.add_argument('--exclude-ids', type=Path, help='JSON list of original immutable_id values')
    args = parser.parse_args()
    if not re.fullmatch(r'[0-9a-fA-F]{40}', args.revision):
        parser.error('--revision must be a full Hugging Face commit SHA, not a mutable branch name')
    from datasets import load_dataset
    rows = load_dataset(args.dataset, args.subset, split='train', revision=args.revision, streaming=True)
    print(json.dumps(export_rows(rows, args.outdir, args.dataset, args.subset, args.revision,
                                args.max_records, json.loads(args.exclude_ids.read_text()) if args.exclude_ids else ()), indent=2))


if __name__ == '__main__':
    main()
