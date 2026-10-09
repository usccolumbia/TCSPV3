import json
from pathlib import Path
import sys

import pytest
from pymatgen.core import Structure

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'scripts'))
from export_lemat_cifs import export_rows, row_to_structure


def record():
    return dict(immutable_id='test-123', nperiodic_dimensions=3, dimension_types=[1, 1, 1],
                lattice_vectors=[[4., 0., 0.], [0., 4., 0.], [0., 0., 4.]],
                cartesian_site_positions=[[0., 0., 0.], [2., 2., 2.]], nsites=2,
                species_at_sites=['sodium_site', 'chlorine_site'],
                species=json.dumps([dict(name='sodium_site', chemical_symbols=['Na'], concentration=[1.]),
                                    dict(name='chlorine_site', chemical_symbols=['Cl'], concentration=[1.])]))


def test_cartesian_coordinates_and_species_labels():
    structure = row_to_structure(record())
    assert structure[1].frac_coords.tolist() == [.5, .5, .5]
    assert [s.specie.symbol for s in structure] == ['Na', 'Cl']


def test_disorder_and_nonbulk_are_rejected():
    row = record()
    row['species'] = [dict(name='sodium_site', chemical_symbols=['Na', 'K'], concentration=[.5, .5])]
    with pytest.raises(ValueError, match='occupancy'):
        row_to_structure(row)
    row = record()
    row['dimension_types'] = [1, 1, 0]
    with pytest.raises(ValueError, match='bulk'):
        row_to_structure(row)


def test_export_records_errors_and_keeps_duplicate_ids_in_separate_files(tmp_path):
    good = record()
    bad = record()
    bad['nsites'] = 3
    excluded = dict(record(), immutable_id='holdout')
    result = export_rows([good, good, bad, excluded], tmp_path / 'export', 'LeMaterial/LeMat-BulkUnique',
                         'unique_pbe', 'a' * 40, exclude_ids=['holdout'])
    assert result['accepted'] == 2 and result['rejected'] == 1 and result['excluded'] == 1
    files = sorted((tmp_path / 'export/cifs').glob('*.cif'))
    assert len(files) == 2 and all(len(Structure.from_file(p)) == 2 for p in files)
    with pytest.raises(FileExistsError):
        export_rows([good], tmp_path / 'export', '', '', '')
