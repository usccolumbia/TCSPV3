import json
from pathlib import Path
import pickle
import subprocess
import sys

from pymatgen.core import Lattice, Structure
from pymatgen.io.cif import CifWriter
from tcspv3.core import Predictor, props

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'scripts'))
from build_template_index import build, index_key


def write_fixture(folder):
    folder.mkdir()
    source = Structure(Lattice.cubic(4), ['Na', 'Cl'], [[0, 0, 0], [.5, .5, .5]])
    source.make_supercell([2, 1, 1])
    (folder / 'donor.cif').write_text(str(CifWriter(source)))
    (folder / 'broken.cif').write_text('not a CIF')
    return source


def test_native_keys_follow_loader_conventions():
    import re
    for formula in ['O2', 'Fe', 'Li2O2', 'Co2Te2', 'Na4Cl4', 'SrTiO3', 'Ca2Yb2B15C']:
        key = index_key(formula)
        decoded = tuple(sorted(float(n or 1) for n in re.findall(r'[A-Z](\d*)', key)))
        assert decoded == props(formula)[1]


def test_primitive_store_loads_and_generates_with_native_predictor(tmp_path):
    import pytest
    source = write_fixture(tmp_path / 'input')
    library = tmp_path / 'library'
    manifest = build(tmp_path / 'input', library)
    assert manifest['accepted'] == 1 and manifest['rejected'] == 1
    metadata = json.loads((library / 'template_metadata.jsonl').read_text())
    assert metadata['original_atoms'] == len(source) == 4
    assert metadata['primitive_atoms'] == 2
    assert 'broken.cif' in (library / 'rejections.jsonl').read_text()
    aux = tmp_path / 'aux'
    (aux / 'data').mkdir(parents=True)
    (aux / 'data/matscholar-embedding.json').write_text(json.dumps({'Na': [0., 0.], 'Cl': [0., 1.], 'K': [1., 0.], 'Br': [1., 1.]}))
    engine = Predictor(aux, library / 'template_index.pkl', cif_store=library / 'data/source_cifs.bin')
    try:
        assert len(engine.source(metadata['template_id'])[0]) == 2
        plan = engine.retrieve('KBr')
        output = engine.generate(plan, 'v3', {}, {}, pool_size=5, topk=2, target_atoms=2)
        assert output['selected'] and output['selected'][0]['nsites'] == 2
        assert not engine.retrieve('NaCl')['formula_count']  # Same composition excluded.
        assert not engine.generate(plan, 'v3', {}, {}, target_atoms=4)['selected']
    finally:
        engine.cif_store.close()
    with pytest.raises(FileExistsError):
        build(tmp_path / 'input', library)


def test_parallel_and_serial_index_and_store_match(tmp_path):
    write_fixture(tmp_path / 'input')
    first, second = tmp_path / 'serial', tmp_path / 'parallel'
    build(tmp_path / 'input', first)
    subprocess.run([sys.executable, str(ROOT / 'scripts/build_template_index.py'),
                    '--cif-dir', str(tmp_path / 'input'), '--outdir', str(second),
                    '--workers', '2'], check=True, capture_output=True, timeout=60)
    for name in ['template_index.pkl', 'data/source_cifs.bin', 'data/source_cifs.offsets.pkl',
                 'template_metadata.jsonl', 'rejections.jsonl']:
        assert (first / name).read_bytes() == (second / name).read_bytes()


def test_exclusions_and_disorder_are_explicit(tmp_path):
    folder = tmp_path / 'input'
    write_fixture(folder)
    (folder / 'skip.cif').write_text((folder / 'donor.cif').read_text())
    disorder = Structure(Lattice.cubic(4), [{'Na': .5, 'K': .5}, 'Cl'], [[0, 0, 0], [.5, .5, .5]])
    (folder / 'disorder.cif').write_text(str(CifWriter(disorder)))
    manifest = build(folder, tmp_path / 'ids_excluded', exclude_ids=['skip'])
    assert manifest['accepted'] == 1 and manifest['excluded'] == 1 and manifest['rejected'] == 2
    import pytest
    with pytest.raises(ValueError, match='No accepted templates'):
        build(folder, tmp_path / 'formula_excluded', exclude_formulas=['Na2Cl2'])
    assert not (tmp_path / 'formula_excluded/BUILD_COMPLETE.json').exists()
