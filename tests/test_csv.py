import tempfile
from pathlib import Path
import pytest
from tcspv3.__main__ import read_formula_csv

@pytest.mark.parametrize('column',['composition','Composition','reduced_formula','FullFormula'])
def test_detection(column):
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'in.csv';p.write_text('\ufeff'+column+',id\nSrTiO3,1\nSrTiO3,2\n,3\n')
        assert read_formula_csv(p)==(column,[(2,'SrTiO3'),(3,'SrTiO3'),(4,'')])

def test_ambiguous_and_override():
    with tempfile.TemporaryDirectory() as d:
        p=Path(d)/'in.csv';p.write_text('formula,reduced_formula\nNaCl,NaCl\n')
        with pytest.raises(ValueError):read_formula_csv(p)
        assert read_formula_csv(p,'formula')[1]==[(2,'NaCl')]
