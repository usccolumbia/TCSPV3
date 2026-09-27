from tcspv3.__main__ import atom_count
from tcspv3.core import props
import pytest

def test_full_vs_reduced_counts():
    assert atom_count('Fe4O6') == 10
    assert atom_count('Fe2O3') == 5
    assert props('Fe4O6')[2] == props('Fe2O3')[2]
    with pytest.raises(ValueError):atom_count('Fe0.5O')
