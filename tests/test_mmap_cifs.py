import pickle
from pathlib import Path
import pytest
from tcspv3.mmap_cifs import MMapCIFStore


def test_mmap_lookup_and_bounds(tmp_path):
    data = tmp_path/'source_cifs.bin'
    data.write_bytes(b'one-two')
    data.with_suffix('.offsets.pkl').write_bytes(pickle.dumps({
        'format':'tcspv3-cif-mmap-v1', 'bytes':7,
        'offsets':{'a':(0,3),'b':(3,4),'bad':(6,9)}}))
    store = MMapCIFStore(data)
    try:
        assert store['a'] == b'one'
        assert store['b'] == b'-two'
        with pytest.raises(ValueError):store['bad']
        with pytest.raises(KeyError):store['absent']
    finally:
        store.close()
