"""Read-only CIF archive shared through the operating system page cache."""
import mmap
import pickle
from pathlib import Path


class MMapCIFStore:
    def __init__(self, path):
        self.path = Path(path)
        offsets_path = self.path.with_suffix('.offsets.pkl')
        with offsets_path.open('rb') as handle:
            info = pickle.load(handle)  # Only locally generated trusted metadata.
        if info.get('format') != 'tcspv3-cif-mmap-v1':
            raise ValueError('Unsupported CIF mmap offset format')
        self.offsets = info['offsets']
        self.handle = self.path.open('rb')
        try:
            self.map = mmap.mmap(self.handle.fileno(), 0, access=mmap.ACCESS_READ)
            if info['bytes'] != len(self.map):
                raise ValueError('CIF mmap archive size does not match offsets')
        except Exception:
            self.handle.close()
            raise

    def __getitem__(self, mid):
        offset, length = self.offsets[mid]
        if offset < 0 or length < 0 or offset + length > len(self.map):
            raise ValueError('CIF mmap offset outside archive')
        return self.map[offset:offset+length]

    def close(self):
        self.map.close()
        self.handle.close()
