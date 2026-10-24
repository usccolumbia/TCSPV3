import hashlib
import importlib.util
import io
from pathlib import Path
import tempfile
from unittest.mock import patch

import pytest

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('install_data', ROOT / 'scripts/install_data.py')
installer = importlib.util.module_from_spec(spec)
spec.loader.exec_module(installer)


def release(data=b'test bundle', url='https://example.org/bundle.tar.gz'):
    return {'filename': 'bundle.tar.gz', 'bytes': len(data), 'sha256': hashlib.sha256(data).hexdigest(), 'download_url': url}


def test_download_verifies_and_existing_good_file_needs_no_network():
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'bundle.tar.gz'
        with patch.object(installer, 'urlopen', return_value=io.BytesIO(b'test bundle')):
            assert installer.download(release(), output) == output
        assert output.read_bytes() == b'test bundle'
        with patch.object(installer, 'urlopen', side_effect=AssertionError('Unexpected download')):
            installer.download(release(), output)


def test_wrong_hash_is_not_kept():
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'bundle.tar.gz'
        with patch.object(installer, 'urlopen', return_value=io.BytesIO(b'bad!! bundle')):
            with pytest.raises(ValueError, match='SHA256'):
                installer.download(release(), output)
        assert list(Path(folder).iterdir()) == []


def test_missing_url_and_different_existing_file_are_refused():
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'bundle.tar.gz'
        with pytest.raises(ValueError, match='download_url'):
            installer.download(release(url=None), output)
        output.write_bytes(b'other bytes')
        with pytest.raises(ValueError, match='differs'):
            installer.download(release(), output)
        assert output.read_bytes() == b'other bytes'


def test_manifest_records_a_pinned_bundle():
    assert installer.RELEASE['filename'] == 'TCSPV3-v3.0.0-template-assets.tar.gz'
    assert installer.RELEASE['asset_count'] == len(installer.EXPECTED) == 10
    assert len(installer.RELEASE['sha256']) == 64
