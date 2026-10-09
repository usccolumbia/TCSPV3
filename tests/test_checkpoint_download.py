import hashlib
import importlib.util
import io
import json
from pathlib import Path
import tempfile
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location('download_mattersim', ROOT / 'scripts/download_mattersim.py')
downloader = importlib.util.module_from_spec(spec)
spec.loader.exec_module(downloader)


def expected(data=b'test checkpoint'):
    return {'url': 'https://example.org/model.pth', 'bytes': len(data),
            'sha256': hashlib.sha256(data).hexdigest()}


def test_download_verifies_and_existing_good_file_needs_no_network():
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'model.pth'
        with patch.object(downloader, 'urlopen', return_value=io.BytesIO(b'test checkpoint')):
            assert downloader.download(expected(), output) == output
        assert output.read_bytes() == b'test checkpoint'
        with patch.object(downloader, 'urlopen', side_effect=AssertionError('Unexpected download')):
            downloader.download(expected(), output)


def test_wrong_hash_is_not_installed():
    import pytest
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'model.pth'
        with patch.object(downloader, 'urlopen', return_value=io.BytesIO(b'bad! checkpoint')):
            with pytest.raises(ValueError, match='SHA256'):
                downloader.download(expected(), output)
        assert not output.exists()
        assert list(Path(folder).iterdir()) == []


def test_existing_different_checkpoint_is_not_overwritten():
    import pytest
    with tempfile.TemporaryDirectory() as folder:
        output = Path(folder) / 'model.pth'
        output.write_bytes(b'other')
        with patch.object(downloader, 'urlopen', side_effect=AssertionError('Unexpected download')):
            with pytest.raises(ValueError, match='size'):
                downloader.download(expected(), output)
        assert output.read_bytes() == b'other'


def test_mattersim_is_external_and_benchmark_model_is_unchanged():
    models = json.loads((ROOT / 'model_manifest.json').read_text())
    assert models['benchmark_model'] == '1M'
    assert models['models']['5M']['url'] == (
        'https://github.com/microsoft/mattersim/raw/refs/heads/main/'
        'pretrained_models/mattersim-v1.0.0-5M.pth')
    assets = json.loads((ROOT / 'asset_manifest.json').read_text())['assets']
    assert len(assets) == 10
    assert not any('mattersim' in name for name in assets)
