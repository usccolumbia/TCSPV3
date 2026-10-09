"""Download an official MatterSim checkpoint, verifying size and SHA256."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import tempfile
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]


def verify(path, expected):
    if path.stat().st_size != expected['bytes']:
        raise ValueError('Checkpoint size mismatch: ' + str(path))
    digest = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(4 * 1024 * 1024), b''):
            digest.update(chunk)
    if digest.hexdigest() != expected['sha256']:
        raise ValueError('Checkpoint SHA256 mismatch: ' + str(path))


def download(expected, destination):
    destination = Path(destination)
    if destination.exists():
        verify(destination, expected)  # Never overwrite a different existing model.
        return destination
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(expected['url'], headers={'User-Agent': 'TCSPV3-checkpoint-downloader'})
    with tempfile.NamedTemporaryFile(dir=destination.parent, suffix='.partial', delete=False) as output:
        temporary = Path(output.name)
        try:
            with urlopen(request, timeout=60) as response:
                total = 0
                while chunk := response.read(4 * 1024 * 1024):
                    total += len(chunk)
                    if total > expected['bytes']:
                        raise ValueError('Checkpoint exceeds expected size')
                    output.write(chunk)
            output.flush()
            verify(temporary, expected)
            # Exclusive publication prevents replacing a concurrently installed model.
            os.link(temporary, destination)
        finally:
            temporary.unlink(missing_ok=True)
    return destination


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--model', choices=['1M', '5M'], default='1M',
                        help='1M reproduces the recorded benchmark; 5M is listed by Matbench Discovery')
    parser.add_argument('--output', type=Path, help='Default: checkpoints/<official filename>')
    args = parser.parse_args()
    expected = json.loads((ROOT / 'model_manifest.json').read_text())['models'][args.model]
    destination = args.output or ROOT / 'checkpoints' / expected['filename']
    print('Verified', download(expected, destination))


if __name__ == '__main__':
    main()
