"""Install a locally obtained versioned release data bundle after SHA256 checks."""
import argparse
import hashlib
import json
import os
import tarfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = json.loads((ROOT/'asset_manifest.json').read_text())['assets']


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(4*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def install(bundle, check_only=False):
    bundle = Path(bundle)
    release = json.loads((ROOT/'release_manifest.json').read_text())
    if bundle.stat().st_size != release['bytes'] or sha256(bundle) != release['sha256']:
        raise ValueError('Release bundle size/SHA256 mismatch')
    seen = set()
    with tarfile.open(bundle, 'r:gz') as archive:
        members = archive.getmembers()
        for member in members:
            name = member.name
            if name not in EXPECTED or not member.isfile() or name in seen:
                raise ValueError(f'Unexpected, duplicate, or unsafe bundle member: {name}')
            seen.add(name)
        if seen != set(EXPECTED):
            raise ValueError(f'Bundle inventory mismatch: missing {sorted(set(EXPECTED)-seen)}')
        for member in members:
            target = ROOT/member.name
            expected = EXPECTED[member.name]
            if target.exists():
                if target.stat().st_size != expected['bytes'] or sha256(target) != expected['sha256']:
                    raise ValueError(f'Existing file differs from release: {target}')
                continue
            if check_only:
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            temporary = target.with_name(target.name+'.partial')
            digest = hashlib.sha256()
            with archive.extractfile(member) as source, temporary.open('wb') as output:
                for chunk in iter(lambda: source.read(4*1024*1024), b''):
                    digest.update(chunk)
                    output.write(chunk)
            if temporary.stat().st_size != expected['bytes'] or digest.hexdigest() != expected['sha256']:
                temporary.unlink(missing_ok=True)
                raise ValueError(f'Asset size/SHA256 mismatch: {member.name}')
            os.replace(temporary, target)
    return len(seen)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--bundle', required=True, help='Path to TCSPV3-v3.0.0-assets.tar.gz')
    parser.add_argument('--check-only', action='store_true', help='Verify bundle inventory and checksum without installing')
    args = parser.parse_args()
    print(f'Verified {install(args.bundle, args.check_only)} release assets')


if __name__ == '__main__':
    main()
