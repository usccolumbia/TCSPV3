"""Download or install the TCSPV3 data bundle, checking SHA256."""
import argparse
import hashlib
import json
import os
import tarfile
import tempfile
from pathlib import Path
from urllib.request import Request, urlopen

ROOT = Path(__file__).resolve().parents[1]
EXPECTED = json.loads((ROOT/'asset_manifest.json').read_text())['assets']
RELEASE = json.loads((ROOT/'release_manifest.json').read_text())


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(4*1024*1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def download(release, destination, url=None):
    """Download the bundle; keep it only if size and SHA256 match."""
    destination = Path(destination)
    if destination.exists():
        if destination.stat().st_size != release['bytes'] or sha256(destination) != release['sha256']:
            raise ValueError(f'Existing file differs from release bundle: {destination}')
        return destination  # already downloaded
    url = url or release.get('download_url')
    if not url:
        raise ValueError('release_manifest.json has no download_url yet; pass --url or use --bundle with a local copy')
    destination.parent.mkdir(parents=True, exist_ok=True)
    request = Request(url, headers={'User-Agent': 'TCSPV3-data-installer'})
    handle, temporary = tempfile.mkstemp(prefix=destination.name+'.', suffix='.partial', dir=destination.parent)
    temporary = Path(temporary)
    try:
        digest, size = hashlib.sha256(), 0
        with os.fdopen(handle, 'wb') as output, urlopen(request, timeout=60) as response:
            for chunk in iter(lambda: response.read(4*1024*1024), b''):
                digest.update(chunk)
                size += len(chunk)
                output.write(chunk)
        if size != release['bytes'] or digest.hexdigest() != release['sha256']:
            raise ValueError(f'Downloaded bundle size/SHA256 mismatch from {url}')
        os.replace(temporary, destination)
    finally:
        temporary.unlink(missing_ok=True)
    return destination


def install(bundle, check_only=False):
    bundle = Path(bundle)
    release = RELEASE
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
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument('--bundle', help=f'Path to a local {RELEASE["filename"]}')
    source.add_argument('--download', action='store_true', help='Download the bundle from the archived release (Zenodo) first')
    parser.add_argument('--url', help='Override the download URL recorded in release_manifest.json')
    parser.add_argument('--cache-dir', type=Path, default=ROOT/'release_artifacts', help='Where --download stores the bundle')
    parser.add_argument('--check-only', action='store_true', help='Verify bundle inventory and checksum without installing')
    args = parser.parse_args()
    bundle = download(RELEASE, args.cache_dir/RELEASE['filename'], args.url) if args.download else args.bundle
    print(f'Verified {install(bundle, args.check_only)} release assets')


if __name__ == '__main__':
    main()
