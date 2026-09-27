"""Build the versioned release asset from a validated local deploy_tcsp folder."""
import argparse
import hashlib
import json
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


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--source', type=Path, required=True, help='Local standalone deployment folder')
    parser.add_argument('--output', type=Path, default=ROOT/'release_artifacts/TCSPV3-v3.0.0-assets.tar.gz')
    args = parser.parse_args()
    for name, expected in EXPECTED.items():
        path = args.source/name
        if path.stat().st_size != expected['bytes'] or sha256(path) != expected['sha256']:
            raise ValueError(f'Asset mismatch: {path}')
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = args.output.with_name(args.output.name+'.partial')
    with tarfile.open(temporary, 'w:gz', compresslevel=3) as archive:
        for name in sorted(EXPECTED):
            archive.add(args.source/name, arcname=name, recursive=False)
    temporary.replace(args.output)
    release = {'version':'3.0.0','filename':args.output.name,'bytes':args.output.stat().st_size,
               'sha256':sha256(args.output),'asset_count':len(EXPECTED)}
    (ROOT/'release_manifest.json').write_text(json.dumps(release,indent=2)+'\n')
    print(json.dumps(release))


if __name__ == '__main__':
    main()
