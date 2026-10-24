# Zenodo release of the data bundle

Published 2026-10-24: record [23627915](https://zenodo.org/records/23627915), DOI [10.5281/zenodo.23627915](https://doi.org/10.5281/zenodo.23627915), MIT licence.

- File: `TCSPV3-v3.0.0-template-assets.tar.gz`, 342,903,364 bytes
- SHA256: `afe73e7940f1a6037313031e303069481bcfdb096230783a8db437f9d2a66157`
- `metadata.json`: the record metadata as uploaded

## New version

1. Build: `python scripts/build_release_bundle.py --source <deployed assets>`
2. On the Zenodo record, click "New version", upload the file, publish.
3. Put the new URL and DOI in `release_manifest.json`, `README.md` and `docs/DATA.md`.
4. Test: `python scripts/install_data.py --download --check-only`
