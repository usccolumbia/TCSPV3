# T180 v3.1: TCSPV3 ranked predictions through Top-20

| Metric | Recovered targets | Rate |
|---|---:|---:|
| Top-1 | 122/180 | 67.8% |
| Top-5 | 154/180 | 85.6% |
| Top-20 | **163/180** | **90.6%** |

This package contains **2,581 available ranked CIFs for all 180 targets** in the leakage-audited T180 release v3.1. Each `predictions/<material-id>/` directory preserves the native chemical ranking through rank 20. Seventy-nine targets have 20 candidates; the other 101 have fewer distinct eligible template candidates. No candidates were duplicated or padded. All 180 targets remain in the metric denominator.

The results use formula-only, integer cell-count conditioning, a 100-candidate instantiation budget, native TCSPV3 chemical ranking, and no physical relaxation. Recovery uses pymatgen StructureMatcher with `ltol=0.2`, `stol=0.3`, `angle_tol=5`, `primitive_cell=True`, `scale=True`, and `attempt_supercell=False`. These are structural-recovery metrics, not stability or experimental-validation claims.

T180 v3.1 replaces `mp-1183076` with `mp-3576`. It is a different cohort and protocol from the repository's historical **original TCSP 180-case** table. Those historical results remain separately documented.

The actual query for `mp-557387` was **Sr4Tm2Nb2O12 (20 atoms)**. The campaign explicitly resolved a 0.01-A pseudo-translation ambiguity using that cell, whereas the release metadata lists Nb1O6Sr2Tm1 (10 atoms). `queries.json` records both the release formula and actual input; the predictions and scores are preserved without changing the historical input. Other targets use their recorded campaign formulas, including the replacement target.

- `per_target.csv` / `per_target.json`: actual input formula, reference formula, pool size, first match rank, hit flags and reference CIF SHA256.
- `file_manifest.json`: every prediction's rank, relative path, byte size and SHA256.
- `queries.json`: release queries plus actual campaign input formulas.
- `summary.json`: counts, protocol and the explicit input exception.
- `verification.json`: independent rematching and verification of all exported CIF hashes and compositions.

Reference IDs identify Materials Project structures. The reference CIF bytes are identified by SHA256 and are not bundled here. Prediction CIFs are computed benchmark artifacts; the software license does not assert new licensing terms for upstream template collections. The full template database and model weights are not included.
