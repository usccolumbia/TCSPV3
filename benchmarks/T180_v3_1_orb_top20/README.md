# T180 v3.1: TCSPV3 default (ORB-v3 relaxed and ranked) predictions through Top-20

| Metric | Recovered targets | Rate |
|---|---:|---:|
| Top-1 | 130/180 | 72.2% |
| Top-5 | 159/180 | 88.3% |
| Top-20 | **163/180** | **90.6%** |

This is the output of the default TCSPV3 configuration (release 3.1.0): the same candidate pools as the [native-order package](../T180_v3_1_top20/README.md), each candidate relaxed with ORB-v3 and ranked by ORB-v3 energy per atom. It holds **2,581 relaxed, ranked CIFs for all 180 targets** (up to 20 each; 79 targets have 20 candidates, the others have fewer distinct eligible template candidates; nothing is padded). All 180 targets stay in the denominator.

Relaxation: `orb-v3-conservative-inf-omat`, FIRE with ExpCellFilter, force threshold 0.01 eV/A, at most 500 steps; candidates that did not converge are kept and ranked by their final energy. Query, candidate budget and matcher are as in the native package (formula-only integer cell-count input, 100 instantiated candidates, pymatgen `StructureMatcher` with `ltol=0.2`, `stol=0.3`, `angle_tol=5`, `primitive_cell=True`, `scale=True`, `attempt_supercell=False`). These are structural-recovery metrics, not stability or experimental-validation claims. The relaxation changes the order, not the pool: native Top-1/5/20 is 122/154/163.

File names are `rank<NN>_<formula>_<candidate id>.cif`, where the rank is the ORB-v3 energy rank and the candidate id (`tcspv3_NNNN`) is the native chemical-score position. Use `file_manifest.json` for paths rather than rebuilding names.

- `file_manifest.json`: every CIF with ORB-v3 rank, candidate id, native rank, formula, ORB-v3 energy (eV/atom), convergence flag, size and SHA256.
- `per_target.csv` / `per_target.json`: input and reference formula, pool size, first-match rank under ORB-v3 ranking and under native ranking, hit flags, reference CIF SHA256.
- `queries.json`, `summary.json`: queries (including the `mp-557387` input exception described in the native package) and counts, protocol and relaxation settings.

Reference IDs identify Materials Project structures; reference CIF bytes are identified by SHA256 and are not bundled. The software license does not assert new licensing terms for upstream template collections.
