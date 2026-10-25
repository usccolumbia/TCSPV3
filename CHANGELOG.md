# Changelog

## 3.1.1 — 2026-10-24

- Template data and BERTOS model download from Zenodo: `python scripts/install_data.py --download` ([10.5281/zenodo.23627915](https://doi.org/10.5281/zenodo.23627915)).
- Parallel CSV runs in `tcsp.py` with `--workers N` or `--workers auto`; serial and parallel runs share one pipeline and give identical candidates. `tcspx.py` is kept as an alias.
- `--orb-device auto` uses the least-busy GPU; `cuda:N` picks one, `cpu` runs without a GPU.
- Cleaner output: version banner (`--version`, `--no-banner`), stage messages and an ORB-v3 progress bar; library warnings hidden unless `--show-warnings`.
- An existing single-formula prediction is overwritten with a warning instead of an error; `--help` shows defaults.
- Python 3.10–3.12 (`dm-tree` does not install on 3.13); adds `nvidia-ml-py` and `tqdm`.
- Shorter README with one performance table against CrystalFlow and MatterGen.
- No change to retrieval, scoring or the ORB-v3 settings; benchmark results are unchanged.

## 3.1.0 — 2026-10-14

- **ORB-v3 relaxation and energy ranking is now the default final stage** (`--ranker orb-v3`). The 20 native candidates are relaxed with ORB-v3 (`orb-v3-conservative-inf-omat`; FIRE with ASE `ExpCellFilter`, at most 500 steps, force threshold 0.01 eV/Å, nonconverged structures kept) and ranked by ORB-v3 energy per atom, the same post-processing used for the deep generative models in the benchmark. The native candidates are saved in `native/`; `--ranker native` keeps the chemical-score order without relaxation. New options: `--orb-device`, `--orb-weights`, `--orb-steps`, `--orb-fmax`.
- Benchmark effect (structure match, Top-1/5/20): T46 34/39/41 → 35/41/41, T180 122/154/163 → 130/159/163, T800 512/575/591 → 518/579/592; Top-20 is unchanged to within one target because the stage adds no candidates. The native results in the table below remain valid for `--ranker native`.
- `ase` and `orb-models>=0.5.5` are now core dependencies. No template library, retrieval, mapping or scoring code changed; the data bundle is still 3.0.0.

## Source updates — 2026-10-09

- Added the current T180 v3.1 Top-20 result package: 2,581 ranked CIFs for 180 targets, per-target metrics and checksums, independent matcher verification, and the documented actual-query exception. Native chemical recovery is 122/154/163 at Top-1/5/20; the historical original-180 results are unchanged.
- Added a LeMat template guide, revision-pinned streaming CIF exporter, and parallel primitive-cell native index/mmap builder with explicit provenance and rejection records.
- Replaced bundled MatterSim weights with official upstream download links and size/SHA256 verification; distinguish the Matbench-listed 5M checkpoint from the historical 1M benchmark model.
- Added the repository URL BibTeX citation with Ying Feng, Lai Wei, and Jianjun Hu as authors.
- No full LeMat template-library build or new benchmark evaluation was performed for these source updates.

## 3.0.0 — prepared 2026-09-27

- Added hybrid embedding/periodic-number retrieval, multiple coefficient-preserving assignments, soft BERTOS neutral-state evidence, coordination-aware radius penalties, and structural deduplication.
- Added optional radius scaling, MatterSim single-point ranking, and limited relaxation.
- Added single-formula, serial CSV, and parallel CSV command-line tools; full/reduced formula modes; three selectable template indexes.
- Added a shared read-only memory-mapped CIF archive and checksum-audited versioned data bundle.
- Documented Testset 46 and exact upstream TCSP 180-case evaluations, including protocol limits relative to published TCSP2.0.

This version records the standalone implementation and does not change the existing TCSP2.0 repository or CSP dispatcher.
