# Changelog

## Source updates — 2026-10-09

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
