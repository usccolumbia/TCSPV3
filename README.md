# TCSPV3

**TCSPV3** is a composition-conditioned, template-based crystal structure predictor. It extends [TCSP2.0](docs/TCSP2.md) with complementary template retrieval, multiple stoichiometry-preserving substitutions, soft oxidation-state evidence, size and coordination checks, structural deduplication, and optional MatterSim ranking. This repository contains the standalone V3 implementation; it does not replace the existing TCSP2.0 repository or connect to a CSP dispatcher.

**Developed by Jianjun Hu, Ying Feng, and Lai Wei at Machine Learning and Evolution Lab at University of South Carolina.

**Version:** `3.1.1`. The template data and model bundle (version 3.0.0) is on [Zenodo](https://doi.org/10.5281/zenodo.23627915). The default final stage is ORB-v3 relaxation and energy ranking. See [the separate V2 page](docs/TCSP2.md), [CHANGELOG.md](CHANGELOG.md), [algorithm details](docs/ALGORITHM.md), [benchmark protocol](docs/BENCHMARKS.md), and [data provenance](docs/DATA.md).

## What is new in TCSP3.0

TCSP2.0 already uses BERTOS oxidation-state predictions, element embeddings and periodic-group similarity, flexible substitution, and space-group voting ([paper](https://arxiv.org/abs/2503.23183); [V2 README](https://github.com/usccolumbia/TCSP/blob/main/README.md)). V3 changes how these signals are used:

| Area | V3 implementation |
|---|---|
| Retrieval | Ranks template formulas by the chemical similarity of their elements and keeps the best candidates under the search budget. |
| Element mapping | Keeps up to three near-optimal assignments that preserve reduced stoichiometric coefficients. |
| Oxidation states | Aggregates BERTOS probabilities per species and softly scores up to four neutral alternatives, rather than letting one uncertain assignment reject a template. |
| Geometry | Uses nearest-shell coordination and compatible Shannon/atomic radii for size penalties and optional bounded lattice scaling. |
| Candidate diversity | Deduplicates with `StructureMatcher` while retaining distinct polymorphs; majority-space-group voting is optional and was not part of the benchmarked default. |
| Final stage | By default relaxes the 20 native candidates with ORB-v3 and ranks them by ORB-v3 energy per atom (`--ranker orb-v3`); `--ranker native` keeps the chemical-score order without relaxation. An optional MatterSim single-point stage is also available; chemical and physical rankings remain separately inspectable. |
| Reproducibility | Checks actual source CIF composition, preserves source hashes, records failures, and supports versioned indexes and a shared read-only memory-mapped CIF archive. |

## Performance

Structure match (pymatgen StructureMatcher) with the reference, shown as matched targets (match rate) among the Top-1, Top-5 and Top-20 candidates; missing predictions count as misses.

| Test set | Top-k | TCSP 3.1 | CrystalFlow | MatterGen |
|---|---|---:|---:|---:|
| T46 (46) | Top-1 | **35 (76.1%)** | 26 (56.5%) | 25 (54.3%) |
|  | Top-5 | **41 (89.1%)** | 27 (58.7%) | 25 (54.3%) |
|  | Top-20 | **41 (89.1%)** | 31 (67.4%) | 27 (58.7%) |
| T180 (180) | Top-1 | **130 (72.2%)** | 115 (63.9%) | 97 (53.9%) |
|  | Top-5 | **159 (88.3%)** | 121 (67.2%) | 110 (61.1%) |
|  | Top-20 | **163 (90.6%)** | 130 (72.2%) | 119 (66.1%) |
| T240 (240) | Top-1 | **149 (62.1%)** | 130 (54.2%) | 103 (42.9%) |
|  | Top-5 | **186 (77.5%)** | 137 (57.1%) | 118 (49.2%) |
|  | Top-20 | **193 (80.4%)** | 149 (62.1%) | 127 (52.9%) |
| T800 (800) | Top-1 | **518 (64.8%)** | 70 (8.8%) | 43 (5.4%) |
|  | Top-5 | **579 (72.4%)** | 74 (9.2%) | 47 (5.9%) |
|  | Top-20 | **592 (74.0%)** | 78 (9.8%) | 59 (7.4%) |

All three methods relax and rank their candidates with ORB-v3. CrystalFlow and MatterGen use MP-20-trained checkpoints (structures up to 20 atoms); T800 targets have up to 52 atoms. Browse the [T180 prediction CIFs](benchmarks/T180_v3_1_orb_top20/README.md).

## Install

TCSPV3 needs Python 3.10–3.12 (the `dm-tree` dependency does not install on 3.13). Clone the repository and install it in editable mode:

```bash
git clone https://github.com/usccolumbia/TCSPV3.git
cd TCSPV3
python3.12 -m venv .env
source .env/bin/activate
python -m pip install -e .
```

If `python3.12` is not available, create the environment another way, then run `python -m pip install -e .` as above:

- With uv, which downloads Python 3.12 for you if it's missing:

  ```bash
  uv venv -p 3.12 --seed .env
  ```

- With conda:

  ```bash
  conda create -y -p .env python=3.12
  ```

  Then activate it with `conda activate ./.env` instead of `source .env/bin/activate`.

Download the template data and BERTOS model from Zenodo ([DOI 10.5281/zenodo.23627915](https://doi.org/10.5281/zenodo.23627915), 343 MB):

```bash
python scripts/install_data.py --download
```

Test the installation. The candidates are relaxed and ranked with ORB-v3 (`orb-models` is installed by the `pip` step above); the first run downloads the ORB-v3 weights.

Run with GPU (default; uses the least-busy GPU, or pick one with `--orb-device cuda:N`):

```bash
python tcsp.py --formula SrTiO3
```

Run without GPU (CPU only, slower):

```bash
python tcsp.py --formula SrTiO3 --orb-device cpu
```

Add `--ranker native` to skip the ORB-v3 stage and keep the chemical-score order.

## Predict one composition

```bash
python tcsp.py --formula SrTiO3                       # default: ORB-v3 relaxation and energy ranking of the 20 candidates
python tcsp.py --formula SrTiO3 --ranker native        # chemical-score order, no relaxation
python tcsp.py --formula Fe4O6 --formula-type full --output results/Fe4O6
python tcsp.py --formula Fe2O3 --formula-type reduced --output results/Fe2O3_multiple_Z
```

`--formula-type full` is the default: it requires the output cell's atom count to equal the count written in the formula. `--formula-type reduced` permits the template cell sizes available for that reduced composition and reports each candidate's atom count and formula-unit count. The latter does not promise primitive cells or arbitrary new values of Z. Historical V3 benchmarks used unrestricted template sizes, corresponding to `--formula-type reduced` in this CLI.

`--mode v3` is the default. The predictor excludes **all templates with the same reduced composition** as the query, and rechecks the actual source CIF. It creates `prediction.json`, `oxidation.json`, and ranked `tcspv3_###.cif` files (ORB-v3-relaxed structures in ORB-v3 energy order by default; the unrelaxed native candidates are kept in `native/`, and `prediction.json` records the final stage and each candidate's native rank). `--topk 20` is a maximum: fewer may be returned after structural deduplication or cell-size filtering. The default pool contains up to 100 valid instantiated candidates; increase `--pool-size` if more distinct structures are needed.

## Batch and parallel CSV prediction

```bash
python tcsp.py --csv examples/formulas.csv --output results/serial
python tcsp.py --csv examples/formulas.csv --workers auto --output results/parallel
python tcsp.py --csv your.csv --column FullFormula --workers 8 --output results/custom
```

The CSV reader finds an exact `composition` header (case-insensitive), or one header containing `formula`; use `--column` when ambiguous. Each row gets `row_<CSV row number>_<formula>` with a filename-safe formula. `batch_summary.csv` and `batch_summary.json` record successes and errors. Invalid rows do not stop the batch.

`--workers N` runs N CPU worker processes that build candidates in parallel, each loading the index and BERTOS once; `--workers auto` uses 60% of idle CPU capacity, capped by available memory. The ORB-v3 stage then relaxes each formula's candidates in the main process on the GPU. Results are identical to a serial run. `tcspx.py` is kept as an alias for `tcsp.py --workers auto`.

Runs default to `--output ./cif_output`, relative to the current working directory. `--out` remains an alias. Use a new output directory for each run. The default index is `tcspv3/731K_index_minus_MP20_test.pkl` (719,450 entries). Use `--index tcspv3/731K_index.pkl` for the full 731,293-entry index or `--index tcspv3/731K_index_minus_MPST52_test.pkl` for the MPTS-52-test-excluded index. See [docs/DATA.md](docs/DATA.md) for exactly what each filtered index excludes.

## Use LeMat or your own CIFs as templates

See [the LeMat template guide](docs/LEMAT_TEMPLATES.md) for obtaining individual
CIFs, building a **parallel primitive-cell index and mmap CIF store**, recording
exclusions/provenance, and using the new library with `--index` and `--cif-store`.
The tested builder is `scripts/build_template_index.py`; auxiliary BERTOS and
embedding assets remain necessary. LeMat results are a separate experiment, not
the existing 731K benchmark scores.

## Development, citation, and license

```bash
python -m pip install -e '.[test]'
python -m pytest tests -q
```

The source is MIT-licensed, following the [existing TCSP2.0 license](https://github.com/usccolumbia/TCSP/blob/main/LICENSE). Large data/model assets retain their own provenance and applicable terms. For attribution, cite this software version and the [TCSP2.0 paper](https://arxiv.org/abs/2503.23183). See [CITATION.cff](CITATION.cff).

## BibTeX citation

```bibtex
@misc{tcspv3,
  author = {Feng, Ying and Wei, Lai and Hu, Jianjun},
  title = {TCSPV3},
  url = {https://github.com/usccolumbia/TCSPV3}
}
```
