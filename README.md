# TCSPV3

**TCSPV3** is a composition-conditioned, template-based crystal structure predictor. It extends [TCSP2.0](docs/TCSP2.md) with complementary template retrieval, multiple stoichiometry-preserving substitutions, soft oxidation-state evidence, size and coordination checks, structural deduplication, and optional MatterSim ranking. This repository contains the standalone V3 implementation; it does not replace the existing TCSP2.0 repository or connect to a CSP dispatcher.

**Developed by Jianjun Hu, Ying Feng, and Lai Wei at Machine Learning and Evolution Lab at University of South Carolina.

**Version:** `3.1.0` (source publication; the data/model bundle, version 3.0.0, is not publicly distributed here). Since 3.1.0 the default final stage is ORB-v3 relaxation and energy ranking; see [CHANGELOG.md](CHANGELOG.md). See [the separate V2 page](docs/TCSP2.md), [CHANGELOG.md](CHANGELOG.md), [algorithm details](docs/ALGORITHM.md), [benchmark protocol](docs/BENCHMARKS.md), and [data provenance](docs/DATA.md).

## What is new in TCSP3.0

TCSP2.0 already uses BERTOS oxidation-state predictions, element embeddings and periodic-group similarity, flexible substitution, and space-group voting ([paper](https://arxiv.org/abs/2503.23183); [V2 README](https://github.com/usccolumbia/TCSP/blob/main/README.md)). V3 changes how these signals are used:

| Area | V3 implementation |
|---|---|
| Retrieval | Unions embedding-based formulas with Villars periodic-number neighborhoods, inspired by [PNcsp+](https://doi.org/10.1021/acs.jctc.6c00044), preserving candidates from both routes under the search budget. |
| Element mapping | Keeps up to three near-optimal assignments that preserve reduced stoichiometric coefficients. |
| Oxidation states | Aggregates BERTOS probabilities per species and softly scores up to four neutral alternatives, rather than letting one uncertain assignment reject a template. |
| Geometry | Uses nearest-shell coordination and compatible Shannon/atomic radii for size penalties and optional bounded lattice scaling. |
| Candidate diversity | Deduplicates with `StructureMatcher` while retaining distinct polymorphs; majority-space-group voting is optional and was not part of the benchmarked default. |
| Final stage | By default relaxes the 20 native candidates with ORB-v3 and ranks them by ORB-v3 energy per atom (`--ranker orb-v3`); `--ranker native` keeps the chemical-score order without relaxation. An optional MatterSim single-point stage is also available; chemical and physical rankings remain separately inspectable. |
| Reproducibility | Checks actual source CIF composition, preserves source hashes, records failures, and supports versioned indexes and a shared read-only memory-mapped CIF archive. |

The periodic-number idea is credited to Oran, Caputo, Villars, and Tekin’s [PNcsp+ paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13085239/). TCSPV3 combines a PN retrieval route with its existing embedding route; it does not implement PNcsp+’s OQMD search, formation-enthalpy filter, or MACE/M3GNet/ALIGNN-FF ensemble.

## Performance and limits

All percentages below are **StructureMatcher** recovery, not stability rates. The exploratory 46-case internal comparison used a repaired soft-OS/dedup baseline, **not** an exact reproduction of TCSP2.0.

| Evaluation | Top-1 | Top-5 | Top-20 |
|---|---:|---:|---:|
| Testset 46 internal repaired baseline | 23/46 (50.0%) | 34/46 (73.9%) | 39/46 (84.8%) |
| Testset 46 V3 chemical ranking | 23/46 (50.0%) | 36/46 (78.3%) | 39/46 (84.8%) |
| Testset 46 V3 + radius-scaled MatterSim single-point | 28/46 (60.9%) | 38/46 (82.6%) | 39/46 (84.8%) |
| Testset 46 V3.1 default (ORB-v3 relax and rank, current T46 cohort) | **35/46 (76.1%)** | **41/46 (89.1%)** | **41/46 (89.1%)** |
| Original TCSP 180-case list, V3 chemical ranking | 115/180 (63.89%) | **153/180 (85.00%)** | not tabulated |
| Original TCSP 180-case list, V3 + radius-scaled MatterSim single-point | **116/180 (64.44%)** | 150/180 (83.33%) | not tabulated |
| T180 v3.1 (leakage-audited), V3.1 default (ORB-v3 relax and rank) | **130/180 (72.2%)** | **159/180 (88.3%)** | **163/180 (90.6%)** |
| Published TCSP2.0 on its 180-case protocol **(not controlled against V3)** | 68.3% | 78.33% | not reported |

The Top-20 for the single-point rows equals the chemical-ranking pool (the same 20 candidates are only reordered). The historical Testset 46 rows are the older exploratory cohort; the V3.1 row uses the current T46 cohort, so rows are not paired across cohorts.

The separate, current **T180 v3.1** benchmark recovers **130/180 (72.2%) at Top-1**, **159/180 (88.3%) at Top-5**, and **163/180 (90.6%) at Top-20** with the default ORB-v3 relaxation and energy ranking (native chemical ranking without relaxation: 122/154/163).

| Cohort (structure match, Top-1/5/20) | Native order, no relaxation | Default: ORB-v3 relax and rank |
|---|---:|---:|
| T46 (46 targets) | 34/39/41 | **35/41/41** |
| T180 v3.1 (180 targets) | 122/154/163 | **130/159/163** |
| T800 (800 targets; 19 without candidates count as misses) | 512/575/591 | **518/579/592** |

The ORB-v3 stage adds no candidates, so Top-20 stays (to within one target) the same; it improves the order. Without relaxation the template stage is fast: the median is 28 s per target on two CPU cores. The ORB-v3 stage adds about 219 s per target on a GPU, including loading the weights for each target. Download/browse the [ORB-v3 relaxed and ranked Top-20 prediction CIFs (default; 130/159/163), with per-target scores and hashes](benchmarks/T180_v3_1_orb_top20/README.md), or the [native-order Top-20 CIFs without relaxation (122/154/163)](benchmarks/T180_v3_1_top20/README.md). This package contains 2,581 available candidates across all 180 targets, with up to 20 per target. Its audited cohort and cell-count conditioning differ from the historical original-180 rows above; the package documents the actual input exception for `mp-557387`.



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
python tcspx.py --csv examples/formulas.csv --output results/parallel
python tcspx.py --csv examples/formulas.csv --workers 4 --formula-type reduced --output results/parallel4
python tcspx.py --csv your.csv --column FullFormula --workers 8 --output results/custom
```

The CSV reader finds an exact `composition` header (case-insensitive), or one header containing `formula`; use `--column` when ambiguous. Each row gets `row_<CSV row number>_<formula>` with a filename-safe formula. `batch_summary.csv` and `batch_summary.json` record successes and errors. Invalid rows do not stop the batch.

`tcspx.py` runs CPU workers, each reusing its loaded index/BERTOS model. By default it samples idle CPU capacity for one second, uses 60% of it, and caps the worker count by available memory; `--workers N` overrides this. Workers share the CIF archive through a read-only memory map, while parsed structures are cached per worker. A large batch can benefit from this; a speedup for short batches has not been established.

Both commands default to `--output ./cif_output`, relative to the current working directory. `--out` remains an alias. Use a new output directory for each run. The default index is `tcspv3/731K_index_minus_MP20_test.pkl` (719,450 entries). Use `--index tcspv3/731K_index.pkl` for the full 731,293-entry index or `--index tcspv3/731K_index_minus_MPST52_test.pkl` for the MPTS-52-test-excluded index. See [docs/DATA.md](docs/DATA.md) for exactly what each filtered index excludes.

## Use LeMat or your own CIFs as templates

See [the LeMat template guide](docs/LEMAT_TEMPLATES.md) for obtaining individual
CIFs, building a **parallel primitive-cell index and mmap CIF store**, recording
exclusions/provenance, and using the new library with `--index` and `--cif-store`.
The tested builder is `scripts/build_template_index.py`; auxiliary BERTOS and
embedding assets remain necessary. LeMat results are a separate experiment, not
the existing 731K benchmark scores.

## Optional MatterSim stage

Install `pip install -e '.[physics]'` in a compatible environment. MatterSim checkpoints are **downloaded from official upstream URLs, not bundled or rehosted here**.

[Matbench Discovery's MatterSim entry](https://github.com/janosh/matbench-discovery/blob/main/models/mattersim/mattersim-v1-5m.yml) lists [MatterSim v1.0.0-5M (direct checkpoint)](https://github.com/microsoft/mattersim/raw/refs/heads/main/pretrained_models/mattersim-v1.0.0-5M.pth). Our recorded benchmarks used **1M**, available from [the official 1M checkpoint URL](https://github.com/microsoft/mattersim/raw/refs/heads/main/pretrained_models/mattersim-v1.0.0-1M.pth). The models are distinct: choosing 5M does not reproduce the reported 1M results.

For the benchmark-compatible default, download 1M with size/SHA256 verification, then run the physical stage on a **single-formula** prediction directory:

```bash
python scripts/download_mattersim.py --model 1M
python -m tcspv3.energy --run cif_output --device cpu
```

To explicitly use the checkpoint linked by Matbench Discovery instead:

```bash
python scripts/download_mattersim.py --model 5M
python -m tcspv3.energy --run cif_output --device cpu --checkpoint checkpoints/mattersim-v1.0.0-5M.pth
```

Use a fresh prediction/output directory when changing the checkpoint. `model_manifest.json` records both official URLs, byte sizes and SHA256 hashes. The earlier locally prepared `TCSPV3-v3.0.0-assets.tar.gz` archive included MatterSim; it remains a historical artifact, not the current template-only bundle.

This separate stage computes raw/scaled single-point energies and limited relaxation; it can take substantially longer than chemical prediction. It does not turn a low within-composition energy into a stability claim. The parallel CSV CLI currently parallelizes chemical prediction only.

## Development, citation, and license

```bash
python -m pip install -e '.[test]'
python -m pytest tests -q
```

The source is MIT-licensed, following the [existing TCSP2.0 license](https://github.com/usccolumbia/TCSP/blob/main/LICENSE). Large data/model assets retain their own provenance and applicable terms. For attribution, cite this software version, the [TCSP2.0 paper](https://arxiv.org/abs/2503.23183), and [PNcsp+](https://doi.org/10.1021/acs.jctc.6c00044) for the periodic-number retrieval idea. See [CITATION.cff](CITATION.cff).

## BibTeX citation

```bibtex
@misc{tcspv3,
  author = {Feng, Ying and Wei, Lai and Hu, Jianjun},
  title = {TCSPV3},
  url = {https://github.com/usccolumbia/TCSPV3}
}
```
