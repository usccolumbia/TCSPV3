# TCSPV3

**TCSPV3** is a composition-conditioned, template-based crystal structure predictor. It extends [TCSP2.0](docs/TCSP2.md) with complementary template retrieval, multiple stoichiometry-preserving substitutions, soft oxidation-state evidence, size and coordination checks, structural deduplication, and optional MatterSim ranking. This repository contains the standalone V3 implementation; it does not replace the existing TCSP2.0 repository or connect to a CSP dispatcher.

**Developed by Jianjun Hu, Ying Feng, and Lai Wei at Machine Learning and Evolution Lab at University of South Carolina.

**Version:** `3.0.0` (source publication; the data/model bundle is not publicly distributed here). See [the separate V2 page](docs/TCSP2.md), [CHANGELOG.md](CHANGELOG.md), [algorithm details](docs/ALGORITHM.md), [benchmark protocol](docs/BENCHMARKS.md), and [data provenance](docs/DATA.md).

## What is new in TCSP3.0

TCSP2.0 already uses BERTOS oxidation-state predictions, element embeddings and periodic-group similarity, flexible substitution, and space-group voting ([paper](https://arxiv.org/abs/2503.23183); [V2 README](https://github.com/usccolumbia/TCSP/blob/main/README.md)). V3 changes how these signals are used:

| Area | V3 implementation |
|---|---|
| Retrieval | Unions embedding-based formulas with Villars periodic-number neighborhoods, inspired by [PNcsp+](https://doi.org/10.1021/acs.jctc.6c00044), preserving candidates from both routes under the search budget. |
| Element mapping | Keeps up to three near-optimal assignments that preserve reduced stoichiometric coefficients. |
| Oxidation states | Aggregates BERTOS probabilities per species and softly scores up to four neutral alternatives, rather than letting one uncertain assignment reject a template. |
| Geometry | Uses nearest-shell coordination and compatible Shannon/atomic radii for size penalties and optional bounded lattice scaling. |
| Candidate diversity | Deduplicates with `StructureMatcher` while retaining distinct polymorphs; majority-space-group voting is optional and was not part of the benchmarked default. |
| Physical stage | Optionally reranks the same candidates with MatterSim single-point energies and limited relaxation; chemical ranking and physical ranking remain separately inspectable. |
| Reproducibility | Checks actual source CIF composition, preserves source hashes, records failures, and supports versioned indexes and a shared read-only memory-mapped CIF archive. |

The periodic-number idea is credited to Oran, Caputo, Villars, and Tekin’s [PNcsp+ paper](https://pmc.ncbi.nlm.nih.gov/articles/PMC13085239/). TCSPV3 combines a PN retrieval route with its existing embedding route; it does not implement PNcsp+’s OQMD search, formation-enthalpy filter, or MACE/M3GNet/ALIGNN-FF ensemble.

## Performance and limits

All percentages below are **StructureMatcher** recovery, not stability rates. The exploratory 46-case internal comparison used a repaired soft-OS/dedup baseline, **not** an exact reproduction of TCSP2.0.

| Evaluation | Top-1 | Top-5 |
|---|---:|---:|
| Testset 46 internal repaired baseline | 23/46 (50.0%) | 34/46 (73.9%) |
| Testset 46 V3 chemical ranking | 23/46 (50.0%) | 36/46 (78.3%) |
| Testset 46 V3 + radius-scaled MatterSim single-point | **28/46 (60.9%)** | **38/46 (82.6%)** |
| Original TCSP 180-case list, V3 chemical ranking | 115/180 (63.89%) | **153/180 (85.00%)** |
| Original TCSP 180-case list, V3 + radius-scaled MatterSim single-point | **116/180 (64.44%)** | 150/180 (83.33%) |
| Published TCSP2.0 on its 180-case protocol **(not controlled against V3)** | 68.3% | 78.33% |

The separate, current **T180 v3.1** benchmark recovers **122/180 (67.8%) at Top-1**, **154/180 (85.6%) at Top-5**, and **163/180 (90.6%) at Top-20** by native chemical ranking. Download/browse the [ranked prediction CIFs, per-target scores and hashes](benchmarks/T180_v3_1_top20/README.md). This package contains 2,581 available candidates across all 180 targets, with up to 20 per target. Its audited cohort and cell-count conditioning differ from the historical original-180 rows above; the package documents the actual input exception for `mp-557387`.



## Install

Python 3.10 or newer is recommended. Clone this source repository and install it in editable mode:

```bash
git clone https://github.com/usccolumbia/TCSPV3.git
cd TCSPV3
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The **versioned template/BERTOS data bundle** is separate from Git history because the CIF archive alone is 1.35 GB. **This source publication does not include a public template/BERTOS bundle download.** If you already have the locally prepared, authorized bundle named `TCSPV3-v3.0.0-template-assets.tar.gz`, install it with checksum verification:

```bash
python scripts/install_data.py --bundle /path/to/TCSPV3-v3.0.0-template-assets.tar.gz
```

The installer verifies the bundle SHA256 and file inventory, and checks installed/extracted asset SHA256 hashes. Source installation and unit tests work without the bundle; actual prediction requires these assets. Public distribution of the bundle is pending review of upstream dataset/model terms. The code and dataset/model licensing are distinct; review [docs/DATA.md](docs/DATA.md) before redistributing the assets.

## Predict one composition

```bash
python tcsp.py --formula SrTiO3
python tcsp.py --formula Fe4O6 --formula-type full --output results/Fe4O6
python tcsp.py --formula Fe2O3 --formula-type reduced --output results/Fe2O3_multiple_Z
```

`--formula-type full` is the default: it requires the output cell's atom count to equal the count written in the formula. `--formula-type reduced` permits the template cell sizes available for that reduced composition and reports each candidate's atom count and formula-unit count. The latter does not promise primitive cells or arbitrary new values of Z. Historical V3 benchmarks used unrestricted template sizes, corresponding to `--formula-type reduced` in this CLI.

`--mode v3` is the default. The predictor excludes **all templates with the same reduced composition** as the query, and rechecks the actual source CIF. It creates `prediction.json`, `oxidation.json`, and ranked `tcspv3_###.cif` files. `--topk 20` is a maximum: fewer may be returned after structural deduplication or cell-size filtering. The default pool contains up to 100 valid instantiated candidates; increase `--pool-size` if more distinct structures are needed.

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
