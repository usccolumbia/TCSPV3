# TCSPV3

**TCSPV3** is a composition-conditioned, template-based crystal structure predictor. It extends [TCSP2.0](https://github.com/usccolumbia/TCSP) with complementary template retrieval, multiple stoichiometry-preserving substitutions, soft oxidation-state evidence, size and coordination checks, structural deduplication, and optional MatterSim ranking. This repository contains the standalone V3 implementation; it does not replace the existing TCSP2.0 repository or connect to a CSP dispatcher.

**Version:** `3.0.0` (release candidate prepared locally; no public V3 release is claimed). See [CHANGELOG.md](CHANGELOG.md), [algorithm details](docs/ALGORITHM.md), [benchmark protocol](docs/BENCHMARKS.md), and [data provenance](docs/DATA.md).

## What changed from TCSP2.0

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

Published TCSP2.0 Figure 1 reports **68.3% Top-1**, and its paper abstract reports **78.33% Top-5** on its 180 cases. Its Top-1 is above V3's best tested **64.44% Top-1**, while the reported Top-5 is below V3's **85.00% chemical Top-5**. The protocols differ in relaxation, accepted reference polymorphs, and template exclusion, so these numbers do **not** establish that either implementation is superior under matched conditions. In V3's original-180 run, MatterSim gave a net gain of one Top-1 recovery but a net loss of three Top-5 recoveries. The 46-case set had been used during method inspection. Full denominators, tolerances, resource settings, frozen upstream commit, and caveats are in [docs/BENCHMARKS.md](docs/BENCHMARKS.md).

## Install

Python 3.10 or newer is recommended. Clone this source repository and install it in editable mode:

```bash
git clone <this-repository-url> TCSPV3
cd TCSPV3
python3 -m venv .venv
source .venv/bin/activate
python -m pip install -e .
```

The **versioned data bundle** is separate from Git history because the CIF archive alone is 1.35 GB. Download the asset named `TCSPV3-v3.0.0-assets.tar.gz` from the corresponding V3 release (or use the locally prepared bundle), then install it with checksum verification:

```bash
python scripts/install_data.py --bundle /path/to/TCSPV3-v3.0.0-assets.tar.gz
```

The installer checks the release SHA256, file inventory, and each asset SHA256. Until the asset bundle is published, this source checkout requires the locally prepared bundle. The code and dataset/model licensing are distinct; review [docs/DATA.md](docs/DATA.md) before redistributing the assets.

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

## Optional MatterSim stage

Install `pip install -e '.[physics]'` in a compatible environment. The data bundle includes the benchmark MatterSim 1M checkpoint. For a **single-formula** prediction directory, run:

```bash
python -m tcspv3.energy --run cif_output --device cpu
```

This separate stage computes raw/scaled single-point energies and limited relaxation; it can take substantially longer than chemical prediction. It does not turn a low within-composition energy into a stability claim. The parallel CSV CLI currently parallelizes chemical prediction only.

## Development, citation, and license

```bash
python -m pip install -e '.[test]'
python -m pytest tests -q
```

The source is MIT-licensed, following the [existing TCSP2.0 license](https://github.com/usccolumbia/TCSP/blob/main/LICENSE). Large data/model assets retain their own provenance and applicable terms. For attribution, cite this software version, the [TCSP2.0 paper](https://arxiv.org/abs/2503.23183), and [PNcsp+](https://doi.org/10.1021/acs.jctc.6c00044) for the periodic-number retrieval idea. See [CITATION.cff](CITATION.cff).
