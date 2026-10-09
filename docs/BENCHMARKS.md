# Benchmark evidence and comparison to TCSP2.0

**Metric.** Top-1/Top-5/Top-20 mean that at least one of the first 1/5/20 ranked candidates matches a reference with pymatgen StructureMatcher (`ltol=0.2`, `stol=0.3`, angle tolerance 5°, primitive-cell reduction and scaling enabled, no supercell matching). All cases remain in the denominator. These are structure-recovery rates, not stability or experimental validation.

## Current leakage-audited T180 v3.1

| Ranking | Top-1 | Top-5 | Top-20 |
|---|---:|---:|---:|
| TCSPV3 native chemical | 122/180 (67.8%) | 154/180 (85.6%) | **163/180 (90.6%)** |

The [complete result package](../benchmarks/T180_v3_1_top20/README.md) contains 2,581 ranked prediction CIFs, all 180 per-target records, checksums and independent verification. It uses integer cell-count inputs and no physical relaxation. Seventy-nine pools contain 20 candidates; the others contain fewer available native templates. The denominator remains 180. The package documents the 20-atom actual query for the `mp-557387` primitive-cell ambiguity and the `mp-1183076` to `mp-3576` reference replacement. This release is separate from the historical original TCSP 180-case cohort below.

## Exploratory Testset 46

| Method | Top-1 | Top-5 |
|---|---:|---:|
| Repaired soft-OS + dedup internal baseline | 23/46 (50.0%) | 34/46 (73.9%) |
| PN-only retrieval | 24/46 (52.2%) | 33/46 (71.7%) |
| Embedding + PN retrieval | 23/46 (50.0%) | 35/46 (76.1%) |
| TCSPV3 chemical ranking | 23/46 (50.0%) | 36/46 (78.3%) |
| TCSPV3 + MatterSim single-point | 26/46 (56.5%) | 36/46 (78.3%) |
| TCSPV3 + radius scaling + single-point | **28/46 (60.9%)** | **38/46 (82.6%)** |

The full V3 pipeline gained 5 Top-1 and 4 Top-5 recoveries relative to the **internal repaired baseline**. That baseline is not a byte-for-byte TCSP2.0 reproduction. The 46 cases had already been used for method inspection, so these are exploratory ablations rather than untouched holdout evidence. Limited relaxation did not improve Top-5 recovery in this experiment. See [report](benchmark_sources/Testset46_REPORT.md) and [protocol](benchmark_sources/Testset46_PROTOCOL.md).

## Exact original TCSP 180-case cohort

The test list is the upstream [`data/180_testdata.csv`](https://github.com/usccolumbia/TCSP/blob/main/data/180_testdata.csv), pinned at commit `4da9ee6eae151732b406cec7df81e16e9d0a7125`. It is **not** the earlier curated local 180-case set, which shares only 31 IDs. References came from the archived original-paper ground-truth ZIP and were verified against all 180 IDs/compositions/space groups. This run used CPU partition, two CPUs per task, 100 instantiated candidates, up to 20 distinct candidates, one reference CIF per target, and no physical relaxation.

| Ranking | Top-1 | Top-5 |
|---|---:|---:|
| TCSPV3 chemical | 115/180 (63.89%) | **153/180 (85.00%)** |
| TCSPV3 radius-scaled MatterSim single-point | **116/180 (64.44%)** | 150/180 (83.33%) |

The physical rerank gains 30 and loses 29 Top-1 cases (net +1), and gains 6 but loses 9 Top-5 cases (net -3). Zero energy-model errors were recorded; 21 cases had fewer than five predictions, with none empty. See [report](benchmark_sources/Original180_REPORT.md), [protocol](benchmark_sources/Original180_PROTOCOL.md), and [input manifest](benchmark_sources/Original180_input_manifest.json).

## Comparison with published TCSP2.0

The [TCSP2.0 paper](https://arxiv.org/abs/2503.23183) reports **68.3% Top-1** in Figure 1 and **78.33% Top-5** structural similarity in its abstract on its 180-case benchmark. V3's best tested Top-1 on the exact upstream list is **64.44%**, 3.86 percentage points lower numerically. This is **not a controlled head-to-head**: the paper relaxes predictions, accepts multiple reference polymorphs, and describes exact-entry template exclusion, whereas this V3 run uses single-point ranking, one archived reference per target, and excludes every target composition from templates. The repository therefore does **not** claim that V3 outperforms published V2.

To establish a V2/V3 improvement, run both versions with the same template pool and exclusions, reference-polymorph set, relaxation budget, candidate budget, and matcher settings. Also evaluate the updated `--formula-type full` default separately: these historical V3 benchmarks used unrestricted template cell sizes, now specified by `--formula-type reduced`. The historical run used the full 731,293-entry index with query-composition exclusion; the current CLI defaults to the MP20-test-excluded index, so pass `--index tcspv3/731K_index.pkl` for a closer input configuration.
