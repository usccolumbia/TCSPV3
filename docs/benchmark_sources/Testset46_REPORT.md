# TCSPV3 — Testset 46 results

Standalone source: `/scratch/projects/apps/tcspv3`. Run: `/scratch/projects/jobs/csp/tcspv3_test46_20260927`. The CSP dispatcher was not modified.

All 46 cases were evaluated. This is an exploratory benchmark on a set already used for method inspection, not an untouched validation result. All 46 reduced compositions were excluded from the template pool and the actual source CIFs were audited.

## Recovery

| Method | SM@1 | SM@5 | SM@20 | SG@5 | Pool oracle |
|---|---:|---:|---:|---:|---:|
| Previous soft-OS + dedup baseline | 23/46 (50.0%) | 34/46 (73.9%) | 39/46 (84.8%) | 38/46 (82.6%) | 41/46 (89.1%) |
| PN-only retrieval | 24/46 (52.2%) | 33/46 (71.7%) | 39/46 (84.8%) | 37/46 (80.4%) | 39/46 (84.8%) |
| Embedding + PN retrieval | 23/46 (50.0%) | 35/46 (76.1%) | 41/46 (89.1%) | 37/46 (80.4%) | 42/46 (91.3%) |
| TCSPV3 chemical ranking | 23/46 (50.0%) | 36/46 (78.3%) | 39/46 (84.8%) | 34/46 (73.9%) | 42/46 (91.3%) |
| TCSPV3 + MatterSim single-point | 26/46 (56.5%) | 36/46 (78.3%) | 39/46 (84.8%) | 36/46 (78.3%) | — |
| TCSPV3 + radius scaling + single-point | 28/46 (60.9%) | 38/46 (82.6%) | 39/46 (84.8%) | 37/46 (80.4%) | — |
| Chemical top 5, limited relaxation | 24/46 (52.2%) | 36/46 (78.3%) | — | 35/46 (76.1%) | — |
| Single-point top 5, limited relaxation | 24/46 (52.2%) | 36/46 (78.3%) | — | 37/46 (80.4%) | — |

SM is pymatgen StructureMatcher recovery with ltol=0.2, stol=0.3, angle_tol=5°, primitive_cell=True, scale=True, attempt_supercell=False. SG uses symprec=0.01 Å and angle_tolerance=5°. All denominators remain 46, including failures and empty searches.

The baseline is the repaired composition-only retrieval baseline, not a byte-for-byte reproduction of TCSP 2.0. The earlier verified SM@5 was 34/46; the new baseline result above is checked against it. Every retrieval arm has a 100-valid-structure cap, but PN adds a formula route, so total retrieval costs are not identical. The oracle examines the candidate pool only during evaluation and does not guide ranking.

Single-point methods rerank the same top-20 distinct TCSPV3 candidates. Scaling uses only query/template radii. The two relaxed arms each spend five candidate relaxations, reusing shared candidates; 100 FIRE steps, fmax=0.05 eV/Å, variable cell, no symmetry constraint. Relaxed SM counts may include unconverged optimizations, so convergence-qualified recovery is reported below.

## Relaxation and failures

| Arm | Converged candidates | Failed candidates | Recovery with a converged candidate at top 5 |
|---|---:|---:|---:|
| Chemical top 5, limited relaxation | 132/222 | 0 | 35/46 (76.1%) |
| Single-point top 5, limited relaxation | 139/222 | 0 | 34/46 (73.9%) |

Unique relaxation tasks: 294. Shared candidates are reused between relaxed arms. Model identity: MatterSim v1.0.0 1M; checkpoint SHA256 `28b0b0b0f13efefee06b47ea4c9105a26bd3e2c8396da193430da96b3b49a8be`.

Energy ranks are within one composition. Low energy is not a stability claim. Model failures retain explicit error records; missing values are not zero. No MACE/M3GNet ensemble claim is made.

Errors by output-arm occurrence: `{}`.

## Paired top-5 changes

| Comparison | Gained cases | Lost cases |
|---|---:|---:|
| baseline -> pn | 1 | 2 |
| baseline -> union | 1 | 0 |
| baseline -> v3 | 2 | 0 |
| v3 -> singlepoint | 1 | 1 |
| singlepoint -> scaled_singlepoint | 2 | 0 |
| v3 -> v3_relaxed | 0 | 0 |
| singlepoint -> singlepoint_relaxed | 0 | 0 |

Exact gained/lost material IDs are in `results/summary.json`. No weights, cutoffs, or match tolerances were selected using these outcomes.

## Audit and artifacts

- 15408 candidate instances validated against the requested composition.
- 5417 unique source CIFs checked for actual composition exclusions and SHA256 identity.
- 731293 template-index rows; index SHA256 `874ae374775b2666ccf6b0e1d746b2d96b586a797c50b153f2ce57df51b8aee6`.
- BERTOS neutral alternatives available for 32/46 queries. All-metal targets intentionally omit ionic penalties.
- Testset SHA256 `dc316c430f614856f05db6987b58f72ad42702544066138fb3395ad4c5fc24be`; frozen protocol SHA256 `95ae40e053cfb71dd7d9bb28cb36e274ba8ea683868f8a72a630e100f512f24a`.
- Generation wall time: 1487.0 seconds; scheduler records and per-case timings provide the physical-stage timing.
- Six unit tests passed before the benchmark; seven passed after an additional allocation-budget regression check. A separate standalone SrTiO3 CLI smoke test produced five candidates with zero generation failures.
- `results/predictions/`, `results/physical/`, and `results/cifs/` preserve candidate manifests and ranked structures; `results/per_case.json` and `results/summary.json` preserve metrics.

## Scope and next validation

Implemented: constrained mapping, multiple mappings, soft neutral oxidation alternatives, exact PN neighborhoods, hybrid retrieval, radius compatibility, structural deduplication, safe substitutions, energy reranking, and limited relaxation. Optional symmetry voting and energy-rank fusion are available but are not validated by this experiment.

Not established: mixed-valence site assignments, bond-valence scoring, rigorous topology clustering, calibrated uncertainty, consistent thermochemical template quality, or absence of pretrained-model overlap. Vacuum-gap detection is only a diagnostic. Parameter selection and final acceptance should use a separate development/holdout split. Thermodynamic and dynamical stability require their own screening.

## Delivered-code robustness check

After generation, the existing permutation-count limit was moved before list allocation to prevent excessive memory use on unusually high-component formulas. No scoring, ordering or benchmark settings changed. The largest Testset 46 mapping space is 24 permutations, below the unchanged limit of 40,320; the guard does not affect these predictions. The exact core used for generation is preserved in `results/source_snapshot/core.py`, with before/after hashes in `results/post_benchmark_robustness_patch.json`.

## Export validation

4133 ranked CIF files passed composition and strict structure round-trip checks. The audit uses `CifParser(frac_tolerance=0)` to compare the coordinates actually stored in each file. The initial audit exposed the default parser rounding one relaxed coordinate toward an ideal rational value; disabling that parser transformation resolved the check without changing any prediction, CIF, or benchmark matching tolerance. The initial failure log is preserved.
