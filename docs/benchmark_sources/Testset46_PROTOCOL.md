# TCSPV3 exploratory Testset 46 protocol

Frozen before the new evaluation, 2026-09-27. This is a standalone package, not the installed `tcsp3` or a CSP dispatcher change.

All 46 reduced compositions are excluded from the full 731,293-template TCSP index. Actual source CIF composition is checked before substitution. Prediction receives composition only; reference CIF, space group and energy are accessible only to the evaluator. No parameter tuning on reference recovery is permitted in this run.

Four retrieval arms use a cap of 100 valid instantiated structures and return at most 20 distinct candidates:

1. `baseline`: reproduce the prior top-200 electronegativity plus top-200 assignment union, single neutral BERTOS assignment, soft OS weight 0.5, and structural deduplication.
2. `pn`: exact Villars PN table; coefficient-constrained radius 4, radius 6 only if the radius-4 formula pool is empty; rank by the existing soft-OS embedding cost.
3. `union`: baseline formula union plus up to 200 PN formulas; same single-OS chemical score. Retrieve extra PN candidates by round-robin merging two independently ranked routes, counting shared candidates once.
4. `v3`: same routes, up to three coefficient-preserving near-optimal mappings (within 0.25 of the best score), probability-weighted neutral BERTOS alternatives (top three labels per species, at most four neutral assignments), plus a PN worst-pair penalty `0.1 * min(max_PN_distance/6, 2)` and a coordination-dependent radius mismatch penalty `0.25 * mean(abs(log(r_query/r_template)))`. Neutral evidence mass scales OS penalties; no feasible neutral prediction means no OS penalty. All-metal compositions disable OS penalties. Geometry comes only from the source template. Radius estimates use Shannon ionic radii when both states and coordination are supported, otherwise atomic radii for both species.

The initial candidate budget precedes geometry reranking. Duplicate equivalence uses StructureMatcher with ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True. Baseline preserves the previous same-SG acceleration for exact reproducibility; new arms check all previously retained candidates. Record source SHA256, mappings, score terms, radius coverage, composition validation, shortest contact, and a vacuum-gap diagnostic. Reject only invalid/disordered structures and contacts below 0.5 A. No hard radius, oxidation, formation-energy, or symmetry-majority filters.

MatterSim single-point reranks the same top-20 distinct candidates returned by `v3`. Compare raw structures with a separate composition-only volume correction (median query/template radius ratio, linear scale clipped to [0.85,1.18]). Score energy per atom only within a composition. Log forces, model exceptions and checkpoint identity. Failed energies sort after successful predictions; no reference-informed fallback.

For `v3` chemical ranking and raw single-point ranking, relax the first five candidates with MatterSim, at most 100 FIRE steps, fmax=0.05 eV/A, variable cell, without symmetry constraint. Shared candidates are relaxed once and reused. Report convergence and recovery separately; unfinished optimization is not convergence. The limited relaxed comparison is SM@1/5, not SM@20. A 180-second per-candidate timeout prevents pathological tasks from stalling the benchmark.

Report all 46 cases, SM@1/5/10/20 and SG@5, candidate-pool oracle recovery, paired gains/losses, runtime, failures, generated/unique counts, and relaxed convergence. The oracle is evaluation only. Earlier Testset 46 inspection makes this exploratory; potential and embedding training overlap is not excluded by template-composition exclusion.

Pre-evaluation implementation check: the separate SrTiO3 physics smoke test showed that ExpCellFilter's optimizer convergence can differ slightly from the raw atomic-force threshold. Store both flags; `converged` requires optimizer convergence AND maximum atomic force <= 0.05 eV/A. This clarification was made before any Testset 46 reference evaluation or physical-stage execution.

Optional production interfaces include symmetry-diversity voting and a second energy calculator; neither is a default improvement claim. A negative-formation-energy filter is unavailable without consistent source metadata and is not fabricated. Mixed-valence site assignment, learned uncertainty calibration, rigorous topology clustering and independent holdout validation remain future work.
