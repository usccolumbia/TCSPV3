# TCSPV3 original TCSP 180-case list — CPU2 results

| Subset | N | Ranking | Top-1 | Top-5 |
|---|---:|---|---:|---:|
| all180 | 180 | chemical | 115 (63.89%) | 153 (85.00%) |
| all180 | 180 | scaled_singlepoint | 116 (64.44%) | 150 (83.33%) |
| additional_ids | 154 | chemical | 98 (63.64%) | 128 (83.12%) |
| additional_ids | 154 | scaled_singlepoint | 96 (62.34%) | 125 (81.17%) |
| overlap_ids | 26 | chemical | 17 (65.38%) | 25 (96.15%) |
| overlap_ids | 26 | scaled_singlepoint | 20 (76.92%) | 25 (96.15%) |

Top-k means at least one StructureMatcher match among the first k ranked predictions; all cases remain in the denominator. Matching tolerances: ltol=0.2, stol=0.3, angle_tol=5, primitive_cell=True, scale=True, attempt_supercell=False.

All 180 benchmark compositions were excluded from retrieval, including actual source-CIF checks. Predictions used composition only, up to 100 instantiated candidates and up to 20 distinct candidates; MatterSim scored radius-scaled structures on CPU without relaxation. Top-5 CIFs were reread and checked against their JSON structures.

This run uses the exact GitHub TCSP data/180_testdata.csv at commit 4da9ee6eae151732b406cec7df81e16e9d0a7125, with all180 reference CIFs from the archived paper ground-truth ZIP. Input/reference_mapping.json records every mapping and hash. Template exclusion does not establish exclusion from BERTOS, embeddings, or MatterSim training. Cohort overlaps with the earlier46 run are reported using their actual counts. Figure1 reports original TCSP2.0 Top1=68.3%; comparison remains affected by relaxation and template-exclusion differences.

Energy errors: 0. Audited source instances: 17883; distinct source CIFs: 14037; exported CIFs: 1710.

See PROTOCOL.md, input_manifest.json, submission.json, per_case.csv and summary.json for frozen settings and provenance.

## Transfer verification

All 180 prediction/physics pairs, 1710 exported CIFs, 180 archived reference CIFs, source snapshots and 540 stage runtime records were verified after transfer. CPU device, two Torch threads and one identical model state were confirmed for every case. 21 cases returned fewer than five candidates; 0 returned none. All stay in the denominator.

## Figure 1 comparison

The source CSV exactly matches all180 IDs in the archived original 180 list. Original TCSP2.0 Figure1 StructureMatcher Top1 is 68.3%; TCSPV3 chemical Top1 is 63.89% and radius-scaled single-point Top1 is 64.44%. The latter is -3.86 percentage points relative to the published value. Top5 must not be compared to a Top1 figure.

The paper relaxes predictions and accepts multiple reference polymorphs, while this frozen TCSPV3 run uses single-point energy scoring and one archived reference per target. All180 compositions are excluded from templates here, whereas the paper describes exact-entry exclusion. Identical target IDs alone do not make these protocols identical.

The previous curated180 result (50% Top1,75% Top5) used a different dataset with only31 shared IDs; keep that result separate.

## Paired ranking changes

Chemical to radius-scaled MatterSim ranking: Top1 gains 30 and losses 29 (net +1/180); Top5 gains 6 and losses 9 (net -3/180). Thus the best tested Top5 ranking on this original 180 set is chemical-only at 153/180 (85.0%); the best tested Top1 ranking is scaled single-point at 116/180 (64.44%). See paired_changes.json for IDs.

Slurm arrays 33438/33439/33440 and finalizer 33441 completed; every task requested partition=cpu, cpus-per-task=2. Finalizer exit code 0:0 and full state are saved in final_verification.json.
