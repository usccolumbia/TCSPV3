# Algorithm and scope

TCSPV3 retains TCSP2.0's composition-to-template workflow, BERTOS, Matscholar embeddings, and 731,293-template source index. It is neither a de novo generator nor a newly trained composition model.

1. Reduce the query for chemistry and retrieve templates with matching reduced stoichiometric coefficients. The exact target reduced composition is excluded from templates, including an actual source-CIF check.
2. Rank formulas by both embedding/group distance and a Villars periodic-number neighborhood, drawing on the PN similarity principle and fourth-/sixth-neighbor search described by [Oran et al., PNcsp+](https://doi.org/10.1021/acs.jctc.6c00044); merge the routes under a fixed candidate budget.
3. Search coefficient-preserving element assignments and retain up to three near-optimal mappings instead of only one.
4. Average BERTOS per-species probabilities, enumerate up to four charge-neutral alternatives, and apply a soft mismatch penalty only when evidence is available. All-metal compositions omit the ionic penalty.
5. Penalize query/template radius incompatibility using nearest-shell coordination and compatible Shannon or atomic radii. Substitute all elements simultaneously, validate composition and geometry, and deduplicate structures with pymatgen StructureMatcher.
6. Optionally score the same candidate set using radius-scaled MatterSim single-point energies or limited relaxation. Lower energy within a composition is not proof of thermodynamic stability.

Default retrieval caps the valid pool at 100 and exports up to 20 distinct structures. `--formula-type full` (default) restricts output to the input formula's atom count; `--formula-type reduced` allows cell sizes represented by retrieved templates. The output `formula_units` refers to the output cell and does not imply a primitive cell. Historical benchmark results used unrestricted sizes, equivalent to `--formula-type reduced` in the current CLI.

TCSPV3 also validates source-index agreement, preserves source CIF SHA256 hashes, records failures, and offers optional symmetry voting and rank fusion. Optional features should not be presented as benchmark-validated. It does not provide calibrated uncertainty, site-resolved mixed valence, bond-valence validation, rigorous topology clustering, or an independently validated multi-potential ensemble.

The PN route is an adaptation of a published chemical-similarity idea, not an independently invented periodic-number scheme or a copy of the PNcsp+ implementation. TCSPV3 retains the TCSP template database and merges PN and embedding retrieval; PNcsp+ uses a different prototype inventory and ML evaluation workflow. See [PNcsp+ full text](https://pmc.ncbi.nlm.nih.gov/articles/PMC13085239/).
