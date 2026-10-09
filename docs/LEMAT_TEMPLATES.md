# Use LeMat as a TCSPV3 template library

LeMat templates replace the **retrieval index and source CIF store**, not the
Matscholar embeddings or BERTOS oxidation-state model. This is index preparation,
not TCSPV3 retraining. The same workflow works for another authorized CIF collection.

## 1. Obtain and freeze the input CIFs

Use the official [LeMat-Bulk](https://huggingface.co/datasets/LeMaterial/LeMat-Bulk)
or [LeMat-BulkUnique](https://huggingface.co/datasets/LeMaterial/LeMat-BulkUnique)
dataset and follow its license/attribution requirements. Record the dataset
revision, subset, export procedure and any filters. These datasets are distributed
as structural records/Parquet, not a ready-made TCSPV3 pickle. Export your chosen
records to **one ordered, three-dimensionally periodic structure per `.cif` file**.
Keep stable upstream IDs in filenames and preserve a separate original-ID/source
manifest. Do not treat a Hugging Face split named `train` as your benchmark's
leakage-controlled training split.

If you already have an authorized individual-CIF folder such as
`/path/to/lematerial_cifs`, start at step 2. A multi-block `lematerial.cif` file
must be split into individual structures first; the builder rejects multi-structure
files rather than silently indexing the first block. CIFs do not generally encode
2D/3D dimensionality: remove non-bulk records during export using dataset metadata.

For the official Parquet records, the repository includes a streaming exporter.
Install its optional dependency, obtain a full commit SHA from the dataset's
Hugging Face history, and start with 1,000 records:

```bash
python -m pip install datasets
python scripts/export_lemat_cifs.py \
  --dataset LeMaterial/LeMat-BulkUnique --subset unique_pbe \
  --revision <40-character-dataset-commit-sha> \
  --max-records 1000 --outdir /path/to/lemat_export_pilot
```

The placeholders must be replaced before execution. The input to step 2 is
`/path/to/lemat_export_pilot/cifs`. Use a new directory and `--max-records 0`
only when explicitly exporting the entire frozen subset; expect substantial
download, disk and runtime costs. For LeMat-Bulk use `--subset compatible_pbe`
with its own pinned revision. `--exclude-ids heldout_source_ids.json` excludes
original `immutable_id` values during export. The exporter retains duplicate IDs
as different record-numbered files, logs failures, and writes
`EXPORT_COMPLETE.json` and `export_metadata.jsonl`; those map generated filenames
back to the upstream records. No dataset records are downloaded during the
builder's synthetic tests, and a full LeMat export has not been qualified here.

For official records, reconstruct geometry from `lattice_vectors`,
`cartesian_site_positions` (**Cartesian**, not fractional), and `species_at_sites`.
Resolve species labels through the record's `species` JSON when necessary; reject
mixed occupancy rather than choosing an element arbitrarily. Do not infer site
counts from the reduced or descriptive formula columns. Consult the
[official data-field descriptions](https://huggingface.co/datasets/LeMaterial/LeMat-BulkUnique/blob/main/README.md).

LeMat-BulkUnique removes duplicates according to its own fingerprint procedure.
It is not a guarantee that every crystal prototype is unique or thermodynamically
stable. Stability/energy filtering must be an explicit, separately documented choice.

## 2. Build the native index and primitive CIF store in parallel

Install TCSPV3's dependencies first (`python -m pip install -e .`). Run from this
repository, using a new output folder **outside** your original CIF tree:

```bash
python scripts/build_template_index.py \
  --cif-dir /path/to/lematerial_cifs \
  --outdir /path/to/lemat_tcspv3_library \
  --workers 16 \
  --primitive-tolerance 0.01 \
  --source-revision 'dataset=LeMaterial/LeMat-BulkUnique; subset=unique_pbe; revision=<your-frozen-revision>'
```

Use a small, copied input subset for an initial pilot. For Slurm, request CPU
resources and set `--workers` no higher than `SLURM_CPUS_PER_TASK`. This builder
does not require GPUs. Large libraries still need RAM for the native index and
offset map, disk for the primitive CIF archive, and time for CIF parsing; parallel
workers do not remove these costs. There is no automatic resume: if interrupted,
retain the partial output as diagnostics and use a new output folder.

The builder checks ordered structures and finite geometry, removes oxidation-state
labels, and derives the actual primitive cell at the recorded tolerance. It does
not idealize symmetry, pad supercells, impose an atom-count cap, relax structures,
or select only stable materials. Original source CIFs remain untouched. Each
accepted template records its original path/ID, source SHA256, exported SHA256,
original atom count and actual primitive formula/site count.

The output contains:

- `template_index.pkl`: native TCSPV3 anonymous-coefficient buckets mapping to
  `(template_id, primitive_cell_formula)` pairs.
- `data/source_cifs.bin` and `data/source_cifs.offsets.pkl`: the matching primitive
  CIF archive and read-only mmap offsets; no millions-of-files output tree is needed.
- `template_metadata.jsonl`: original-to-exported template provenance.
- `rejections.jsonl`: parse/geometry errors and explicit exclusions.
- `BUILD_COMPLETE.json`: accepted/rejected/excluded counts, settings and hashes.

**Use the repository's builder, not LeMat's anonymous-formula strings directly.**
Its keys follow `tcspv3.core.props()` and the predictor's native key parser,
including pymatgen's elemental/diatomic/peroxide conventions. The index is not
a CIF archive; it must be used with its matching store. Load only trusted,
locally generated pickle files. Structural prototype deduplication is not performed
by this builder; repeated prototypes/polymorphs are retained, while TCSPV3 performs
candidate deduplication during prediction.

Optional benchmark exclusions are JSON lists:

```bash
python scripts/build_template_index.py \
  --cif-dir /path/to/approved_training_validation_cifs \
  --outdir /path/to/lemat_excluded_library --workers 16 \
  --exclude-ids /path/to/heldout_source_ids.json \
  --exclude-formulas /path/to/heldout_compositions.json
```

ID exclusions match original filename stems. Formula exclusions remove the entire
reduced composition, not just one polymorph. These controls do not automatically
discover cross-database duplicates or prove temporal independence. Apply the
study's declared leakage/time cutoff rules before indexing and audit them separately.

## 3. Predict using the new library

Keep the existing `tcspv3/database` auxiliary assets installed: Matscholar
embeddings plus BERTOS tokenizer/config/weights are still necessary. Explicit
`--cif-store` overrides the old default CIF archive, so LeMat IDs resolve against
the **new** store rather than the original 731K database:

```bash
python tcsp.py --formula SrTiO3 --formula-type full \
  --database tcspv3/database \
  --index /path/to/lemat_tcspv3_library/template_index.pkl \
  --cif-store /path/to/lemat_tcspv3_library/data/source_cifs.bin \
  --topk 20 --output results/lemat_SrTiO3

python tcspx.py --csv examples/formulas.csv --workers 4 --formula-type full \
  --database tcspv3/database \
  --index /path/to/lemat_tcspv3_library/template_index.pkl \
  --cif-store /path/to/lemat_tcspv3_library/data/source_cifs.bin \
  --output results/lemat_batch
```

For known targets, feed their **actual primitive-cell formulas**, for example
`Co2Te2` when the desired primitive cell contains four atoms, not `CoTe`.
`--formula-type full` enforces the specified cell atom count; it does not infer Z.
Use `--formula-type reduced` only when deliberately allowing the library's different
cell sizes. TCSPV3 excludes same-reduced-composition templates at inference and
rechecks source CIF compositions. Fewer than 20 distinct eligible candidates can
be returned; do not pad or relabel a short pool as full coverage.

The prediction records index/source hashes and template IDs. Keep the matching
`BUILD_COMPLETE.json` and metadata next to your results. Optional MatterSim ranking
is a separate stage; use the [official checkpoint-download instructions](../README.md#optional-mattersim-stage).
Changing the template library is a new experiment: existing 731K benchmark
accuracies are not LeMat results, and BERTOS/MatterSim pretraining provenance is
not changed merely by using a new template set.
