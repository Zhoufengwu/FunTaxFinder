# FunTaxFinder

FunTaxFinder is a Python command-line tool for screening functional microbial
ASVs from PICRUSt2 per-sequence KEGG Orthology (KO) prediction tables.

It is designed for workflows where a researcher wants to identify ASVs that
carry one or more target functional markers, extract their abundance table,
attach taxonomy, and optionally compare the matched taxa with FAPROTAX output.

## Main features

- Screen ASVs using one or multiple KO identifiers.
- Support `any`, `all`, and `at_least_n` matching rules.
- Extract a matched ASV abundance table from the original ASV table.
- Write matched KO matrices, binary hit matrices, hit summaries, and ASV ID lists.
- Attach taxonomy annotations and split the matched ASV table by rank.
- Optionally run FAPROTAX for cross-validation only.
- Optionally run PICRUSt2 when no existing KO prediction table is provided.

## Expected PICRUSt2 input

FunTaxFinder expects ASV-by-KO prediction files, for example:

- `combined_KO_predicted.tsv`
- `combined_KO_predicted.tsv.gz`
- `bac_KO_predicted.tsv`
- `bac_KO_predicted.tsv.gz`
- `arc_KO_predicted.tsv`
- `arc_KO_predicted.tsv.gz`

Do not use `pred_metagenome_unstrat.tsv.gz` as the primary input for this
tool, because that file is usually KO-by-sample rather than ASV-by-KO.

## Installation

For local development:

```bash
python -m pip install -e .
```

After installation:

```bash
funtaxfinder --help
funtaxfinder --version
```

## Minimal example

```bash
funtaxfinder \
  --asv-table examples/asv_table.tsv \
  --ko K02586 \
  --ko-predicted examples/combined_KO_predicted.tsv \
  --outdir demo_output
```

## Example with multiple KOs and taxonomy splitting

```bash
funtaxfinder \
  --asv-table examples/asv_table.tsv \
  --ko K02586 K02588 \
  --ko-predicted examples/combined_KO_predicted.tsv \
  --tax-table examples/taxonomy.tsv \
  --split-rank genus \
  --mode any \
  --outdir demo_output
```

## Example using a KO list file

```bash
funtaxfinder \
  --asv-table examples/asv_table.tsv \
  --ko-list-file examples/nitrogen_fixation_kos.txt \
  --ko-predicted examples/combined_KO_predicted.tsv \
  --mode at_least_n \
  --at-least-n 2 \
  --outdir demo_output
```

## Output files

The output directory contains:

| File | Description |
| --- | --- |
| `target_ko_matrix_all_asvs.tsv` | Target KO values for all ASVs found in the prediction file(s). |
| `matched_asv_ko_matrix.tsv` | Target KO values for matched ASVs only. |
| `matched_asv_ko_binary.tsv` | Binary KO hit matrix for matched ASVs. |
| `matched_asv_summary.tsv` | Matched KOs and hit counts for each ASV. |
| `matched_asv_ids.txt` | One matched ASV ID per line. |
| `sub_asv_table.tsv` | Original ASV abundance table filtered to matched ASVs. |
| `matched_asv_taxonomy.tsv` | Matched ASVs with taxonomy, when `--tax-table` is supplied. |
| `tax_split/` | Rank-specific ASV subtables, when `--split-rank` is supplied. |
| `run_summary.txt` | Key run settings and output paths. |

## Missing KO policy

By default, FunTaxFinder raises an error if any requested KO is absent from all
prediction files:

```bash
--missing-ko-policy error
```

This strict default is important for `all` and `at_least_n` modes. Otherwise,
an absent KO could be silently ignored and produce false positive matches.

If absent KOs should be treated as zero-valued columns, use:

```bash
--missing-ko-policy zero
```

## Optional PICRUSt2 run

If no existing KO prediction file is provided, FunTaxFinder can run PICRUSt2:

```bash
funtaxfinder \
  --asv-table asv_table.tsv \
  --rep-seqs rep_seqs.fasta \
  --ko K02586 K02588 \
  --outdir output
```

PICRUSt2 must already be installed and available in the environment.

## Optional FAPROTAX cross-validation

FAPROTAX is optional and does not affect the primary KO-based screening result.

```bash
funtaxfinder \
  --asv-table asv_table.tsv \
  --ko K02586 \
  --ko-predicted combined_KO_predicted.tsv.gz \
  --tax-table taxonomy.tsv \
  --cross-validate-faprotax \
  --faprotax-script /path/to/collapse_table.py \
  --faprotax-db /path/to/FAPROTAX.txt \
  --outdir output
```


## License

MIT License.
