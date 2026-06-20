"""Command-line interface for FunTaxFinder."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from funtaxfinder import __version__
from funtaxfinder.core import (
    attach_taxonomy_to_hits,
    build_faprotax_input,
    build_hit_summary,
    build_target_ko_matrix,
    count_asv_rows,
    eprint,
    ensure_exists,
    filter_asv_table,
    filter_hit_asvs,
    find_prediction_files,
    load_tax_table,
    read_ko_list,
    run_faprotax_cross_validation,
    run_picrust2,
    split_sub_asv_by_tax,
    write_kv_summary,
)


def build_parser() -> argparse.ArgumentParser:
    """Create the FunTaxFinder argument parser."""

    parser = argparse.ArgumentParser(
        prog="funtaxfinder",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        description=(
            "Screen functional microbial ASVs from PICRUSt2 ASV-by-KO prediction "
            "tables, extract matched ASV abundance tables, and optionally attach "
            "taxonomy or run FAPROTAX cross-validation."
        ),
    )
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    input_group = parser.add_argument_group("Input")
    input_group.add_argument(
        "--asv-table",
        "--asv_table",
        dest="asv_table",
        required=True,
        help="ASV abundance table. The first column must contain ASV IDs.",
    )
    input_group.add_argument(
        "--ko",
        nargs="+",
        help="Target KO identifiers, for example K02586 K02588.",
    )
    input_group.add_argument(
        "--ko-list-file",
        "--ko_list_file",
        dest="ko_list_file",
        help="Text file with one target KO identifier per line.",
    )
    input_group.add_argument(
        "--ko-predicted",
        "--ko_predicted",
        dest="ko_predicted",
        nargs="+",
        help=(
            "Existing ASV-by-KO prediction file(s), such as "
            "combined_KO_predicted.tsv.gz or bac/arc_KO_predicted.tsv.gz."
        ),
    )
    input_group.add_argument(
        "--picrust2-dir",
        "--picrust2_dir",
        dest="picrust2_dir",
        help=(
            "Existing PICRUSt2 output directory. FunTaxFinder will search for "
            "combined_KO_predicted or bac/arc_KO_predicted files."
        ),
    )

    picrust_group = parser.add_argument_group("Optional PICRUSt2 run")
    picrust_group.add_argument(
        "--rep-seqs",
        "--rep_seqs",
        dest="rep_seqs",
        help=(
            "Representative sequence FASTA. Required only when no existing KO "
            "prediction file is provided and PICRUSt2 should be run automatically."
        ),
    )
    picrust_group.add_argument(
        "--picrust2-cmd",
        "--picrust2_cmd",
        dest="picrust2_cmd",
        default="picrust2_pipeline.py",
        help="PICRUSt2 command.",
    )
    picrust_group.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Thread count for the optional PICRUSt2 run.",
    )
    picrust_group.add_argument(
        "--run-picrust2-outdir",
        "--run_picrust2_outdir",
        dest="run_picrust2_outdir",
        default=None,
        help="Output directory for the optional PICRUSt2 run.",
    )

    screen_group = parser.add_argument_group("Screening")
    screen_group.add_argument(
        "--min-value",
        "--min_value",
        dest="min_value",
        type=float,
        default=0.0,
        help="KO hit threshold. A KO is considered present when its value is greater than this.",
    )
    screen_group.add_argument(
        "--mode",
        choices=["any", "all", "at_least_n"],
        default="any",
        help="Combination rule when multiple target KOs are provided.",
    )
    screen_group.add_argument(
        "--at-least-n",
        "--at_least_n",
        dest="at_least_n",
        type=int,
        default=1,
        help="Required number of hit KOs when --mode at_least_n is used.",
    )
    screen_group.add_argument(
        "--missing-ko-policy",
        "--missing_ko_policy",
        dest="missing_ko_policy",
        choices=["error", "zero"],
        default="error",
        help=(
            "How to handle requested KOs absent from all prediction files. "
            "'error' avoids silent false positives in all/at_least_n screening."
        ),
    )

    tax_group = parser.add_argument_group("Optional taxonomy")
    tax_group.add_argument(
        "--tax-table",
        "--tax_table",
        dest="tax_table",
        help="Optional taxonomy table used for matched-ASV annotation and rank splitting.",
    )
    tax_group.add_argument(
        "--tax-id-col",
        "--tax_id_col",
        dest="tax_id_col",
        default=None,
        help="ASV ID column name in the taxonomy table. By default, the first column is used.",
    )
    tax_group.add_argument(
        "--tax-column",
        "--tax_column",
        dest="tax_column",
        default=None,
        help="Column containing a full taxonomy lineage string.",
    )
    tax_group.add_argument(
        "--split-rank",
        "--split_rank",
        dest="split_rank",
        default=None,
        choices=["kingdom", "phylum", "class", "order", "family", "genus", "species"],
        help="Taxonomic rank used to split the matched ASV abundance table.",
    )

    faprotax_group = parser.add_argument_group("Optional FAPROTAX cross-validation")
    faprotax_group.add_argument(
        "--cross-validate-faprotax",
        "--cross_validate_faprotax",
        dest="cross_validate_faprotax",
        action="store_true",
        help="Run FAPROTAX after KO-based screening. This does not affect the primary screen.",
    )
    faprotax_group.add_argument(
        "--faprotax-python",
        "--faprotax_python",
        dest="faprotax_python",
        default=sys.executable,
        help="Python interpreter used to run FAPROTAX.",
    )
    faprotax_group.add_argument(
        "--faprotax-script",
        "--faprotax_script",
        dest="faprotax_script",
        default=None,
        help="Path to the FAPROTAX collapse_table.py script.",
    )
    faprotax_group.add_argument(
        "--faprotax-db",
        "--faprotax_db",
        dest="faprotax_db",
        default=None,
        help="Path to the FAPROTAX database file.",
    )
    faprotax_group.add_argument(
        "--faprotax-taxonomy-column-name",
        "--faprotax_taxonomy_column_name",
        dest="faprotax_taxonomy_column_name",
        default="taxonomy",
        help="Taxonomy column name written to the generated FAPROTAX input table.",
    )
    faprotax_group.add_argument(
        "--faprotax-comment-char",
        "--faprotax_comment_char",
        dest="faprotax_comment_char",
        default="#",
        help="FAPROTAX -c comment character.",
    )
    faprotax_group.add_argument(
        "--faprotax-report",
        "--faprotax_report",
        dest="faprotax_report",
        action="store_true",
        help="Write an additional FAPROTAX report file.",
    )

    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--outdir",
        required=True,
        help="Output directory.",
    )

    return parser


def run(args: argparse.Namespace) -> None:
    """Run the FunTaxFinder workflow from parsed CLI arguments."""

    ensure_exists(args.asv_table, "ASV abundance table")
    target_kos = read_ko_list(args.ko, args.ko_list_file)

    if args.mode == "at_least_n" and args.at_least_n < 1:
        raise ValueError("--at-least-n must be >= 1.")
    if args.mode == "at_least_n" and args.at_least_n > len(target_kos):
        raise ValueError(
            "--at-least-n cannot exceed the number of requested KO identifiers "
            f"({len(target_kos)})."
        )

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    prediction_files: list[str] = []

    if args.ko_predicted:
        prediction_files = [str(Path(x)) for x in args.ko_predicted]
        for pred_file in prediction_files:
            ensure_exists(pred_file, "KO prediction file")
        eprint(f"[INFO] Using user-provided KO prediction file(s): {prediction_files}")
    else:
        if args.picrust2_dir:
            prediction_files = find_prediction_files(args.picrust2_dir)
            if prediction_files:
                eprint(
                    "[INFO] Found KO prediction file(s) in the existing PICRUSt2 "
                    f"directory: {prediction_files}"
                )

        if not prediction_files:
            if not args.rep_seqs:
                raise ValueError(
                    "No KO prediction file was provided or found. Provide --rep-seqs "
                    "to run PICRUSt2 automatically."
                )

            run_out = args.run_picrust2_outdir or str(outdir / "picrust2_run")
            prediction_files = run_picrust2(
                asv_table=args.asv_table,
                rep_seqs=args.rep_seqs,
                outdir=run_out,
                threads=args.threads,
                picrust2_cmd=args.picrust2_cmd,
                stratified=True,
                verbose=True,
            )
            eprint(f"[INFO] PICRUSt2 produced KO prediction file(s): {prediction_files}")

    ko_matrix, found_kos, missing_kos = build_target_ko_matrix(
        prediction_files,
        target_kos,
        missing_ko_policy=args.missing_ko_policy,
    )
    eprint(f"[INFO] Requested KO identifiers: {target_kos}")
    eprint(f"[INFO] Found KO identifiers: {found_kos}")
    if missing_kos:
        eprint(f"[WARNING] Missing KO identifiers treated as zero: {missing_kos}")

    all_target_matrix_file = outdir / "target_ko_matrix_all_asvs.tsv"
    ko_matrix.to_csv(all_target_matrix_file, sep="\t", float_format="%.6f")

    hit_df = filter_hit_asvs(
        ko_matrix=ko_matrix,
        min_value=args.min_value,
        mode=args.mode,
        at_least_n=args.at_least_n,
    )

    hit_matrix_file = outdir / "matched_asv_ko_matrix.tsv"
    hit_df.to_csv(hit_matrix_file, sep="\t", float_format="%.6f")

    hit_binary = (hit_df > args.min_value).astype(int)
    hit_binary_file = outdir / "matched_asv_ko_binary.tsv"
    hit_binary.to_csv(hit_binary_file, sep="\t")

    hit_summary = build_hit_summary(hit_df, min_value=args.min_value)
    hit_summary_file = outdir / "matched_asv_summary.tsv"
    hit_summary.to_csv(hit_summary_file, sep="\t", index=False)

    hit_ids = list(hit_df.index)
    hit_ids_file = outdir / "matched_asv_ids.txt"
    with open(hit_ids_file, "w", encoding="utf-8") as handle:
        for asv_id in hit_ids:
            handle.write(f"{asv_id}\n")

    sub_asv_table_file = outdir / "sub_asv_table.tsv"
    written_n, written_ids = filter_asv_table(
        args.asv_table,
        keep_ids=hit_ids,
        output_file=sub_asv_table_file,
    )

    missing_from_asv_table = sorted(set(hit_ids) - written_ids)
    missing_asv_file = outdir / "matched_asvs_missing_from_abundance_table.txt"
    with open(missing_asv_file, "w", encoding="utf-8") as handle:
        for asv_id in missing_from_asv_table:
            handle.write(f"{asv_id}\n")
    if missing_from_asv_table:
        eprint(
            "[WARNING] "
            f"{len(missing_from_asv_table)} matched ASV(s) were absent from the abundance table."
        )

    total_asv_n = count_asv_rows(args.asv_table)

    matched_tax_df = None
    matched_tax_file = None
    tax_split_files: list[str] = []

    if args.tax_table:
        tax_df = load_tax_table(args.tax_table, tax_id_col=args.tax_id_col)
        matched_tax_df = attach_taxonomy_to_hits(
            hit_ids=hit_ids,
            tax_df=tax_df,
            tax_column=args.tax_column,
            split_rank=args.split_rank,
        )

        matched_tax_file = outdir / "matched_asv_taxonomy.tsv"
        matched_tax_df.to_csv(matched_tax_file, sep="\t", index=False)

        if args.split_rank:
            tax_split_files = split_sub_asv_by_tax(
                sub_asv_table_file=sub_asv_table_file,
                matched_tax_df=matched_tax_df,
                outdir=outdir,
                split_rank=args.split_rank,
            )
    elif args.split_rank:
        raise ValueError("--split-rank requires --tax-table.")

    faprotax_input_file = None
    faprotax_output_file = None
    faprotax_report_file = None
    faprotax_input_n = None

    if args.cross_validate_faprotax:
        if not args.tax_table:
            raise ValueError("--cross-validate-faprotax requires --tax-table.")
        if not args.faprotax_script:
            raise ValueError("--cross-validate-faprotax requires --faprotax-script.")
        if not args.faprotax_db:
            raise ValueError("--cross-validate-faprotax requires --faprotax-db.")

        faprotax_input_file = outdir / "faprotax_input.tsv"
        faprotax_input_file, faprotax_input_n = build_faprotax_input(
            sub_asv_table_file=sub_asv_table_file,
            matched_tax_df=matched_tax_df,
            output_file=faprotax_input_file,
            taxonomy_col_name=args.faprotax_taxonomy_column_name,
        )

        faprotax_output_file = outdir / "faprotax_output.tsv"
        if args.faprotax_report:
            faprotax_report_file = outdir / "faprotax_report.txt"

        run_faprotax_cross_validation(
            python_exe=args.faprotax_python,
            script_path=args.faprotax_script,
            db_path=args.faprotax_db,
            input_file=faprotax_input_file,
            output_file=faprotax_output_file,
            taxonomy_col=args.faprotax_taxonomy_column_name,
            comment_char=args.faprotax_comment_char,
            verbose=True,
            report_file=faprotax_report_file,
        )

        write_kv_summary(
            outdir / "faprotax_cross_validation_summary.txt",
            {
                "mode": "cross_validation_only",
                "input": str(faprotax_input_file),
                "input_asv_count": faprotax_input_n,
                "output": str(faprotax_output_file),
                "report": str(faprotax_report_file) if faprotax_report_file else "NO",
                "note": "FAPROTAX is used only for cross-validation and does not affect KO-based screening.",
            },
        )

    summary_file = outdir / "run_summary.txt"
    write_kv_summary(
        summary_file,
        {
            "software": "FunTaxFinder",
            "version": __version__,
            "input_asv_table": args.asv_table,
            "prediction_files_used": "; ".join(prediction_files),
            "target_kos_requested": ", ".join(target_kos),
            "target_kos_found": ", ".join(found_kos),
            "target_kos_missing": ", ".join(missing_kos) if missing_kos else "NO",
            "missing_ko_policy": args.missing_ko_policy,
            "min_value": args.min_value,
            "mode": args.mode,
            "at_least_n": args.at_least_n if args.mode == "at_least_n" else "NA",
            "total_asv_in_input_table": total_asv_n,
            "hit_asv_in_prediction": len(hit_ids),
            "hit_asv_written_to_sub_table": written_n,
            "matched_asv_missing_from_abundance_table": len(missing_from_asv_table),
            "tax_table_used": args.tax_table if args.tax_table else "NO",
            "tax_column_used": args.tax_column if args.tax_column else "AUTO",
            "split_rank": args.split_rank if args.split_rank else "NO",
            "tax_split_file_count": len(tax_split_files),
            "faprotax_cross_validation": "ON" if args.cross_validate_faprotax else "OFF",
            "faprotax_input_asv_count": faprotax_input_n if faprotax_input_n is not None else "NO",
            "output_sub_asv_table": str(sub_asv_table_file),
            "output_hit_ids": str(hit_ids_file),
            "output_hit_matrix": str(hit_matrix_file),
            "output_hit_binary": str(hit_binary_file),
            "output_hit_summary": str(hit_summary_file),
            "output_missing_asv_ids": str(missing_asv_file),
            "output_matched_taxonomy": str(matched_tax_file) if matched_tax_file else "NO",
            "output_faprotax_input": str(faprotax_input_file) if faprotax_input_file else "NO",
            "output_faprotax_output": str(faprotax_output_file) if faprotax_output_file else "NO",
            "output_faprotax_report": str(faprotax_report_file) if faprotax_report_file else "NO",
        },
    )

    eprint("")
    eprint("[INFO] FunTaxFinder finished successfully.")
    eprint(f"[INFO] Input ASV count: {total_asv_n}")
    eprint(f"[INFO] Matched ASV count in prediction table(s): {len(hit_ids)}")
    eprint(f"[INFO] Matched ASV count written to subtable: {written_n}")
    if args.split_rank:
        eprint(f"[INFO] Taxonomy split files: {len(tax_split_files)}")
    if args.cross_validate_faprotax:
        eprint(f"[INFO] FAPROTAX input ASV count: {faprotax_input_n}")
        eprint("[INFO] FAPROTAX cross-validation finished.")
    eprint(f"[INFO] Output directory: {outdir}")


def main(argv: list[str] | None = None) -> int:
    """CLI entry point."""

    parser = build_parser()
    args = parser.parse_args(argv)

    try:
        run(args)
    except Exception as exc:
        eprint(f"[ERROR] {exc}")
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
