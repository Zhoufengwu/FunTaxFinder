"""Core functions for FunTaxFinder.

FunTaxFinder screens amplicon sequence variants (ASVs) using PICRUSt2
per-sequence KEGG Orthology (KO) predictions. The primary expected prediction
tables are ASV-by-KO files such as ``combined_KO_predicted.tsv.gz``,
``bac_KO_predicted.tsv.gz``, or ``arc_KO_predicted.tsv.gz``.
"""

from __future__ import annotations

import csv
import gzip
import re
import shlex
import subprocess
import sys
from pathlib import Path
from typing import Iterable

import pandas as pd


RANK_ORDER = ["kingdom", "phylum", "class", "order", "family", "genus", "species"]

RANK_ALIASES = {
    "kingdom": ["kingdom", "domain", "superkingdom", "k", "d"],
    "phylum": ["phylum", "p"],
    "class": ["class", "c"],
    "order": ["order", "o"],
    "family": ["family", "f"],
    "genus": ["genus", "g"],
    "species": ["species", "s"],
}

RANK_PREFIXES = {
    "kingdom": ["k__", "d__", "k_", "d_", "kingdom__", "domain__", "superkingdom__"],
    "phylum": ["p__", "p_", "phylum__"],
    "class": ["c__", "c_", "class__"],
    "order": ["o__", "o_", "order__"],
    "family": ["f__", "f_", "family__"],
    "genus": ["g__", "g_", "genus__"],
    "species": ["s__", "s_", "species__"],
}


def eprint(*args: object, **kwargs: object) -> None:
    """Print a message to standard error."""

    print(*args, file=sys.stderr, **kwargs)


def ensure_exists(path: str | Path | None, description: str) -> None:
    """Raise ``FileNotFoundError`` when ``path`` does not exist."""

    if not path or not Path(path).exists():
        raise FileNotFoundError(f"{description} does not exist: {path}")


def open_text_auto(path: str | Path, mode: str = "rt"):
    """Open plain-text or gzip-compressed text files."""

    path = str(path)
    if path.endswith(".gz"):
        return gzip.open(path, mode, encoding="utf-8", newline="")
    return open(path, mode, encoding="utf-8", newline="")


def unique_preserve_order(items: Iterable[str]) -> list[str]:
    """Return unique values while preserving their first-seen order."""

    seen: set[str] = set()
    output: list[str] = []
    for item in items:
        if item not in seen:
            seen.add(item)
            output.append(item)
    return output


def normalize_ko_name(ko: object) -> str:
    """Normalize KO identifiers to uppercase bare identifiers such as K02586."""

    value = str(ko).strip()
    if value.lower().startswith("ko:"):
        value = value[3:]
    return value.upper()


def sanitize_filename(text: object, replacement: str = "_") -> str:
    """Make a taxonomy label safe for use in a file name."""

    value = str(text).strip()
    if not value:
        return "Unassigned"
    value = re.sub(r"[^\w\-.]+", replacement, value, flags=re.UNICODE)
    value = re.sub(r"_+", "_", value).strip("_")
    return value if value else "Unassigned"


def clean_tax_value(value: object) -> str:
    """Remove common rank prefixes and normalize empty taxonomy labels."""

    text = str(value).strip()
    text = re.sub(r"^(k|d|p|c|o|f|g|s)__", "", text, flags=re.IGNORECASE)
    text = re.sub(r"^(k|d|p|c|o|f|g|s)_", "", text, flags=re.IGNORECASE)
    text = re.sub(
        r"^(kingdom|domain|superkingdom|phylum|class|order|family|genus|species)__",
        "",
        text,
        flags=re.IGNORECASE,
    )
    text = re.sub(
        r"^(kingdom|domain|superkingdom|phylum|class|order|family|genus|species)_",
        "",
        text,
        flags=re.IGNORECASE,
    )
    text = text.strip()

    empty_labels = {
        "",
        "NA",
        "NaN",
        "nan",
        "None",
        "__",
        "_",
        "unclassified",
        "uncultured",
        "unidentified",
    }
    if text in empty_labels:
        return "Unassigned"
    return text


def count_asv_rows(asv_table: str | Path) -> int:
    """Count feature rows in an ASV table, excluding the header."""

    n_rows = 0
    with open_text_auto(asv_table, "rt") as handle:
        reader = csv.reader(handle, delimiter="\t")
        next(reader, None)
        for _ in reader:
            n_rows += 1
    return n_rows


def write_kv_summary(path: str | Path, data: dict[str, object]) -> None:
    """Write a simple key-value summary file."""

    with open(path, "w", encoding="utf-8") as handle:
        for key, value in data.items():
            handle.write(f"{key}: {value}\n")


def read_ko_list(
    ko_args: Iterable[str] | None = None,
    ko_file: str | Path | None = None,
) -> list[str]:
    """Read target KO identifiers from CLI arguments and/or a text file."""

    kos: list[str] = []

    if ko_args:
        kos.extend(str(x) for x in ko_args)

    if ko_file:
        ensure_exists(ko_file, "KO list file")
        with open(ko_file, "r", encoding="utf-8") as handle:
            for line in handle:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                kos.append(line)

    kos = [normalize_ko_name(x) for x in kos if str(x).strip()]
    kos = unique_preserve_order(kos)

    if not kos:
        raise ValueError("No target KO was provided. Use --ko or --ko-list-file.")

    return kos


def find_prediction_files(picrust2_dir: str | Path) -> list[str]:
    """Find supported PICRUSt2 ASV-by-KO prediction files in a directory."""

    root = Path(picrust2_dir)
    if not root.exists():
        return []

    all_files = [path for path in root.rglob("*") if path.is_file()]

    def lname(path: Path) -> str:
        return path.name.lower()

    combined = [
        path
        for path in all_files
        if lname(path) in {"combined_ko_predicted.tsv", "combined_ko_predicted.tsv.gz"}
    ]
    if combined:
        combined = sorted(combined, key=lambda path: (len(path.parts), str(path)))
        return [str(combined[0])]

    bac = [
        path
        for path in all_files
        if lname(path) in {"bac_ko_predicted.tsv", "bac_ko_predicted.tsv.gz"}
    ]
    arc = [
        path
        for path in all_files
        if lname(path) in {"arc_ko_predicted.tsv", "arc_ko_predicted.tsv.gz"}
    ]

    files = sorted(bac + arc, key=lambda path: (len(path.parts), str(path)))
    return [str(path) for path in files]


def run_picrust2(
    asv_table: str | Path,
    rep_seqs: str | Path,
    outdir: str | Path,
    threads: int = 1,
    picrust2_cmd: str = "picrust2_pipeline.py",
    stratified: bool = True,
    verbose: bool = True,
) -> list[str]:
    """Run PICRUSt2 and return detected ASV-by-KO prediction files."""

    ensure_exists(asv_table, "ASV abundance table")
    ensure_exists(rep_seqs, "Representative sequence file")

    outdir = Path(outdir)
    outdir.parent.mkdir(parents=True, exist_ok=True)

    if outdir.exists():
        existing = find_prediction_files(outdir)
        if existing:
            eprint(f"[INFO] Reusing existing PICRUSt2 KO prediction file(s): {existing}")
            return existing
        raise FileExistsError(
            "The PICRUSt2 output directory already exists, but no reusable KO "
            f"prediction file was found: {outdir}. Remove the directory or use "
            "a different --run-picrust2-outdir."
        )

    cmd = shlex.split(picrust2_cmd)
    cmd.extend(
        [
            "-s",
            str(rep_seqs),
            "-i",
            str(asv_table),
            "-o",
            str(outdir),
            "-p",
            str(threads),
        ]
    )
    if stratified:
        cmd.append("--stratified")
    if verbose:
        cmd.append("--verbose")

    eprint("[INFO] Running PICRUSt2:")
    eprint("       " + " ".join(shlex.quote(x) for x in cmd))

    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError("PICRUSt2 failed. Check the environment, database, and inputs.")

    prediction_files = find_prediction_files(outdir)
    if not prediction_files:
        raise FileNotFoundError(
            "PICRUSt2 finished, but no supported KO prediction file was found "
            f"in: {outdir}"
        )

    return prediction_files


def inspect_prediction_header(pred_file: str | Path) -> tuple[str, dict[str, str], list[str]]:
    """Inspect a PICRUSt2 prediction header and map normalized KO names."""

    header = pd.read_csv(pred_file, sep="\t", nrows=0, compression="infer").columns.tolist()

    if len(header) < 2:
        raise ValueError(f"Prediction file has too few columns: {pred_file}")

    id_col = header[0]
    ko_cols = header[1:]

    ko_map: dict[str, str] = {}
    for col in ko_cols:
        normalized = normalize_ko_name(col)
        if normalized not in ko_map:
            ko_map[normalized] = col

    return id_col, ko_map, ko_cols


def load_selected_ko_matrix(
    pred_file: str | Path,
    target_kos: list[str],
) -> tuple[pd.DataFrame | None, list[str], list[str]]:
    """Load only the requested KO columns from one prediction file."""

    id_col, ko_map, ko_cols = inspect_prediction_header(pred_file)
    found_kos = [ko for ko in target_kos if ko in ko_map]

    if not found_kos:
        preview = [normalize_ko_name(x) for x in ko_cols[:15]]
        return None, [], preview

    usecols = [id_col] + [ko_map[ko] for ko in found_kos]

    df = pd.read_csv(pred_file, sep="\t", compression="infer", usecols=usecols)

    rename_dict = {id_col: "ASV_ID"}
    for ko in found_kos:
        rename_dict[ko_map[ko]] = ko

    df = df.rename(columns=rename_dict)
    df["ASV_ID"] = df["ASV_ID"].astype(str)
    df = df.set_index("ASV_ID")

    for ko in found_kos:
        df[ko] = pd.to_numeric(df[ko], errors="coerce").fillna(0)

    return df, found_kos, [normalize_ko_name(x) for x in ko_cols[:15]]


def build_target_ko_matrix(
    pred_files: Iterable[str | Path],
    target_kos: list[str],
    missing_ko_policy: str = "error",
) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Build a merged ASV-by-target-KO matrix from one or more prediction files.

    ``missing_ko_policy`` controls requested KOs that are absent from all
    prediction files:

    - ``error``: raise an error.
    - ``zero``: add missing KO columns filled with zero.
    """

    if missing_ko_policy not in {"error", "zero"}:
        raise ValueError("--missing-ko-policy must be either 'error' or 'zero'.")

    frames: list[pd.DataFrame] = []
    found_union: list[str] = []
    preview_info: dict[str, list[str]] = {}

    for pred_file in pred_files:
        df, found_kos, preview = load_selected_ko_matrix(pred_file, target_kos)
        preview_info[str(pred_file)] = preview
        if df is not None and not df.empty:
            frames.append(df)
            for ko in found_kos:
                if ko not in found_union:
                    found_union.append(ko)

    missing_kos = [ko for ko in target_kos if ko not in found_union]

    if not frames:
        msg = ["None of the requested KO identifiers was found in the prediction file(s)."]
        for pred_file, preview in preview_info.items():
            msg.append(f"File: {pred_file}")
            msg.append(f"Example KO columns: {preview}")
        raise ValueError("\n".join(msg))

    if missing_kos and missing_ko_policy == "error":
        msg = [
            "The following requested KO identifiers were not found in any prediction file:",
            ", ".join(missing_kos),
            "Use --missing-ko-policy zero only if absent KOs should be treated as zero.",
        ]
        raise ValueError("\n".join(msg))

    merged = pd.concat(frames, axis=0, join="outer").fillna(0)

    if merged.index.duplicated().any():
        merged = merged.groupby(level=0).max()

    for ko in missing_kos:
        merged[ko] = 0.0

    merged = merged[target_kos]

    return merged, found_union, missing_kos


def filter_hit_asvs(
    ko_matrix: pd.DataFrame,
    min_value: float = 0.0,
    mode: str = "any",
    at_least_n: int = 1,
) -> pd.DataFrame:
    """Filter ASVs using KO abundance thresholds and combination logic."""

    if ko_matrix.empty:
        return ko_matrix.iloc[0:0].copy()

    binary = ko_matrix > min_value

    if mode == "any":
        mask = binary.any(axis=1)
    elif mode == "all":
        mask = binary.all(axis=1)
    elif mode == "at_least_n":
        mask = binary.sum(axis=1) >= at_least_n
    else:
        raise ValueError(f"Unknown screening mode: {mode}")

    return ko_matrix.loc[mask].copy()


def build_hit_summary(hit_df: pd.DataFrame, min_value: float = 0.0) -> pd.DataFrame:
    """Build a compact table listing matched KOs for each hit ASV."""

    if hit_df.empty:
        return pd.DataFrame(columns=["ASV_ID", "hit_KOs", "hit_count"])

    binary = (hit_df > min_value).astype(int)
    records = []

    for asv_id, row in binary.iterrows():
        hit_kos = [ko for ko, value in row.items() if value == 1]
        records.append(
            {
                "ASV_ID": asv_id,
                "hit_KOs": ",".join(hit_kos),
                "hit_count": len(hit_kos),
            }
        )

    summary_df = pd.DataFrame(records)
    summary_df = summary_df.sort_values(
        by=["hit_count", "ASV_ID"], ascending=[False, True]
    ).reset_index(drop=True)
    return summary_df


def filter_asv_table(
    asv_table: str | Path,
    keep_ids: Iterable[str],
    output_file: str | Path,
) -> tuple[int, set[str]]:
    """Write an ASV abundance subtable containing only ``keep_ids``."""

    keep_ids = {str(x) for x in keep_ids}
    written_ids: set[str] = set()
    written_n = 0

    with open_text_auto(asv_table, "rt") as fin, open(
        output_file, "w", encoding="utf-8", newline=""
    ) as fout:
        reader = csv.reader(fin, delimiter="\t")
        writer = csv.writer(fout, delimiter="\t")

        header = next(reader, None)
        if header is None:
            raise ValueError(f"ASV table is empty: {asv_table}")
        writer.writerow(header)

        for row in reader:
            if not row:
                continue
            asv_id = str(row[0]).strip()
            if asv_id in keep_ids:
                writer.writerow(row)
                written_ids.add(asv_id)
                written_n += 1

    return written_n, written_ids


def load_tax_table(tax_table: str | Path, tax_id_col: str | None = None) -> pd.DataFrame:
    """Load a taxonomy table and index it by ASV ID."""

    ensure_exists(tax_table, "Taxonomy table")

    df = pd.read_csv(tax_table, sep="\t", dtype=str, compression="infer")
    df = df.fillna("")

    if df.shape[1] < 2:
        raise ValueError("The taxonomy table must contain at least two columns.")

    if tax_id_col:
        if tax_id_col not in df.columns:
            raise ValueError(f"--tax-id-col was not found in the taxonomy table: {tax_id_col}")
        id_col = tax_id_col
    else:
        id_col = df.columns[0]

    df = df.rename(columns={id_col: "ASV_ID"})
    df["ASV_ID"] = df["ASV_ID"].astype(str)
    df = df.set_index("ASV_ID")

    return df


def infer_tax_string_column(tax_df: pd.DataFrame, preferred: str | None = None) -> str | None:
    """Infer the column containing a full taxonomy lineage string."""

    if preferred and preferred in tax_df.columns:
        return preferred

    candidates = [
        "taxonomy",
        "Taxonomy",
        "tax",
        "Tax",
        "Taxon",
        "taxon",
        "lineage",
        "Lineage",
        "Consensus.Lineage",
        "consensus_taxonomy",
        "classification",
    ]
    for col in candidates:
        if col in tax_df.columns:
            return col
    return None


def resolve_rank_column(tax_df: pd.DataFrame, split_rank: str | None) -> str | None:
    """Resolve a taxonomy rank to an existing column name."""

    if split_rank is None:
        return None

    split_rank = split_rank.lower()
    lower_map = {col.lower(): col for col in tax_df.columns}

    if split_rank in lower_map:
        return lower_map[split_rank]

    aliases = RANK_ALIASES.get(split_rank, [])
    for alias in aliases:
        if alias.lower() in lower_map:
            return lower_map[alias.lower()]

    return None


def build_taxonomy_string(
    row: pd.Series,
    tax_df: pd.DataFrame,
    tax_column: str | None = None,
) -> str:
    """Build a full taxonomy string from either a lineage or rank columns."""

    tax_col = infer_tax_string_column(tax_df, preferred=tax_column)
    if tax_col and tax_col in row.index:
        return str(row[tax_col]).strip()

    lower_map = {col.lower(): col for col in tax_df.columns}
    values: list[str] = []
    for rank in RANK_ORDER:
        found_col = None
        for alias in RANK_ALIASES[rank]:
            if alias.lower() in lower_map:
                found_col = lower_map[alias.lower()]
                break
        if found_col:
            value = str(row[found_col]).strip()
            if value:
                values.append(value)

    if values:
        return "; ".join(values)

    return ""


def parse_taxonomy_string_to_rank(tax_string: object, split_rank: str | None) -> str:
    """Extract one rank label from a semicolon- or pipe-delimited lineage."""

    if split_rank is None:
        return "Unassigned"

    split_rank = split_rank.lower()
    if not tax_string or not str(tax_string).strip():
        return "Unassigned"

    parts = [x.strip() for x in re.split(r"[;|]", str(tax_string)) if x.strip()]
    if not parts:
        return "Unassigned"

    prefixes = RANK_PREFIXES.get(split_rank, [])
    for part in parts:
        lowered = part.lower()
        for prefix in prefixes:
            if lowered.startswith(prefix.lower()):
                value = clean_tax_value(part[len(prefix) :].strip())
                return value if value else "Unassigned"

    if split_rank in RANK_ORDER:
        idx = RANK_ORDER.index(split_rank)
        if idx < len(parts):
            return clean_tax_value(parts[idx])

    return "Unassigned"


def attach_taxonomy_to_hits(
    hit_ids: Iterable[str],
    tax_df: pd.DataFrame,
    tax_column: str | None = None,
    split_rank: str | None = None,
) -> pd.DataFrame:
    """Attach taxonomy rows and split-group labels to matched ASVs."""

    subset = tax_df.reindex(list(hit_ids)).copy()
    subset = subset.fillna("")

    subset["taxonomy_string"] = subset.apply(
        lambda row: build_taxonomy_string(row, tax_df, tax_column=tax_column),
        axis=1,
    )

    if split_rank:
        rank_col = resolve_rank_column(tax_df, split_rank)
        if rank_col:
            subset["split_group"] = subset[rank_col].apply(clean_tax_value)
        else:
            subset["split_group"] = subset["taxonomy_string"].apply(
                lambda value: parse_taxonomy_string_to_rank(value, split_rank)
            )
    else:
        subset["split_group"] = "NA"

    subset.insert(0, "ASV_ID", subset.index.astype(str))
    return subset.reset_index(drop=True)


def split_sub_asv_by_tax(
    sub_asv_table_file: str | Path,
    matched_tax_df: pd.DataFrame | None,
    outdir: str | Path,
    split_rank: str,
) -> list[str]:
    """Split a matched ASV abundance table by a selected taxonomy rank."""

    if matched_tax_df is None or matched_tax_df.empty:
        return []

    split_dir = Path(outdir) / "tax_split"
    split_dir.mkdir(parents=True, exist_ok=True)

    sub_df = pd.read_csv(sub_asv_table_file, sep="\t")
    asv_col = sub_df.columns[0]

    group_map = dict(
        zip(
            matched_tax_df["ASV_ID"].astype(str),
            matched_tax_df["split_group"].astype(str),
        )
    )

    sub_df["__split_group__"] = sub_df[asv_col].astype(str).map(group_map).fillna("Unassigned")

    output_files: list[str] = []
    summary_records: list[dict[str, object]] = []

    for group, group_df in sub_df.groupby("__split_group__", dropna=False):
        clean_group = clean_tax_value(group)
        group_file = split_dir / f"sub_asv.{split_rank}.{sanitize_filename(clean_group)}.tsv"
        out_df = group_df.drop(columns=["__split_group__"]).copy()
        out_df.to_csv(group_file, sep="\t", index=False)
        output_files.append(str(group_file))
        summary_records.append(
            {
                "split_rank": split_rank,
                "group": clean_group,
                "asv_count": out_df.shape[0],
            }
        )

    summary_df = pd.DataFrame(summary_records).sort_values(
        by=["asv_count", "group"], ascending=[False, True]
    )
    summary_df.to_csv(split_dir / "tax_split_summary.tsv", sep="\t", index=False)

    return output_files


def build_faprotax_input(
    sub_asv_table_file: str | Path,
    matched_tax_df: pd.DataFrame | None,
    output_file: str | Path,
    taxonomy_col_name: str = "taxonomy",
) -> tuple[Path, int]:
    """Build a classical FAPROTAX input table from matched ASVs and taxonomy."""

    if matched_tax_df is None or matched_tax_df.empty:
        raise ValueError("FAPROTAX cross-validation requires taxonomy information.")

    sub_df = pd.read_csv(sub_asv_table_file, sep="\t")
    asv_col = sub_df.columns[0]

    tax_map = dict(
        zip(
            matched_tax_df["ASV_ID"].astype(str),
            matched_tax_df["taxonomy_string"].astype(str),
        )
    )

    sub_df[taxonomy_col_name] = sub_df[asv_col].astype(str).map(tax_map).fillna("")
    sub_df[taxonomy_col_name] = sub_df[taxonomy_col_name].astype(str).str.strip()

    sub_df = sub_df[sub_df[taxonomy_col_name] != ""].copy()

    if sub_df.empty:
        raise ValueError(
            "No matched ASV has a non-empty taxonomy string, so FAPROTAX cannot run."
        )

    cols = [col for col in sub_df.columns if col != taxonomy_col_name] + [taxonomy_col_name]
    sub_df = sub_df[cols]

    output_file = Path(output_file)
    sub_df.to_csv(output_file, sep="\t", index=False)
    return output_file, sub_df.shape[0]


def run_faprotax_cross_validation(
    python_exe: str,
    script_path: str | Path,
    db_path: str | Path,
    input_file: str | Path,
    output_file: str | Path,
    taxonomy_col: str = "taxonomy",
    comment_char: str = "#",
    verbose: bool = True,
    report_file: str | Path | None = None,
) -> Path:
    """Run FAPROTAX as an optional cross-validation step."""

    ensure_exists(script_path, "FAPROTAX script")
    ensure_exists(db_path, "FAPROTAX database")
    ensure_exists(input_file, "FAPROTAX input table")

    cmd = [
        python_exe,
        str(script_path),
        "-i",
        str(input_file),
        "-o",
        str(output_file),
        "-g",
        str(db_path),
        "-d",
        str(taxonomy_col),
        "-c",
        str(comment_char),
    ]
    if verbose:
        cmd.append("-v")
    if report_file:
        cmd.extend(["-r", str(report_file)])

    eprint("[INFO] Running FAPROTAX cross-validation:")
    eprint("       " + " ".join(shlex.quote(x) for x in cmd))

    result = subprocess.run(cmd, check=False)
    if result.returncode != 0:
        raise RuntimeError(
            "FAPROTAX cross-validation failed. Check the script, database, and taxonomy column."
        )

    return Path(output_file)
