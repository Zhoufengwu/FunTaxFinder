from pathlib import Path

import pandas as pd
import pytest

from funtaxfinder.core import (
    build_hit_summary,
    build_target_ko_matrix,
    filter_hit_asvs,
    normalize_ko_name,
    parse_taxonomy_string_to_rank,
)


DATA = Path(__file__).parent / "data"


def test_normalize_ko_name():
    assert normalize_ko_name("ko:k02586") == "K02586"
    assert normalize_ko_name(" K02588 ") == "K02588"


def test_build_target_ko_matrix_strict_missing_ko():
    with pytest.raises(ValueError, match="not found"):
        build_target_ko_matrix(
            [DATA / "combined_KO_predicted.tsv"],
            ["K02586", "K99999"],
            missing_ko_policy="error",
        )


def test_build_target_ko_matrix_zero_missing_ko():
    matrix, found_kos, missing_kos = build_target_ko_matrix(
        [DATA / "combined_KO_predicted.tsv"],
        ["K02586", "K99999"],
        missing_ko_policy="zero",
    )
    assert found_kos == ["K02586"]
    assert missing_kos == ["K99999"]
    assert matrix["K99999"].sum() == 0


def test_filter_hit_asvs_all_mode():
    matrix = pd.DataFrame(
        {
            "K02586": [1.0, 1.0, 0.0],
            "K02588": [1.0, 0.0, 1.0],
        },
        index=["ASV1", "ASV2", "ASV3"],
    )
    hits = filter_hit_asvs(matrix, min_value=0, mode="all")
    assert list(hits.index) == ["ASV1"]


def test_build_hit_summary():
    matrix = pd.DataFrame(
        {
            "K02586": [1.0, 0.0],
            "K02588": [1.0, 2.0],
        },
        index=["ASV1", "ASV2"],
    )
    summary = build_hit_summary(matrix, min_value=0)
    assert summary.loc[0, "ASV_ID"] == "ASV1"
    assert summary.loc[0, "hit_KOs"] == "K02586,K02588"
    assert summary.loc[0, "hit_count"] == 2


def test_parse_taxonomy_string_to_rank():
    lineage = "d__Bacteria; p__Proteobacteria; c__Gammaproteobacteria; g__Klebsiella"
    assert parse_taxonomy_string_to_rank(lineage, "phylum") == "Proteobacteria"
    assert parse_taxonomy_string_to_rank(lineage, "genus") == "Klebsiella"
