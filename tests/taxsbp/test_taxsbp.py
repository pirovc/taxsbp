import os

import pandas as pd
import pytest

from taxsbp.taxsbp import taxsbp

base_dir = os.path.dirname(__file__)
sample_input = f"{base_dir}/data/sample.tsv"
sample_tax = f"{base_dir}/data/sample.tax"


def to_dataframe(bins):
    return pd.DataFrame(bins, columns=["uid", "weigth", "node", "binid"])


def get_sorted_list_uids(bins_df):
    return sorted(
        ["".join(sorted(b)) for b in bins_df.groupby("binid").sum().uid.to_list()]
    )


def sanity_check(sample_input, bins_df, stats):
    # All entries of the input should have a bin assignment
    input_uids = set()
    input_lines = 0
    with open(sample_input, "r") as si:
        for line in si:
            input_uids.add(line.split("\t")[0])
            input_lines += 1
    assert input_lines == bins_df.shape[0]
    assert input_uids == set(bins_df.uid)

    # Bin with max. weigth is below target (or minimum, if smaller then target) on stats
    # target_w = max(stats["target_weigth"], stats["weigths"]["min"])
    # assert stats["weigths"]["max"] <= target_w
    # Double check the output file
    # assert bins_df.groupby(["binid"]).sum()["weigth"].max() <= target_w

    return True


def test_standard():
    bins, stats = taxsbp(input_file=sample_input, taxonomy_file=sample_tax, stats=True)
    assert sanity_check(sample_input, to_dataframe(bins), stats)


@pytest.mark.parametrize(
    "bin_len, expected_bins_content",
    [
        (
            10,
            ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        ),  # One bin for each entry
        (100, ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]),
        (199, ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]),
        (
            200,
            ["AB", "CD", "EF", "GH", "IJ", "KL", "M"],
        ),  # Entries are grouped into bins
        (300, ["AB", "CD", "EFG", "HIM", "JKL"]),
        (500, ["ABCD", "EFGHI", "JKLM"]),
        (700, ["ABCDJKL", "EFGHIM"]),
        (1300, ["ABCDEFGHIJKLM"]),
        (2000, ["ABCDEFGHIJKLM"]),
    ],
)
def test_bin_len(bin_len, expected_bins_content):
    bins, stats = taxsbp(
        input_file=sample_input, taxonomy_file=sample_tax, bin_len=bin_len, stats=True
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)


@pytest.mark.parametrize(
    "n_bins, expected_bins_content",
    [
        (
            1,
            ["ABCDEFGHIJKLM"],
        ),
        (13, ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]),
    ],
)
def test_n_bins(n_bins, expected_bins_content):
    bins, stats = taxsbp(
        input_file=sample_input, taxonomy_file=sample_tax, n_bins=n_bins, stats=True
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)


@pytest.mark.parametrize(
    "pre_cluster, bin_len, expected_bins_content",
    [
        (
            "rank-4",
            100,
            ["AB", "C", "D", "EFGH", "I", "JKL", "M"],
        ),
        (
            "rank-4",
            200,
            ["AB", "CD", "EFGH", "IM", "JKL"],
        ),
        (
            "rank-4",
            400,
            ["ABCD", "EFGH", "IJKL", "M"],
        ),
        (
            "rank-2",
            400,
            ["ABCD", "EFGHI", "JKLM"],
        ),
        (
            "rank-1",
            400,
            ["ABCDEFGHIJKLM"],
        ),
        (
            "rank-5",
            100,
            ["A", "B", "C", "D", "EF", "G", "H", "I", "J", "K", "L", "M"],
        ),
        (
            "leaves",
            100,
            ["AB", "C", "D", "EF", "G", "H", "I", "JKL", "M"],
        ),
        (
            "leaves",
            200,
            ["AB", "CD", "EF", "GH", "IM", "JKL"],
        ),
    ],
)
def test_pre_cluster(pre_cluster, bin_len, expected_bins_content):
    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        pre_cluster=pre_cluster,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)


@pytest.mark.parametrize(
    "pre_cluster_file, bin_len, expected_bins_content",
    [
        (
            ["ABKL"],
            100,
            ["ABKL", "C", "D", "E", "F", "G", "H", "I", "J", "M"],
        ),
        (
            ["ABKL"],
            700,
            ["ABCDKL", "EFGHIJM"],
        ),
        (
            ["ABC", "DEF", "GHI", "JKLM"],
            200,
            ["ABC", "DEF", "GHI", "JKLM"],
        ),
        (
            ["AC", "BD"],
            300,
            ["AC", "BD", "EFG", "HIM", "JKL"],
        ),
        (
            ["AC", "BD"],
            600,
            ["ABCD", "EFGHI", "JKLM"],
        ),
        (
            ["ABCDEFGHIJKLM"],
            600,
            ["ABCDEFGHIJKLM"],
        ),
        (
            ["A", "B", "C", "J"],
            300,
            ["AB", "CD", "EFG", "HIM", "JKL"],
        ),
    ],
)
def test_pre_cluster_file(pre_cluster_file, bin_len, expected_bins_content):

    pc_file = f"{base_dir}/pre_cluster_file.tmp"
    with open(pc_file, "w") as pcfile:
        for line in pre_cluster_file:
            print(*list(line), sep="\t", file=pcfile)

    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        pre_cluster_file=pc_file,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)
    os.remove(pc_file)


@pytest.mark.parametrize(
    "bin_exclusive, bin_len, expected_bins_content",
    [
        (
            "rank-4",
            100,
            ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        ),
        (
            "rank-4",
            200,
            ["AB", "C", "D", "EF", "GH", "I", "J", "KL", "M"],
        ),
        (
            "rank-4",
            800,
            ["AB", "C", "D", "EFGH", "I", "JKL", "M"],
        ),
        (
            "rank-2",
            200,
            ["AB", "CD", "E", "FG", "HI", "JM", "KL"],
        ),
        (
            "rank-5",
            1000,
            ["ABCDHIJKLM", "EF", "G"],
        ),
        (
            "leaves",
            100,
            ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"],
        ),
        (
            "leaves",
            800,
            ["AB", "C", "D", "EF", "G", "H", "I", "JKL", "M"],
        ),
    ],
)
def test_bin_exclusive(bin_exclusive, bin_len, expected_bins_content):
    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        bin_exclusive=bin_exclusive,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)


@pytest.mark.parametrize(
    "bin_exclusive_file, bin_len, expected_bins_content",
    [
        (
            ["ABIJ"],
            200,
            ["AB", "CD", "EF", "GH", "IJ", "KL", "M"],
        ),
        (
            ["ABIJ"],
            400,
            ["ABIJ", "CDKL", "EFGH", "M"],
        ),
        (
            ["ABIJ"],
            900,
            ["ABIJ", "CDEFGHKLM"],
        ),
        (
            ["AEJ", "MDB", "I", "K"],
            400,
            ["AEJ", "BDM", "C", "FGHL", "I", "K"],
        ),
    ],
)
def test_bin_exclusive_file(bin_exclusive_file, bin_len, expected_bins_content):

    be_file = f"{base_dir}/bin_exclusive_file.tmp"
    with open(be_file, "w") as befile:
        for line in bin_exclusive_file:
            print(*list(line), sep="\t", file=befile)

    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        bin_exclusive_file=be_file,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)
    os.remove(be_file)


@pytest.mark.parametrize(
    "pre_cluster, bin_exclusive, bin_len, expected_bins_content",
    [
        (
            "rank-4",
            "rank-3",
            100,
            ["AB", "C", "D", "EFGH", "I", "JKL", "M"],
        ),
        (
            "rank-4",
            "rank-3",
            900,
            ["AB", "CD", "EFGH", "IJKLM"],
        ),
        (
            "rank-4",
            "rank-2",
            900,
            ["ABCD", "EFGHI", "JKLM"],
        ),
        (
            "rank-1",
            "rank-4",  # # cannot be bin exclusive, partial pre-cluster overlap
            500,
            ["ABCDEFGHIJKLM"],
        ),
        (
            "rank-5",
            "rank-4",
            500,
            ["AB", "C", "D", "EFGH", "I", "JKL", "M"],
        ),
    ],
)
def test_pre_cluster_bin_exclusive(
    pre_cluster, bin_exclusive, bin_len, expected_bins_content
):
    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        pre_cluster=pre_cluster,
        bin_exclusive=bin_exclusive,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)


@pytest.mark.parametrize(
    "pre_cluster_file, bin_exclusive_file, bin_len, expected_bins_content",
    [
        (
            ["AB", "CD"],
            ["ABCDE"],
            700,
            ["ABCDE", "FGHIJKL", "M"],
        ),
        (
            ["ABCD"],
            ["AB"],  # cannot be bin exclusive, partial overlap with pre-cluster
            200,
            ["ABCD", "EF", "GH", "IJ", "KL", "M"],
        ),
        (
            ["EG", "FH"],
            ["ABCDM", "FH"],
            800,
            ["ABCDM", "EGIJKL", "FH"],
        ),
        (
            ["AB"],
            ["AB", "CD"],
            300,
            ["AB", "CD", "EFG", "HIM", "JKL"],
        ),
        (
            ["AB"],
            ["BC"],  # cannot be bin exclusive since AB is pre-clustered
            200,
            ["AB", "CD", "EF", "GH", "IJ", "KL", "M"],
        ),
        (
            ["AB"],
            ["ABC"],
            200,
            ["AB", "C", "DM", "EF", "GH", "IJ", "KL"],
        ),
        (
            ["AB"],
            ["ABIM"],
            200,
            ["AB", "CD", "EF", "GH", "IM", "J", "KL"],
        ),
        (
            ["AB", "CD", "EF"],
            ["CDEFJK"],
            800,
            ["ABGHILM", "CDEFJK"],
        ),
    ],
)
def test_pre_cluster_file_bin_exclusive_file(
    pre_cluster_file, bin_exclusive_file, bin_len, expected_bins_content
):

    pc_file = f"{base_dir}/pre_cluster_file.tmp"
    with open(pc_file, "w") as pcfile:
        for line in pre_cluster_file:
            print(*list(line), sep="\t", file=pcfile)

    be_file = f"{base_dir}/bin_exclusive_file.tmp"
    with open(be_file, "w") as befile:
        for line in bin_exclusive_file:
            print(*list(line), sep="\t", file=befile)

    bins, stats = taxsbp(
        input_file=sample_input,
        taxonomy_file=sample_tax,
        pre_cluster_file=pc_file,
        bin_exclusive_file=be_file,
        bin_len=bin_len,
        stats=True,
    )
    bins_df = to_dataframe(bins)
    assert sanity_check(sample_input, bins_df, stats)
    assert get_sorted_list_uids(bins_df) == expected_bins_content
    assert stats["total_bins"] == len(expected_bins_content)
    os.remove(pc_file)
    os.remove(be_file)
