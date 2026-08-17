from taxsbp.taxsbp import taxsbp
import os
import pandas as pd

base_dir = os.path.dirname(__file__)
sample_input = f"{base_dir}/data/sample.tsv"
sample_tax   = f"{base_dir}/data/sample.tax"

def to_dataframe(bins):
    return pd.DataFrame(bins, columns = ["uid", "weigth", "node", "binid"]) 


def sanity_check(sample_input, bins_df, stats):
    
    # All entries got a bin
    input_uids = set()
    input_lines = 0
    with open(sample_input, "r") as si:
        for line in si:
            input_uids.add(line.split("\t")[0])
            input_lines += 1
    assert input_lines == bins_df.shape[0]
    assert input_uids == set(bins_df.uid)

    # Bin with max. weigth is below target
    assert stats["weigths"]["max"] <= stats["target_weigth"]

    return True
    
def test_standard():
    bins, stats = taxsbp(input_file=sample_input, taxonomy_file=sample_tax, stats=True)
    assert sanity_check(sample_input, to_dataframe(bins), stats)

def test_bin_len():
    # One bin for each entry
    bins, stats = taxsbp(input_file=sample_input, taxonomy_file=sample_tax, bin_len=100, stats=True)
    assert sanity_check(sample_input,  to_dataframe(bins), stats)
    assert stats["total_bins"] == 13

    bins, stats = taxsbp(input_file=sample_input, taxonomy_file=sample_tax, bin_len=199, stats=True)
    assert sanity_check(sample_input,  to_dataframe(bins), stats)
    # One bin for each entry
    assert stats["total_bins"] == 13