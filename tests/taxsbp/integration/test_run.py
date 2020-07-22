import unittest
import taxsbp.taxsbp
import pandas as pd

class TestRun(unittest.TestCase):
    
    default_config = {"input_file": "sample_data/seqinfo.tsv",
                        "nodes_file": "sample_data/nodes.dmp",
                        "output_file": "tests/output_test_run1.tsv"}

    def test_basic(self):
        cfg = Config(**self.default_config)
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP finished successfuly")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

    def test_fragment(self):
        cfg = Config(**self.default_config)
        cfg.fragment_len=50
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP finished successfuly")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test
        self.assertTrue(outf.length.max()<=cfg.fragment_len, "Fragment bigger than expected")

    def test_fragment_overlap(self):
        cfg = Config(**self.default_config)
        cfg.fragment_len=22
        cfg.overlap_len=7
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP finished successfuly")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test
        self.assertTrue(outf.length.max()<=cfg.fragment_len+cfg.overlap_len, "Fragment+overlap bigger than expected")

def parse_files(cfg):
    inf = parse_input(cfg.input_file)
    outf = parse_output(cfg.output_file)
    return [inf, outf]

def parse_input(file):
    colums=['seqid', 'length', 'taxid', 'specialization']
    types={'seqid': 'str', 'length': 'uint64', 'taxid': 'str', 'specialization': 'str'}
    return pd.read_csv(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def parse_output(file):
    colums=['seqid', 'start',' end', 'length', 'taxid', 'binid', 'specialization']
    types={'seqid': 'str', 'start': 'uint64', 'end': 'uint64', 'length': 'uint64', 'taxid': 'str', 'binid': 'uint64', 'specialization': 'str'}
    return pd.read_csv(file, sep='\t', header=None, skiprows=0, names=colums, dtype=types)

def sanity_check(cfg, inf, outf):
    if inf.empty: 
        print("Input file is empty")
        return False
    if outf.empty: 
        print("Output file is empty")
        return False
    if not cfg.overlap_len:
        # check if input and output are the same
        if not inf.groupby(["seqid"]).sum()["length"].equals(outf.groupby(["seqid"]).sum()["length"]):
            print("seqids/lengths not consistent between input and output")
            return False
    else:
        # with overlap, length is increased: check seqids and lenghts separated
        if not inf["seqid"].isin(outf["seqid"]).all():
            print("seqids not consistent between input and output")
            return False
        extra_length = (outf.shape[0] - inf.shape[0]) * cfg.overlap_len
        if (outf.length.sum()-extra_length != inf.length.sum()):
            print("lengths not consistent between input and output")
            return False
    return True

class Config():

    def __init__(self,
        bin_exclusive: str=None, 
        bin_len: int=0, 
        bins: int=0, 
        fragment_len: int=0, 
        input_file: str=None,
        merged_file: str=None,
        nodes_file: str=None,
        overlap_len: int=0,
        output_file: str=None,
        pre_cluster: str=None,
        specialization: str=None,
        update_file: str=None,
        output_unique_seqid: bool=False):

        self.bin_exclusive=bin_exclusive
        self.bin_len=bin_len
        self.bins=bins
        self.fragment_len=fragment_len
        self.input_file=input_file
        self.merged_file=merged_file
        self.nodes_file=nodes_file
        self.overlap_len=overlap_len
        self.output_file=output_file
        self.pre_cluster=pre_cluster
        self.specialization=specialization
        self.update_file=update_file
        self.output_unique_seqid=output_unique_seqid

if __name__ == '__main__':
    unittest.main()