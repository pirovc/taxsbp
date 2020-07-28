import unittest
import taxsbp.taxsbp
import pandas as pd
import os, shutil
from utils import *

class TestResults(unittest.TestCase):
    
    base_dir = "tests/taxsbp/integration/"
    results_dir = base_dir + "results/test_results/"
    default_config = {"input_file": base_dir + "data/seqinfo.tsv",
                        "nodes_file": base_dir + "data/nodes.dmp"}
    
    @classmethod
    def setUpClass(self):
        shutil.rmtree(self.results_dir, ignore_errors=True)
        os.makedirs(self.results_dir)

    def test_bin_exclusive(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_exclusive.tsv"
        cfg.bin_exclusive="rank-4"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific check
        # all nodes reported are from rank-4 (and 1 for the entry M)
        self.assertTrue(outf["taxid"].isin(["4.1","4.2","4.3","4.4","4.5","4.6","1"]).all() , "Wrong node reported")

    def test_bin_exclusive_specialization(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_exclusive_specialization.tsv"
        cfg.specialization="spec"
        cfg.bin_len=200

        for bin_exclusive in ["rank-2","rank-3","rank-4","rank-5","spec"]:
        	cfg.bin_exclusive=bin_exclusive
	        # run check
	        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
	        inf, outf = parse_files(cfg)
	        # sanity check
	        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
	        # specific check
	        unique_binid = outf[["binid","specialization" if bin_exclusive=="spec" else "taxid"]].drop_duplicates()
	        # bin ids should be uniquely matched (not mixed)
	        self.assertEqual(unique_binid.shape[0], unique_binid.binid.max()+1, "Bins are not rank exclusive")

if __name__ == '__main__':
    unittest.main()