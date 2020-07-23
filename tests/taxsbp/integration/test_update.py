import unittest
import taxsbp.taxsbp
import pandas as pd
import os, shutil
from utils import *

class TestUpdate(unittest.TestCase):
    
    base_dir = "tests/taxsbp/integration/"
    results_dir = base_dir + "results/test_update/"
    default_config = {"input_file": base_dir + "data/seqinfo.tsv",
                        "nodes_file": base_dir + "data/nodes.dmp"}
    
    @classmethod
    def setUpClass(self):
        shutil.rmtree(self.results_dir, ignore_errors=True)
        os.makedirs(self.results_dir)
    
    def test_basic(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_basic.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.bin_len=300
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        outf = pd.concat([updf,outf], ignore_index=True)

        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # specific check - 3 entries should be on bin 0
        self.assertEqual(sum(outf["binid"]==0), 3, "Update failed")
    
if __name__ == '__main__':
    unittest.main()