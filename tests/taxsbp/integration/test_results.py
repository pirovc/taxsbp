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
    
    def test_invalid_specialization(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_invalid_specialization.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_duplicated_spec.tsv"
        cfg.specialization="spec"
        cfg.bin_len=200

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf, missing_entries=3), "Input/Output files are inconsistent")
        # specific check - no duplicated specialization per taxid
        unique_taxid_spec = outf[["taxid","specialization"]].drop_duplicates()
        self.assertFalse(unique_taxid_spec.specialization.duplicated().any(), "Invalid specialization")

    def test_tax_missing(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_tax_missing.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_missing_tax.tsv"
        cfg.bin_len=200

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, outf, missing_entries=3), "Input/Output files are inconsistent")

    def test_tax_spec_missing_update(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_tax_spec_missing_update.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_missing_tax.tsv"
        cfg.update_file=self.base_dir+"data/bins_M-L_tax_spec_missing.tsv"
        cfg.specialization="spec"
        cfg.bin_len=200

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        #updf = parse_output(cfg.update_file)
        # Join update file on output
        #mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, outf, missing_entries=12), "Input/Output files are inconsistent")

    def test_update_not_merging_bins(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_update_merging_bins.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo.tsv"
        cfg.update_file=self.base_dir+"data/bins_LKM.tsv"
        cfg.bin_len=200

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        
        # specific check - binid 0 and 1 have to have 2 sequences each
        self.assertEqual(sum(mergedf["binid"]==0), 2, "Update failed")
        self.assertEqual(sum(mergedf["binid"]==1), 2, "Update failed")

    def test_update_merging_bins(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_update_merging_bins.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo.tsv"
        cfg.update_file=self.base_dir+"data/bins_LKM.tsv"
        cfg.bin_len=400
        # Output is complete when using allow_merge
        cfg.allow_merge=True

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)

        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
    
        # specific check - binid 0 is updated with contents from bin 1
        self.assertEqual(sum(outf["binid"]==0), 4, "Merge failed")
        self.assertTrue(outf[outf["binid"]==0]["seqid"].isin(["L","K","M","J"]).all(), "Merge failed")


if __name__ == '__main__':
    unittest.main()