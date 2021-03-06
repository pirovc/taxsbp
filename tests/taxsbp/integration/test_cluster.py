import unittest
import taxsbp.taxsbp
import pandas as pd
import os, shutil
from utils import *

class TestCluster(unittest.TestCase):
    
    base_dir = "tests/taxsbp/integration/"
    results_dir = base_dir + "results/test_cluster/"
    default_config = {"input_file": base_dir + "data/seqinfo.tsv",
                        "nodes_file": base_dir + "data/nodes.dmp"}
    
    @classmethod
    def setUpClass(self):
        shutil.rmtree(self.results_dir, ignore_errors=True)
        os.makedirs(self.results_dir)

    def test_basic(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_basic.tsv"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        
    def test_input_table(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_input_table.tsv"
        inf = parse_input(cfg.input_file)
        cfg.input_table=inf
        cfg.input_file=None

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        outf = parse_output(cfg.output_file)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        
    def test_nodes_missing(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_nodes_missing.tsv"
        cfg.nodes_file=self.base_dir+"data/nodes_without_4.6.dmp"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check should not pass
        self.assertTrue(sanity_check(cfg, inf, outf, missing_entries=3), "Input/Output files are inconsistent")
        # specific test
        self.assertTrue(outf.shape[0]<=inf.shape[0], "Did not skip missing nodes")
    
    def test_merged_file(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_merged_file.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_merged_nodes.tsv"
        cfg.merged_file=self.base_dir+"data/merged.dmp"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check 
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        
        # Run standard
        cfg_default = Config(**self.default_config)
        cfg_default.output_file=self.results_dir+"test_merged_file_default.tsv"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg_default)), "TaxSBP fails to run")
        inf_default, outf_default = parse_files(cfg_default)
        # sanity check 
        self.assertTrue(sanity_check(cfg_default, inf_default, outf_default), "Input/Output files are inconsistent")
        # specific test - compare if merged fixed first run and it's equal a normal run
        self.assertTrue(outf.equals(outf_default), "Results are different between runs using merged.dmp")

    def test_fragment(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_fragment.tsv"
        cfg.fragment_len=50
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - check if biggest fragment on output <= from requested fragment len.
        self.assertTrue(outf.length.max()<=cfg.fragment_len, "Fragment bigger than expected")

    def test_fragment_overlap(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_fragment_overlap.tsv"
        cfg.fragment_len=22
        cfg.overlap_len=7
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - check if biggest fragment on output <= from requested fragment len. + overlap len.
        self.assertTrue(outf.length.max()<=cfg.fragment_len+cfg.overlap_len, "Fragment+overlap bigger than expected")
    
    def test_specialization(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_specialization.tsv"
        cfg.specialization="spec"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - check if all specializations are on output
        self.assertTrue(inf["specialization"].isin(outf["specialization"]).all() , "Specialization not reported")

    def test_pre_cluster(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_pre_cluster.tsv"
        cfg.bin_len=200
        cfg.pre_cluster="rank-2"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - pre-clustered can and should be bigger than bin length
        self.assertTrue(outf.groupby(["binid"]).sum()["length"].max()>cfg.bin_len , "Pre-clustering failed")

    def test_bin_exclusive(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_exclusive.tsv"
        cfg.bin_len=1300
        cfg.bin_exclusive="rank-2"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - even when possible (bin_len=1300), using bin_exclusive should split the inputs into more bins
        self.assertTrue(outf["binid"].max()>0 , "Bin-exclusive clustering failed")
        unique_binid = outf[["binid","taxid"]].drop_duplicates()
        self.assertEqual(unique_binid.shape[0], unique_binid.binid.max()+1, "Bins are not rank exclusive")

    def test_bins(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bins.tsv"
        cfg.bins=13
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - matching number of bins with chosen value (not always but works in this case)
        self.assertEqual(outf["binid"].max()+1, cfg.bins, "Number of bins do not match")

    def test_bin_len(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_len.tsv"
        cfg.bin_len=215
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - check if bin_len parameters is limiting the size of the bins
        max_bin_len = outf.groupby(["binid"]).sum()["length"].max()
        self.assertTrue(max_bin_len<=cfg.bin_len, "Bin is bigger than requested")


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
            cfg.output_file=self.results_dir+"test_bin_exclusive_specialization_" + bin_exclusive + ".tsv"
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


    def test_all(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_all.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_merged_nodes.tsv"
        cfg.merged_file=self.base_dir+"data/merged.dmp"
        cfg.bin_exclusive="rank-4"
        cfg.bin_len=347 
        cfg.fragment_len=58 
        cfg.overlap_len=9
        cfg.pre_cluster="rank-2"
        cfg.specialization="strain"
  
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
    
    def test_pre_cluster_invalid_rank(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_pre_cluster_invalid_rank.tsv"
        cfg.bin_len=200
        cfg.pre_cluster="XY2418721Z"

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # specific test - pre-clustered can and should be bigger than bin length
        self.assertTrue(outf.groupby(["binid"]).sum()["length"].max()>cfg.bin_len , "Pre-clustering failed")
    
    def test_bin_exclusive_invalid_rank(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_exclusive.tsv"
        cfg.bin_len=1300
        cfg.bin_exclusive="XY421Z91241"

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # specific test - even when possible (bin_len=1300), using bin_exclusive should split the inputs into more bins
        self.assertTrue(outf["binid"].max()>0 , "Bin-exclusive clustering failed")
        unique_binid = outf[["binid","taxid"]].drop_duplicates()
        self.assertEqual(unique_binid.shape[0], unique_binid.binid.max()+1, "Bins are not rank exclusive")

    def test_missing_spec(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_missing_spec.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_missing_spec.tsv"
        cfg.bin_len=500
        cfg.specialization="strain"
        cfg.bin_exclusive="strain"

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # should have "specialiation-*" entries 
        self.assertTrue(outf["specialization"].map(lambda x: x.startswith("specialization-")).any() , "Specialization replacement failed")
        
if __name__ == '__main__':
    unittest.main()