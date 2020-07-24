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
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        outf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")


    def test_nodes_missing(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_nodes_missing.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.nodes_file=self.base_dir+"data/nodes_without_5.2.dmp"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check should not pass
        self.assertFalse(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test
        self.assertTrue(mergedf.shape[0]<=inf.shape[0], "Did not skip missing nodes")
    
    def test_merged_file(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_merged_file.tsv"
        cfg.input_file=self.base_dir+"data/seqinfo_merged_nodes.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.merged_file=self.base_dir+"data/merged.dmp"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check 
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        
        # Run standard
        cfg_default = Config(**self.default_config)
        cfg_default.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg_default.output_file=self.results_dir+"test_merged_file_default.tsv"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg_default)), "TaxSBP fails to run")
        inf_default, outf_default = parse_files(cfg_default)
        updf_default = parse_output(cfg_default.update_file)
        # Join update file on output
        mergedf_default = pd.concat([updf_default,outf_default], ignore_index=True)
        # sanity check 
        self.assertTrue(sanity_check(cfg_default, inf_default, mergedf_default), "Input/Output files are inconsistent")
        # specific test - compare if merged fixed first run and it's equal a normal run
        self.assertTrue(mergedf.equals(mergedf_default), "Results are different between runs using merged.dmp")

    def test_fragment(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_fragment.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.fragment_len=33
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - check if biggest fragment on output <= from requested fragment len.
        self.assertTrue(outf.length.max()<=cfg.fragment_len, "Fragment bigger than expected")

    def test_fragment_overlap(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_fragment_overlap.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.fragment_len=67
        cfg.overlap_len=8
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - check if biggest fragment on output <= from requested fragment len. + overlap len.
        self.assertTrue(outf.length.max()<=cfg.fragment_len+cfg.overlap_len, "Fragment+overlap bigger than expected")

    def test_specialization(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_specialization.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.specialization="spec"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - check if all specializations are on output
        self.assertTrue(inf["specialization"].isin(mergedf["specialization"]).all() , "Specialization not reported")

    def test_pre_cluster(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_pre_cluster.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.bin_len=200
        cfg.pre_cluster="rank-2"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - pre-clustered can and should be bigger than bin length
        self.assertTrue(mergedf.groupby(["binid"]).sum()["length"].max()>cfg.bin_len , "Pre-clustering failed")
    
    def test_bin_exclusive(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_exclusive.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.bin_len=1300
        cfg.bin_exclusive="rank-2"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - even when possible (bin_len=1300), using bin_exclusive should split the inputs into more bins
        self.assertTrue(mergedf["binid"].max()>0 , "Bin-exclusive clustering failed")

    def test_bins(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bins.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
        cfg.bins=12
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific test - matching number of bins with chosen value (not always but works in this case)
        # 12 bins instead of 13 because 2 entries are already clustered in one bin on the update_file
        self.assertEqual(mergedf["binid"].max()+1, cfg.bins, "Number of bins do not match")

    def test_bin_len(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_bin_len300.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"

        # Extended bin len
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

        # Test with smaller bin len
        cfg.output_file=self.results_dir+"test_bin_len200.tsv"
        cfg.bin_len=200
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")
        # specific check - No new entry fits on bin 0
        self.assertEqual(sum(mergedf["binid"]==0), 2, "Update failed")
    
    def test_all(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_all.tsv"
        cfg.update_file=self.base_dir+"data/bins_LJ.tsv"
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
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")

if __name__ == '__main__':
    unittest.main()