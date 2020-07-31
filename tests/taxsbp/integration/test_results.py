import unittest
import taxsbp.taxsbp
import pandas as pd
import os, shutil, gzip
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
        decompress_gzip(self.base_dir+"data/20181219_abfv_refseq_cg.tsv.gz")
        decompress_gzip(self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp.gz")
        # Slice input for update tests
        head_file(self.base_dir+"data/20181219_abfv_refseq_cg.tsv", self.base_dir+"data/20181219_abfv_refseq_cg_first15k_lines.tsv", 15000)

    def test_real_data(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_real_data.tsv"
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg.tsv"
        cfg.nodes_file=self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp"
        
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)

        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
    
    def test_real_data_bin_exclusive(self):
        cfg = Config(**self.default_config)
        
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg.tsv"
        cfg.nodes_file=self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp"

        for bin_exclusive in ["superkingdom", "phylum","class","order","family","genus", "species", "taxid", "spec"]:
            if bin_exclusive=="spec":
                cfg.specialization="spec"
            cfg.output_file=self.results_dir+"test_real_data_bin_exclusive_"+bin_exclusive+".tsv"
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
    
    def test_real_data_pre_cluster(self):
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_real_data_pre_cluster.tsv"
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg.tsv"
        cfg.nodes_file=self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp"
        cfg.pre_cluster="superkingdom"
        cfg.bin_len=1

        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)

        # sanity check - fails because some sequences are invalidated
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        
        # specific check - should contain only one bin for each superkingdom (total 5) despite the short bin length
        self.assertEqual(outf.binid.max()+1, 5, "Input/Output files are inconsistent")
    
    def test_real_data_update(self):
        # Cluster first 15k lines and update with the rest
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_real_data_update.tsv"
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg_first15k_lines.tsv"
        cfg.nodes_file=self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp"
        cfg.bin_len=5000000
        cfg.fragment_len=999500
        cfg.overlap_len=500

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # Update - use full input, will ignore used lines
        cfg.update_file=cfg.output_file
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg.tsv"
        cfg.output_file=self.results_dir+"test_real_data_update2.tsv"
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)
        # Join update file on output
        mergedf = pd.concat([updf,outf], ignore_index=True)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, mergedf), "Input/Output files are inconsistent")

        # specific test - Check if bins did not pass bin_len
        max_bin_len = outf.groupby(["binid"]).sum()["length"].max()
        self.assertTrue(max_bin_len<=cfg.bin_len, "Bin is bigger than requested")
    
    def test_real_data_update_merge(self):
        # Cluster first 15k lines and update with the rest
        cfg = Config(**self.default_config)
        cfg.output_file=self.results_dir+"test_real_data_update_merge.tsv"
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg_first15k_lines.tsv"
        cfg.nodes_file=self.base_dir+"data/20181219_abfv_refseq_cg_nodes.dmp"
        cfg.bin_len=2000000
        cfg.fragment_len=999500
        cfg.overlap_len=500

        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")

        # Update - use full input, will ignore used lines
        # Increase bin_len so many bins will be merged
        cfg.update_file=cfg.output_file
        cfg.input_file=self.base_dir+"data/20181219_abfv_refseq_cg.tsv"
        cfg.output_file=self.results_dir+"test_real_data_update_merge2.tsv"
        cfg.bin_len=10000000
        cfg.allow_merge=True
 
        # run check
        self.assertTrue(taxsbp.taxsbp.pack(**vars(cfg)), "TaxSBP fails to run")
        inf, outf = parse_files(cfg)
        updf = parse_output(cfg.update_file)

        # sanity check
        self.assertTrue(sanity_check(cfg, inf, outf), "Input/Output files are inconsistent")
        # specific test - Check if bins did not pass bin_len
        max_bin_len = outf.groupby(["binid"]).sum()["length"].max()
        self.assertTrue(max_bin_len<=cfg.bin_len, "Bin is bigger than requested")

if __name__ == '__main__':
    unittest.main()