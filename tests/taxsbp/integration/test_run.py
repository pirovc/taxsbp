import unittest
import taxsbp.taxsbp

class TestRun(unittest.TestCase):
    def test_run(self):
        """
        Test if taxsbp runs
        """
        ret = taxsbp.taxsbp.pack(input_file="sample_data/seqinfo.tsv", nodes_file="sample_data/nodes.dmp")

        # check if ran okay
        self.assertTrue(ret, "TaxSBP finished successfuly")
       
if __name__ == '__main__':
    unittest.main()