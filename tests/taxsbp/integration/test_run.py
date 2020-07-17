import unittest
import taxsbp.taxsbp

class TestRun(unittest.TestCase):
    def test_run(self):
        """
        Test if taxsbp runs
        """
        ret = taxsbp.taxsbp.main("")

        # check if ran okay
        self.assertFalse(ret, "TaxSBP show help")
       
if __name__ == '__main__':
    unittest.main()