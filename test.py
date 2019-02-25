import TaxSBP
import unittest

class TestTaxSBP(unittest.TestCase):

    def test_help(self):
        self.assertFalse(TaxSBP.main(["cluster","-h"]))


if __name__ == '__main__':
    unittest.main()
    