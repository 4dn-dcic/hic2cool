#!/usr/bin/env python
# Test code for hic2cool
# Requires pip installation of cooler to run
# Must be run from this directory (/hic2cool/test/)


import unittest
from hic2cool import hic2cool_convert
from hic2cool import __version__
import cooler
import os
import warnings


class TestRunHic(unittest.TestCase):
    infile_name = 'test_data/test_hic.hic'
    outfile_name = 'test_data/test_cool.cool'
    binsize = 250000

    def test_run(self):
        with warnings.catch_warnings(record=True) as w:
            warnings.simplefilter("always")
            # this should trigger a warning, because test file is missing chrMT
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize)
            # verify some things about the warning
            assert len(w) == 49
            assert issubclass(w[-1].category, UserWarning)
        self.assertTrue(os.path.isfile(self.outfile_name))


class TestWithCooler(unittest.TestCase):
    outfile_name = 'test_data/test_cool.cool'
    binsize = 250000

    def test_cooler(self):
        cool = cooler.Cooler(self.outfile_name)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name, cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 8)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 25)
        self.assertEqual(self.binsize, cool.info['bin-size'])
        matrix_res = cool.matrix().fetch('chr1:25000000-25250000')
        self.assertEqual(matrix_res.shape, (1,1))
        self.assertEqual(matrix_res[0][0], 4.0)


if __name__ == '__main__':
    unittest.main()
