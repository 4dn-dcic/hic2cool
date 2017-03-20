#!/usr/bin/env python
# Test code for hic2cool
# Requires pip installation of cooler to run
# Must be run from this directory (./test)

# temp
import sys
sys.path.insert(0, '/Users/carl/Documents/GitHub/hic2cool/hic2cool')

import unittest
import hic2cool
import cooler
import os
import codecs
import numpy as np


class TestRunHic(unittest.TestCase):
    infile_name = 'test_data/test_hic.hic'
    outfile_name = 'test_data/test_cool.cool'
    binsize = 250000

    def test_run(self):
        # can take 10 to 15 mins to run locally
        hic2cool.hic2cool('KR', self.infile_name, 'BP', self.binsize, self.outfile_name)
        self.assertTrue(os.path.isfile(self.outfile_name))

class TestWithCooler(unittest.TestCase):
    outfile_name = 'test_data/test_cool.cool'
    binsize = 250000
    hic2cool_version = open('../hic2cool/_version.py').readlines()[-1].split()[-1].strip("\"'")

    def test_cooler(self):
        cool = cooler.Cooler(self.outfile_name)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name, cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 8)
        self.assertTrue(self.hic2cool_version in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 25)
        self.assertEqual(self.binsize, cool.info['bin-size'])
        matrix_res = cool.matrix().fetch('chr17:25260000-25265000')
        self.assertEqual(matrix_res.shape, (1,1))
        self.assertEqual(matrix_res[0][0], 79.0)


if __name__ == '__main__':
    unittest.main()
