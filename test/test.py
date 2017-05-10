#!/usr/bin/env python
# Test code for hic2cool
# Requires pip installation of cooler to run
# Must be run from this directory (/hic2cool/test/)

from __future__ import absolute_import, print_function
import unittest
# will only work with pip installed package
from hic2cool import hic2cool_convert, __version__
import cooler
import os
import h5py
import sys
from contextlib import contextmanager
try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO


# thanks to Rob Kennedy on S.O. for this bit of code
@contextmanager
def captured_output():
    new_out, new_err = StringIO(), StringIO()
    old_out, old_err = sys.stdout, sys.stderr
    try:
        sys.stdout, sys.stderr = new_out, new_err
        yield sys.stdout, sys.stderr
    finally:
        sys.stdout, sys.stderr = old_out, old_err


class TestRunHic(unittest.TestCase):
    infile_name = 'test/test_data/test_hic.hic'
    outfile_name = 'test/test_data/test_cool_250000.cool'
    outfile_name2 = 'test/test_data/test_cool_1000000.cool'
    binsize = 250000
    binsize2 = 1000000
    normalization = 'KR'

    def test_run_MT_error(self):
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            with self.assertRaises(SystemExit):
                hic2cool_convert(self.infile_name, self.outfile_name, self.binsize, self.normalization)
        read_err = err.getvalue().strip()
        self.assertTrue('ERROR' in read_err)

    def test_run_exclude_MT(self):
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize, self.normalization, True)
        read_err = err.getvalue().strip()
        self.assertFalse('ERROR' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_name))

    def test_run_exclude_MT_1000000(self):
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            hic2cool_convert(self.infile_name, self.outfile_name2, self.binsize2, self.normalization, True)
        read_err = err.getvalue().strip()
        self.assertFalse('ERROR' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_name))

class TestWithCooler(unittest.TestCase):
    outfile_name = 'test/test_data/test_cool_250000.cool'
    outfile_name2 = 'test/test_data/test_cool_1000000.cool'
    binsize = 250000
    binsize2 = 1000000

    def test_cooler_250000(self):
        h5file = h5py.File(self.outfile_name, 'r')
        res_data = h5file['resolutions/'+str(self.binsize)]
        cool = cooler.Cooler(res_data)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name.encode('utf-8'), cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 8)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 24) #MT is excluded
        self.assertEqual(self.binsize, cool.info['bin-size'])
        matrix_res = cool.matrix().fetch('chr1:25000000-25250000')
        self.assertEqual(matrix_res.shape, (1,1))
        self.assertEqual(matrix_res[0][0], 4.0)

    def test_cooler_1000000(self):
        h5file = h5py.File(self.outfile_name2, 'r')
        res_data = h5file['resolutions/'+str(self.binsize2)]
        cool = cooler.Cooler(res_data)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name2.encode('utf-8'), cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 8)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 24)
        self.assertEqual(self.binsize2, cool.info['bin-size'])
        matrix_res = cool.matrix().fetch('chr1:0-10000000')
        self.assertEqual(matrix_res.shape, (10,10))
        self.assertEqual(matrix_res[9][9], 15.0)


if __name__ == '__main__':
    unittest.main()
