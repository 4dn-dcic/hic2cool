#!/usr/bin/env python
# Test code for hic2cool
# Requires pip installation of cooler to run
# Must be run from this directory (/hic2cool/test/)

from __future__ import absolute_import, print_function, unicode_literals
import unittest
# will only work with pip installed package
from hic2cool import hic2cool_convert, __version__
import cooler
import os
import h5py
import sys
import subprocess
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
    outfile_name = 'test/test_data/test_cool_100000.cool'
    outfile_name2 = 'test/test_data/test_cool_2500000.cool'
    outfile_name_all = 'test/test_data/test_cool_multi_res.multi.cool'
    binsize = 100000
    binsize2 = 2500000


    def test_run_with_warnings(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize, True)
        read_out = out.getvalue().strip()
        self.assertTrue('WARNING' in read_out)


    def test_run_exclude_missing_100000(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize)
        read_out = out.getvalue().strip()
        self.assertFalse('WARNING' in read_out)
        self.assertTrue(os.path.isfile(self.outfile_name))


    def test_run_exclude_missings_2500000(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name2, self.binsize2)
        read_out = out.getvalue().strip()
        self.assertFalse('WARNING' in read_out)
        self.assertTrue(os.path.isfile(self.outfile_name2))


    def test_run_exclude_missing_multi_res(self):
        # run hic2cool for all resolutions in the hic file
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            hic2cool_convert(self.infile_name, self.outfile_name_all, 0)
        read_out = out.getvalue().strip()
        self.assertFalse('WARNING' in read_out)
        self.assertTrue(os.path.isfile(self.outfile_name_all))


class TestWithCooler(unittest.TestCase):
    outfile_name = 'test/test_data/test_cool_100000.cool'
    outfile_name2 = 'test/test_data/test_cool_2500000.cool'
    outfile_name_all = 'test/test_data/test_cool_multi_res.multi.cool'
    binsize = 100000
    binsize2 = 2500000


    def test_cooler_100000(self):
        h5file = h5py.File(self.outfile_name, 'r')
        cool = cooler.Cooler(h5file)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name.encode('utf-8'), cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 10)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(cool.info['nchroms'], 25)
        self.assertEqual(len(cool.chromnames), 25) # 'all' excluded
        self.assertEqual(self.binsize, cool.info['bin-size'])
        matrix_res = cool.matrix(balance=False).fetch('chr1:25000000-25100000')
        self.assertEqual(matrix_res.shape, (1,1))
        self.assertEqual(matrix_res[0][0], 2)
        # expect a ValueError, since we haven't balanced yet
        # currently this is not converging
        with self.assertRaises(ValueError):
            matrix_res = cool.matrix(balance=True).fetch('chr1:25000000-25100000')
        # subprocess.call(['cooler', 'balance', '-p', '10', self.outfile_name])
        # # re-open balanced cool
        # h5file = h5py.File(self.outfile_name, 'r')
        # cool = cooler.Cooler(h5file)
        # matrix_res = cool.matrix(balance=True).fetch('chr1:25000000-25250000')
        # self.assertEqual(matrix_res.shape, (1,1))
        # self.assertEqual(round(matrix_res[0][0],3), 3.043)


    def test_cooler_2500000(self):
        h5file = h5py.File(self.outfile_name2, 'a')
        cool = cooler.Cooler(h5file)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name2.encode('utf-8'), cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 10)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 25)
        self.assertEqual(self.binsize2, cool.info['bin-size'])
        matrix_res = cool.matrix(balance=False).fetch('chr1:0-25000000')
        self.assertEqual(matrix_res.shape, (10,10))
        self.assertEqual(matrix_res[9][9], 40)
        with self.assertRaises(ValueError):
            matrix_res = cool.matrix(balance=True).fetch('chr1:0-25000000')
        subprocess.call(['cooler', 'balance', self.outfile_name2])
        h5file.close()
        # re-open balanced cool
        h5file = h5py.File(self.outfile_name2, 'a')
        cool = cooler.Cooler(h5file)
        matrix_res = cool.matrix(balance=True).fetch('chr1:0-25000000')
        self.assertEqual(matrix_res.shape, (10,10))
        self.assertEqual(round(matrix_res[9][9],3), 0.613)


    def test_cooler_multi_res(self):
        h5file = h5py.File(self.outfile_name_all, 'r')
        # since this is multi-res, hdf5 structure is different
        # expect the following 9 resultions to be present:
        # [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]
        self.assertEqual(len(h5file['resolutions'].keys()), 9)
        # take a resolution and check the metadata
        res250kb = h5file['resolutions']['250000']
        cool = cooler.Cooler(res250kb)
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name_all.encode('utf-8'), cool_file)
        # cooler info has 8 entries
        self.assertEqual(len(cool.info), 10)
        self.assertTrue(__version__ in cool.info['generated-by'])
        self.assertEqual(len(cool.chromnames), 25)
        self.assertEqual(250000, cool.info['bin-size'])
        matrix_res = cool.matrix(balance=False).fetch('chr1:25000000-25250000')
        self.assertEqual(matrix_res.shape, (1,1))
        self.assertEqual(matrix_res[0][0], 4)
        # make sure the updated cooler multi-res syntax is working
        cool = cooler.Cooler(self.outfile_name_all + '::resolutions/100000')
        cool_file = cool.filename.encode('utf-8')
        self.assertEqual(self.outfile_name_all.encode('utf-8'), cool_file)
        cool_res = int(cool.info['bin-size'])
        self.assertEqual(100000, cool_res)


    def test_check_norms(self):
        NORMS = ["VC", "VC_SQRT", "KR"]
        h5file = h5py.File(self.outfile_name, 'r')
        bins = h5file['bins']
        for norm in NORMS:
            assert norm in bins.keys()


if __name__ == '__main__':
    unittest.main()
