#!/usr/bin/env python
# Test code for hic2cool
# Requires pip installation of cooler to run
# Must be run from the root directory

from __future__ import absolute_import, print_function, unicode_literals
import unittest
import cooler
import os
import h5py
import subprocess
import sys
import math
import hashlib
import shutil
import numpy as np
from hic2cool import (
    hic2cool_convert,
    hic2cool_update,
    hic2cool_extractnorms,
    hic2cool_print_stderr,
    hic2cool_force_exit,
    __version__
)
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


class TestRunConvertAndExtractNorms(unittest.TestCase):
    """
    Test methods are named "test_#_ ..." because unittest orders methods
    by sorted string name
    """
    infile_name = 'test_data/test_hic.hic'
    infile_no_norms = 'test_data/test_hic_no_norms.hic'
    outfile_name = 'test_data/test_cool_100000.cool'
    outfile_name2 = 'test_data/test_cool_2500000.cool'
    outfile_name_all = 'test_data/test_cool_multi_res.mcool'
    outfile_no_norms = 'test_data/test_cool_no_norms.mcool'
    binsize = 100000
    binsize2 = 2500000

    def get_norms_content_md5(self, fname, use_resolutions=None, bin_lens={}):
        """
        Return hexidecimal digest of md5 from all res in use_resolutions (list)
        Also return int bin length found in the file
        """
        str_val = ''
        with h5py.File(fname) as h5:
            if not use_resolutions:
                use_resolutions = [r for r in h5['resolutions']]
            for res in use_resolutions:
                res = str(res)
                if res not in h5['resolutions']:
                    continue
                for norm in ['KR', 'VC', 'VC_SQRT']:
                    if not bin_lens.get(res):
                        bin_lens[res] = len(h5['resolutions'][res]['bins'][norm])
                    str_val += ','.join(repr(item) for item in
                                        h5['resolutions'][res]['bins'][norm][:bin_lens[res]])
        str_val = str_val.encode()
        md5 = hashlib.md5()
        md5.update(str_val)
        return md5.hexdigest(), bin_lens

    def test_0_run_with_warnings(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize, True)
        read_err = err.getvalue().strip()
        self.assertTrue('WARNING' in read_err)

    def test_1_run_exclude_missing_100000(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name, self.binsize)
        read_err = err.getvalue().strip()
        self.assertFalse('WARNING' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_name))

    def test_2_run_exclude_missings_2500000(self):
        with captured_output() as (out, err):
            hic2cool_convert(self.infile_name, self.outfile_name2, self.binsize2)
        read_err = err.getvalue().strip()
        self.assertFalse('WARNING' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_name2))

    def test_3_run_exclude_missing_multi_res_no_norms(self):
        # run hic2cool for all resolutions in the hic file
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            hic2cool_convert(self.infile_no_norms, self.outfile_no_norms, 0)
        read_err = err.getvalue().strip()
        self.assertTrue('WARNING. No normalization vectors' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_no_norms))

    def test_4_run_exclude_missing_multi_res(self):
        # run hic2cool for all resolutions in the hic file
        with captured_output() as (out, err):
            # this should fail, because test file is missing chrMT
            # and excludeMT was not specified
            hic2cool_convert(self.infile_name, self.outfile_name_all, 0)
        read_err = err.getvalue().strip()
        self.assertFalse('WARNING' in read_err)
        self.assertTrue(os.path.isfile(self.outfile_name_all))

    def test_5_no_norms_in_hic(self):
        with captured_output() as (out, err):
            hic2cool_extractnorms(self.infile_no_norms, self.outfile_no_norms,
                                      exclude_mt=True, show_warnings=True)
        read_out = out.getvalue().strip()
        self.assertTrue('Normalizations:  []' in read_out)
        # no mt in this hic file to begin with
        self.assertTrue('No chromosome found when attempting to exclude MT' in read_out)

    def test_6_add_norms_to_outfile_no_norm(self):
        # infile: [2500000, 1000000, 500000, 250000, 100000, 50000, 25000, 10000, 5000]
        # outfile: [1000000, 16000000, 2000000, 4000000, 500000, 8000000]
        # so, only expected 1000000 and 500000 to be shared
        # copy outfile twice
        copy1_name = 'test_data/test_cool_no_norms_copy1.mcool'
        shutil.copy(self.outfile_no_norms, copy1_name)
        shared = [1000000, 500000]
        not_shared = [2500000, 250000, 100000, 50000, 25000, 10000, 5000]
        with captured_output() as (out, err):
            hic2cool_extractnorms(self.infile_name, copy1_name,
                                  exclude_mt=True, show_warnings=True)
        read_out = out.getvalue().strip()
        for share_res in shared:
            self.assertTrue('normalization vector at %s BP' % share_res in read_out)
        for not_share_res in not_shared:
            self.assertTrue('Skip resolution %s' % not_share_res in read_out)
        self.assertTrue('Excluding chromosome chrMT with index 25' in read_out)
        # ensure that running the function again results the same norms
        # run without exclude_mt, which is actually missing anyways
        init_md5, bin_lens = self.get_norms_content_md5(copy1_name, shared)
        with captured_output() as (out2, err2):
            hic2cool_extractnorms(self.infile_name, copy1_name, show_warnings=True)
        read_err = err2.getvalue().strip()
        self.assertTrue('Normalization vector KR does not exist for chrMT' in read_err)
        fin_md5, _ = self.get_norms_content_md5(copy1_name, shared, bin_lens)
        self.assertEqual(init_md5, fin_md5)
        # ensure that norm values match between different runs of extract-norms
        diff_md5, _ = self.get_norms_content_md5(self.outfile_name_all, shared, bin_lens)
        self.assertEqual(init_md5, diff_md5)
        # ensure that extracting norms again doesn't change the mult-res cooler
        # this time use all resolutions
        init2_md5, bin_lens2 = self.get_norms_content_md5(self.outfile_name_all)
        hic2cool_extractnorms(self.infile_name, self.outfile_name_all, silent=True)
        fin2_md5, bin_lens2_2 = self.get_norms_content_md5(self.outfile_name_all)
        self.assertEqual(init2_md5, fin2_md5)
        self.assertEqual(bin_lens2, bin_lens2_2)


class TestWithCooler(unittest.TestCase):
    outfile_name = 'test_data/test_cool_100000.cool'
    outfile_name2 = 'test_data/test_cool_2500000.cool'
    outfile_name_all = 'test_data/test_cool_multi_res.mcool'
    outfile_no_norms = 'test_data/test_cool_no_norms.mcool'
    binsize = 100000
    binsize2 = 2500000

    def test_cooler_100000(self):
        with h5py.File(self.outfile_name, 'r') as h5file:
            cool = cooler.Cooler(h5file)
            cool_file = cool.filename.encode('utf-8')
            self.assertEqual(self.outfile_name.encode('utf-8'), cool_file)
            # cooler info has 8 entries
            self.assertEqual(len(cool.info), 10)
            self.assertTrue(__version__ in cool.info['generated-by'])
            self.assertEqual(cool.info['nchroms'], 25)
            self.assertEqual(len(cool.chromnames), 25) # 'all' excluded
            self.assertEqual(self.binsize, cool.info['bin-size'])
            # check normalization values and balanced/unbalanced counts
            matrix_res = cool.matrix(balance=False).fetch('chr1:25000000-25100000')
            bin_idx = 250 # chr1:25000000-25100000 corresponds to bin 0 at this res
            bin_raw_val = matrix_res[0][0]
            self.assertEqual(matrix_res.shape, (1,1))
            self.assertEqual(bin_raw_val, 2)
            # expect a ValueError -- not balanced and there is no 'weights' column
            with self.assertRaises(ValueError):
                matrix_res = cool.matrix(balance=True).fetch('chr1:25000000-25100000')
            # check a few norms
            bin_info = cool.bins()[bin_idx]
            for norm in ['KR', 'VC', 'VC_SQRT']:
                norm = str(norm) # needed for python 2
                self.assertTrue(norm in bin_info)
                bin_norm_val = bin_info[norm][bin_idx]  # value for weight in this bin
                norm_matrix_res = cool.matrix(balance=norm).fetch('chr1:25000000-25100000')
                self.assertEqual(norm_matrix_res.shape, (1,1))
                # handle nan
                if math.isnan(bin_norm_val):
                    self.assertTrue(math.isnan(norm_matrix_res[0][0]))
                else:
                    # balanced value is equal to count * weight * weight
                    calc_balanced_val = round(bin_raw_val * bin_norm_val * bin_norm_val, 4)
                    found_balanced_val = round(norm_matrix_res[0][0], 4)
                    self.assertEqual(calc_balanced_val, found_balanced_val)
        # lastly check cooler dump raw count
        v = subprocess.check_output(['cooler', 'dump', self.outfile_name, '-t',
                                     'pixels', '-r', 'chr1:25000000-25100000'])
        formatted_v = [int(vl) for vl in v.decode().strip().split('\t')]
        # output corresponds to bin1, bin2, count
        self.assertEqual(formatted_v, [bin_idx, bin_idx, bin_raw_val])

    def test_cooler_2500000(self):
        with h5py.File(self.outfile_name2, 'a') as h5file:
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
            bin_idx = 9 # chr1:22500000-25000000 corresponds to bin 9 at this res
            bin_raw_val = matrix_res[9][9]
            self.assertEqual(bin_raw_val, 40)
            # expect a ValueError -- not balanced and there is no 'weights' column
            with self.assertRaises(ValueError):
                matrix_res = cool.matrix(balance=True).fetch('chr1:0-25000000')
            # check a few norms
            bin_info = cool.bins()[bin_idx]
            for norm in ['KR', 'VC', 'VC_SQRT']:
                norm = str(norm) # needed for python 2
                self.assertTrue(norm in bin_info)
                bin_norm_val = bin_info[norm][bin_idx]  # value for weight in this bin
                norm_matrix_res = cool.matrix(balance=norm).fetch('chr1:0-25000000')
                self.assertEqual(norm_matrix_res.shape, (10,10))
                # handle nan specially
                if math.isnan(bin_norm_val):
                    self.assertTrue(math.isnan(norm_matrix_res[9][9]))
                else:
                    # balanced value is equal to count * weight * weight
                    calc_balanced_val = round(bin_raw_val * bin_norm_val * bin_norm_val, 4)
                    found_balanced_val = round(norm_matrix_res[9][9], 4)
                    self.assertEqual(calc_balanced_val, found_balanced_val)
        # lastly check cooler dump raw count
        v = subprocess.check_output(['cooler', 'dump', self.outfile_name2, '-t',
                                     'pixels', '-r', 'chr1:22500000-25000000'])
        formatted_v = [int(vl) for vl in v.decode().strip().split('\t')]
        # output corresponds to bin1, bin2, count
        self.assertEqual(formatted_v, [bin_idx, bin_idx, bin_raw_val])

    def test_cooler_multi_res(self):
        with h5py.File(self.outfile_name_all, 'r') as h5file:
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
        with h5py.File(self.outfile_name, 'r') as h5file:
            bins = h5file['bins']
            for norm in NORMS:
                self.assertTrue(norm in bins.keys())

    def test_no_norms(self):
        """
        Added support for missing norm vectors in the hic file
        outfile_no_norms is a mcool
        """
        NORMS = ["VC", "VC_SQRT", "KR"]
        with h5py.File(self.outfile_no_norms, 'r') as h5file:
            # resolutions are ['1000000', '16000000', '2000000', '4000000', '500000', '8000000']
            self.assertEqual(len(h5file['resolutions'].keys()), 6)
            res500kb_bins = h5file['resolutions']['500000']['bins']
            for norm in NORMS:
                self.assertTrue(norm not in res500kb_bins.keys())


class TestRunUpdate(unittest.TestCase):
    infile_name = 'test_data/old_version_single_res.cool'
    outfile_name = 'test_data/new_version_single_res.cool'

    def norm_convert(self, val):
        if val != 0.0:
            return 1 / val
        else:
            return np.nan

    def test_run_update(self):
        """
        Test to ensure the update functionality works and changes generated
        by tags. This specific update will re-invert normalization vectors,
        so ensure that also works.
        """
        # first check some stuff on the input file, which is single-res
        original_kr_data = []  # first 10 values in KR normalization
        final_kr_data = []
        original_creation_date = ''
        update_date = ''
        with h5py.File(self.infile_name, 'r') as h5file:
            self.assertEqual(h5file.attrs.get('generated-by'), 'hic2cool-0.4.2')
            self.assertEqual(h5file.attrs.get('update-date'), None)
            original_creation_date = h5file.attrs.get('creation-date')
            original_kr_data = h5file['bins/KR'][:100]
        hic2cool_update(self.infile_name, self.outfile_name, silent=True)
        self.assertTrue(os.path.isfile(self.outfile_name))
        # ensure that the new file has a new version and update-date
        expected_version = 'hic2cool-' + __version__
        with h5py.File(self.outfile_name, 'r') as h5file:
            self.assertEqual(h5file.attrs.get('generated-by'), expected_version)
            self.assertEqual(h5file.attrs.get('creation-date'), original_creation_date)
            self.assertTrue(h5file.attrs.get('update-date') is not None)
            update_date = h5file.attrs.get('update-date')
            final_kr_data = h5file['bins/KR'][:100]
        # make sure the norms got updated correctly
        for comp in zip(original_kr_data, final_kr_data):
            if math.isnan(comp[0]):  # special case
                self.assertTrue(math.isnan(comp[1]))
                self.assertTrue(math.isnan(self.norm_convert(comp[1])))
            else:
                self.assertEqual(comp[0], self.norm_convert(comp[1]))
        # make sure running it again does nothing (update-date unchanged)
        hic2cool_update(self.outfile_name, silent=True)
        with h5py.File(self.outfile_name, 'r') as h5file:
            self.assertEqual(h5file.attrs.get('generated-by'), expected_version)
            self.assertEqual(h5file.attrs.get('update-date'), update_date)


class TestUtilities(unittest.TestCase):
    infile_name = 'test_data/test_hic.hic'

    def test_print_stderr(self):
        with captured_output() as (out, err):
            hic2cool_print_stderr('error message!')
        read_err = err.getvalue().strip()
        self.assertTrue('error message!' in read_err)

    def test_force_exit(self):
        """
        Open a dummy file and make sure it is closed
        """
        req = open(self.infile_name, 'rb')
        with captured_output() as (out, err):
            with self.assertRaises(SystemExit) as exc:
                hic2cool_force_exit('fatal error!', req)
            self.assertEqual(exc.exception.code, 1)
        read_err = err.getvalue().strip()
        self.assertTrue('fatal error!' in read_err)
        self.assertTrue(req.closed)  # file closed by force_exit


if __name__ == '__main__':
    unittest.main()
