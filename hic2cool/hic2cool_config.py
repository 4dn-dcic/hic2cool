"""
This module contains global hic/cooler attributes shared across other modules
in hic2cool
"""
from __future__ import (
    absolute_import,
    division,
    print_function,
    unicode_literals
)
import numpy as np

# Cooler metadata
COOLER_FORMAT = "HDF5::Cooler"
COOLER_FORMAT_VERSION = 3
MCOOL_FORMAT = "HDF5::MCOOL"
MCOOL_FORMAT_VERSION = 2
URL = "https://github.com/4dn-dcic/hic2cool"
CHROM_DTYPE = np.dtype('S32')
CHROMID_DTYPE = np.int32
CHROMSIZE_DTYPE = np.int32
COORD_DTYPE = np.int32
BIN_DTYPE = np.int64
COUNT_DTYPE = np.int32
CHROMOFFSET_DTYPE = np.int64
BIN1OFFSET_DTYPE = np.int64
NORM_DTYPE = np.float64
PIXEL_FIELDS = ('bin1_id', 'bin2_id', 'count')
CHUNK_DTYPE = {'names':['bin1_id','bin2_id','count'], 'formats':[BIN_DTYPE, BIN_DTYPE, COUNT_DTYPE]}
H5OPTS = {'compression': 'gzip', 'compression_opts': 6, 'shuffle': True}
