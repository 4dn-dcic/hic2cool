"""
hic2cool
--------

CLI for a converter between .hic files (from juicer) and .cool files (for
cooler).  Both hic and cool files describe Hi-C contact matrices. Intended to be
lightweight, this is a stand-alone Python file written by Carl Vitzthum of the
HMS DBMI Park lab. Originally published 1/26/17.

This code was originally based off of the straw project by Neva C. Durand and
Yue Wu (https://github.com/theaidenlab/straw). The cooler file writing was based
off of much of the CLI code contained in this repo:
https://github.com/mirnylab/cooler.

The most import method is hic2cool_convert, which can be imported if the package
is pip installed.

Usage on the command line is: python hic2cool.py infile outfile resolution
normalization include_MT If an invalid resolution is given (i.e. one that does
not correspond to a hic resolution, the program will terminate and prompt you to
enter a valid one) See the README or main() function in this file for more usage
information.

"""
from __future__ import absolute_import, division, print_function, unicode_literals
from collections import OrderedDict
import sys
import time
import struct
import zlib
import math
from datetime import datetime

import numpy as np
import h5py
from ._version import __version__

# Global hic normalization types used
NORMS = []
# determine size of counts chunk
CHUNK_SIZE = 1000
# Cooler metadata
MAGIC = "HDF5::Cooler"
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


# read function
def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return buf.encode("utf-8", errors="ignore")
            return buf.decode("utf-8")
        elif b == "":
            raise EOFError("Buffer unexpectedly empty while trying to read null-terminated string")
        else:
            buf += b


def read_header(infile):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms
    contained in file) and values are [chr idx (int), chr name (str), chrom
    length (str)]. Returns the masterindex used by the file as well as the open
    file object.

    """
    req = open(infile, 'rb')
    chrs = {}
    resolutions = []
    magic_string = struct.unpack(b'<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        print('This does not appear to be a HiC file; '
              'magic string is incorrect')
        sys.exit()
    global version
    version = struct.unpack(b'<i', req.read(4))[0]
    masterindex = struct.unpack(b'<q', req.read(8))[0]
    genome = b""
    c = req.read(1)
    while (c != b'\0'):
        genome += c
        c = req.read(1)
    genome = genome.decode('ascii')
    nattributes = struct.unpack(b'<i', req.read(4))[0]
    for x in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
    nChrs = struct.unpack(b'<i', req.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(req)
        length = struct.unpack(b'<i', req.read(4))[0]
        if name and length:
            formatted_name = ('chr' + name if ('all' not in name.lower() and
                              'chr' not in name.lower()) else name.lower())
            formatted_name = ('chrM' if formatted_name == 'chrMT'
                              else formatted_name)
            chrs[i] = [i, formatted_name, length]
    nBpRes = struct.unpack(b'<i', req.read(4))[0]
    # find bp delimited resolutions supported by the hic file
    for x in range(0, nBpRes):
        res = struct.unpack(b'<i', req.read(4))[0]
        resolutions.append(res)
    return req, chrs, resolutions, masterindex, genome


def read_footer(req, master, unit, resolution):
    """
    Takes in an open hic file and generates two dictionaries. pair_footer_info
    contains the file position of info for any given chromosome pair (formatted
    as a string). Chr_footer_info gives chromosome-level size and position info
    relative to the file for each normalization type. This way, this function
    only has to run once. This is code from straw.

    """
    pair_footer_info = {}
    chr_footer_info = {}
    req.seek(master)
    nBytes = struct.unpack(b'<i', req.read(4))[0]
    nEntries = struct.unpack(b'<i', req.read(4))[0]
    found = False
    for i in range(nEntries):
        stri = readcstr(req)
        fpos = struct.unpack(b'<q', req.read(8))[0]
        sizeinbytes = struct.unpack(b'<i', req.read(4))[0]
        pair_footer_info[stri] = fpos
    nExpectedValues = struct.unpack(b'<i', req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        binSize = struct.unpack(b'<i', req.read(4))[0]
        nValues = struct.unpack(b'<i', req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack(b'<d', req.read(8))[0]
        nNormalizationFactors = struct.unpack(b'<i', req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack(b'<i', req.read(4))[0]
            v = struct.unpack(b'<d', req.read(8))[0]
    nExpectedValues = struct.unpack(b'<i', req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        str_ = readcstr(req)
        binSize = struct.unpack(b'<i', req.read(4))[0]
        nValues = struct.unpack(b'<i', req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack(b'<d', req.read(8))[0]
        nNormalizationFactors = struct.unpack(b'<i', req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack(b'<i', req.read(4))[0]
            v = struct.unpack(b'<d', req.read(8))[0]
    nEntries = struct.unpack(b'<i', req.read(4))[0]
    for i in range(nEntries):
        normtype = readcstr(req)
        if normtype not in NORMS:
            NORMS.append(normtype)
        if normtype not in chr_footer_info:
            chr_footer_info[normtype] = {}
        chrIdx = struct.unpack(b'<i', req.read(4))[0]
        unit1 = readcstr(req)
        resolution1 = struct.unpack(b'<i', req.read(4))[0]
        filePosition = struct.unpack(b'<q', req.read(8))[0]
        sizeInBytes = struct.unpack(b'<i', req.read(4))[0]
        if (unit1 == unit and resolution1 == resolution):
            chr_footer_info[normtype][chrIdx] = {
                'position': filePosition,
                'size': sizeInBytes
            }

    return req, pair_footer_info, chr_footer_info


# FUN(fin, unit, resolution, blockMap)
# Return [storeBlockData, myBlockBinCount, myBlockColumnCount]
def read_matrix_zoom_data(req, myunit, mybinsize, blockMap):
    unit = readcstr(req)
    temp = struct.unpack(b'<i', req.read(4))[0]
    temp2 = struct.unpack(b'<f', req.read(4))[0]
    temp2 = struct.unpack(b'<f', req.read(4))[0]
    temp2 = struct.unpack(b'<f', req.read(4))[0]
    temp2 = struct.unpack(b'<f', req.read(4))[0]
    binSize = struct.unpack(b'<i', req.read(4))[0]
    blockBinCount = struct.unpack(b'<i', req.read(4))[0]
    blockColumnCount = struct.unpack(b'<i', req.read(4))[0]
    storeBlockData = False
    # for the initial
    myBlockBinCount = -1
    myBlockColumnCount = -1
    if (myunit == unit and mybinsize == binSize):
        myBlockBinCount = blockBinCount
        myBlockColumnCount = blockColumnCount
        storeBlockData = True
    nBlocks = struct.unpack(b'<i', req.read(4))[0]
    for b in range(nBlocks):
        blockNumber = struct.unpack(b'<i', req.read(4))[0]
        filePosition = struct.unpack(b'<q', req.read(8))[0]
        blockSizeInBytes = struct.unpack(b'<i', req.read(4))[0]
        entry = dict()
        entry['size'] = blockSizeInBytes
        entry['position'] = filePosition
        if (storeBlockData):
            blockMap[blockNumber] = entry
    return [storeBlockData, myBlockBinCount, myBlockColumnCount]


# FUN(fin, myFilePos, unit, binsize, blockMap)
# Return [blockBinCount, blockColumnCount]
def read_matrix(req, myFilePos, unit, binsize, blockMap):
    req.seek(myFilePos)
    c1 = struct.unpack(b'<i', req.read(4))[0]
    c2 = struct.unpack(b'<i', req.read(4))[0]
    nRes = struct.unpack(b'<i', req.read(4))[0]
    i = 0
    found = False
    blockBinCount = -1
    blockColumnCount = -1
    while (i < nRes and (not found)):
        list1 = read_matrix_zoom_data(req, unit, binsize, blockMap)
        found = list1[0]
        if (list1[1] != -1 and list1[2] != -1):
            blockBinCount = list1[1]
            blockColumnCount = list1[2]
        i = i + 1
    if (not found):
        print("Error finding block data\n")
        return -1
    return [blockBinCount, blockColumnCount]


# FUN(regionIndices, blockBinCount, blockColumnCount, intra)
# Return blocksSet
def get_block_numbers_for_region_from_bin_position(
        regionIndices, blockBinCount, blockColumnCount, intra):
    col1 = int(regionIndices[0] / blockBinCount)
    col2 = int((regionIndices[1] + 1) / blockBinCount)
    row1 = int(regionIndices[2] / blockBinCount)
    row2 = int((regionIndices[3] + 1) / blockBinCount)
    blocksSet = set()
    for r in range(row1, row2 + 1):
        for c in range(col1, col2 + 1):
            blockNumber = r * blockColumnCount + c
            blocksSet.add(blockNumber)
    if intra:
        for r in range(col1, col2 + 1):
            for c in range(row1, row2 + 1):
                blockNumber = r * blockColumnCount + c
                blocksSet.add(blockNumber)
    return blocksSet


# FUN(fin, blockNumber, blockMap)
def read_block(req, blockNumber, blockMap):
    idx = dict()
    if(blockNumber in blockMap):
        idx = blockMap[blockNumber]
    else:
        idx['size'] = 0
        idx['position'] = 0
    if (idx['size'] == 0):
        return []
    req.seek(idx['position'])
    compressedBytes = req.read(idx['size'])
    uncompressedBytes = zlib.decompress(compressedBytes)
    nRecords = struct.unpack(b'<i', uncompressedBytes[0:4])[0]
    v = []
    global version
    if (version < 7):
        for i in range(nRecords):
            binX = struct.unpack(b'<i', uncompressedBytes[(12*i):(12*i+4)])[0]
            binY = struct.unpack(b'<i', uncompressedBytes[(12*i+4):(12*i+8)])[0]
            counts = struct.unpack(b'<f', uncompressedBytes[(12*i+8):(12*i+12)])[0]
            record = dict()
            record['binX'] = binX
            record['binY'] = binY
            record['counts'] = counts
            v.append(record)
    else:
        binXOffset = struct.unpack(b'<i', uncompressedBytes[4:8])[0]
        binYOffset = struct.unpack(b'<i', uncompressedBytes[8:12])[0]
        useShort = struct.unpack(b'<b', uncompressedBytes[12:13])[0]
        type_ = struct.unpack(b'<b', uncompressedBytes[13:14])[0]
        index = 0
        if (type_ == 1):
            rowCount = struct.unpack(b'<h', uncompressedBytes[14:16])[0]
            temp = 16
            for i in range(rowCount):
                y = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
                temp = temp+2
                binY = y + binYOffset
                colCount = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
                temp = temp+2
                for j in range(colCount):
                    x = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
                    temp = temp+2
                    binX = binXOffset + x
                    if (useShort==0):
                        c = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
                        temp = temp+2
                        counts = c
                    else:
                        counts = struct.unpack(b'<f', uncompressedBytes[temp:(temp+4)])[0]
                        temp = temp+4
                    record = dict()
                    record['binX'] = binX
                    record['binY'] = binY
                    record['counts'] = counts
                    v.append(record)
                    index = index + 1
        elif (type_== 2):
            temp = 14
            nPts = struct.unpack(b'<i', uncompressedBytes[temp:(temp+4)])[0]
            temp = temp+4
            w = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
            temp = temp+2
            for i in range(nPts):
                row = int(i / w)
                col = i - row * w
                bin1 = int(binXOffset + col)
                bin2 = int(binYOffset + row)
                if (useShort == 0):
                    c = struct.unpack(b'<h', uncompressedBytes[temp:(temp+2)])[0]
                    temp = temp+2
                    if (c != -32768):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = c
                        v.append(record)
                        index = index + 1
                else:
                    counts = struct.unpack(b'<f',uncompressedBytes[temp:(temp+4)])[0]
                    temp=temp+4
                    if (counts != 0x7fc00000):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = counts
                        v.append(record)
                        index = index + 1
    return v


# FUN(fin, entry)
# Return Norm
def read_normalization_vector(req, entry):
    req.seek(entry['position'])
    nValues = struct.unpack(b'<i', req.read(4))[0]
    value = []
    for i in range(nValues):
        d = struct.unpack(b'<d', req.read(8))[0]
        value.append(d)
    return value


def parse_hic(req, outfile, chr1, chr2, unit, binsize, covered_chr_pairs,
              pair_footer_info, chr_offset_map, low_mem):
    """
    Adapted from the straw() function in the original straw package.
    Mainly, since all chroms are iterated over, the read_header and read_footer
    functions were placed outside of straw() and made to be reusable across
    any chromosome pair.
    As of version 0.4.0, counts and normalization vectors are written to the
    output cool files in chunks to save memory.
    As of version 0.3.6, all normalization vectors are automatically returned
    to be written in the output cool file.
    """
    blockMap = {}
    magic_string = ""
    if unit not in ["BP", "FRAG"]:
        print("Unit specified incorrectly, must be one of <BP/FRAG>")
        force_exit(warn_string, req)
    chr1ind = chr1[0]
    chr2ind = chr2[0]
    c1pos1 = 0
    c1pos2 = int(chr1[2])
    c2pos1 = 0
    c2pos2 = int(chr2[2])
    c1 = min(chr1ind, chr2ind)
    c2 = max(chr1ind, chr2ind)
    chr_key = str(c1) + "_" + str(c2)
    origRegionIndices = []
    regionIndices = []
    if (chr1ind > chr2ind):
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        regionIndices.append(int(c2pos1 / binsize))
        regionIndices.append(int(c2pos2 / binsize))
        regionIndices.append(int(c1pos1 / binsize))
        regionIndices.append(int(c1pos2 / binsize))
    else:
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        regionIndices.append(int(c1pos1 / binsize))
        regionIndices.append(int(c1pos2 / binsize))
        regionIndices.append(int(c2pos1 / binsize))
        regionIndices.append(int(c2pos2 / binsize))
    try:
        pair_footer_info[chr_key]
    except KeyError:
        warn_string = (
            'ERROR. There is a discrepancy between the chrs declared in the ' +
            'infile header and the actual information it contains.\nThe ' +
            'intersection between ' + chr1[1] + ' and ' + chr2[1] + ' could ' +
            'not be found in the file.')
        force_exit(warn_string, req)
    myFilePos = pair_footer_info[chr_key]
    list1 = read_matrix(req, myFilePos, unit, binsize, blockMap)
    blockBinCount = list1[0]
    blockColumnCount = list1[1]
    blockNumbers = get_block_numbers_for_region_from_bin_position(
        regionIndices, blockBinCount, blockColumnCount, c1 == c2)
    # join chunk is used in the non low_mem version. it is an np array
    # that will contain the entire list of c1/c2 counts
    join_chunk = np.zeros(shape=0, dtype=CHUNK_DTYPE) if not low_mem else None
    for chunk in generate_counts_chunk(req, c1, c2, blockNumbers, blockMap,
                                origRegionIndices, chr_offset_map):
        if low_mem:
            # write this chunk to a temporary dataset
            write_chunk(outfile, binsize, chunk)
        else:
            join_chunk = build_chunk(join_chunk, chunk)
    covered_chr_pairs.append(chr_key)
    return join_chunk


def generate_counts_chunk(req, c1, c2, blockNumbers, blockMap, origRegionIndices, chr_offset_map):
    """
    Generator function that assembles chunks of size CHUNK_SIZE, which are
    2d arrays of size CHUNK_SIZE x 4, in order [bin1, bin2, count, norm values]
    """
    counts_chunk = []
    for i_set in (blockNumbers):
        records = read_block(req, i_set, blockMap)
        for j in range(len(records)):
            rec = records[j]
            # x and y are bin indices
            x = rec['binX']
            y = rec['binY']
            c = rec['counts']
            # bins with offsets considered
            bin1 = chr_offset_map[c1] + x
            bin2 = chr_offset_map[c2] + y
            if ((x >= origRegionIndices[0] and x <= origRegionIndices[1] and
                 y >= origRegionIndices[2] and y <= origRegionIndices[3]) or
                 ((c1 == c2) and
                   y >= origRegionIndices[0] and y <= origRegionIndices[1] and
                   x >= origRegionIndices[2] and x <= origRegionIndices[3])):
                counts_chunk.append([bin1, bin2, c])
            if len(counts_chunk) >= CHUNK_SIZE:
                yield chunk_to_array(counts_chunk)
                counts_chunk = []
    yield chunk_to_array(counts_chunk)


def chunk_to_array(counts_chunk):
    """
    Takes the 2d array structure of the chunks and converts to a named
    2d numpy array with the correct datatypes. This step is important
    because the structured array is used to sort the datasets in the pixels
    group in the finalize_resolution_cool function.
    """
    if len(counts_chunk) == 0:
        return np.zeros(shape=0, dtype=CHUNK_DTYPE)
    intermediate = np.array(counts_chunk)
    named_arr = np.zeros(len(counts_chunk), dtype=CHUNK_DTYPE)
    for idx, name in enumerate(CHUNK_DTYPE['names']):
        named_arr[name] = intermediate[:,idx]
    # the normalization values
    return named_arr


def build_chunk(total_chunk, counts_chunk):
    """
    The non hdf5-storing method of extending the chunks that will be sorted
    and written to the pixels table in the prepare_chunks_dataset fxn
    """
    this_chunk_size = len(counts_chunk)
    if this_chunk_size == 0:
        return total_chunk
    prev_size = total_chunk.shape[0]
    total_chunk = np.resize(total_chunk, prev_size + this_chunk_size)
    total_chunk[prev_size:] = counts_chunk
    return total_chunk


def write_chunk(outfile, resolution, counts_chunk):
    """
    Write a counts chunk generated by generate_counts_chunk function
    to the temporary chunks dataset for the given resolution. This temporary
    dataset is used to generate the final 'pixels' table down the line.
    """
    this_chunk_size = len(counts_chunk)
    if this_chunk_size == 0:
        return
    with h5py.File(outfile, "r+") as h5file:
        h5resolution = h5file[str(resolution)] if multi_res else h5file
        dataset = h5resolution['chunks']
        # dsets indices line up with counts_chunk indices
        prev_size = dataset.shape[0]
        dataset.resize(prev_size + this_chunk_size, axis=0)
        dataset[prev_size:] = counts_chunk


def initialize_cool(outfile, req, chr_info, genome, resolution, chr_footer_info, low_mem):
    """
    Use various information to initialize a cooler file in HDF5 format. Tables
    included are chroms, bins, pixels, and indexes.
    For an in-depth explanation of the data structure, please see:
    http://cooler.readthedocs.io/en/latest/datamodel.html
    """
    h5opts = dict(compression='gzip', compression_opts=6)
    with h5py.File(outfile, "a") as h5file:
        if str(resolution) in h5file.keys():
            del h5file[str(resolution)]
        if multi_res:
            h5resolution = h5file.create_group(str(resolution))
        else:
            h5resolution = h5file
        for key in h5resolution.keys():
            del h5resolution[key]
        # chroms
        grp = h5resolution.create_group('chroms')
        chr_names = write_chroms(grp, chr_info, h5opts)
        # bins
        bin_table, chr_offsets, chr_n_bins, by_chr_offset_map = create_bins(chr_info, resolution)
        grp = h5resolution.create_group('bins')
        write_bins(grp, req, chr_names, bin_table, chr_n_bins, chr_footer_info, h5opts)
        # chunks (a temporary group)
        n_bins = len(bin_table)
        if low_mem:
            prepare_chunks_dataset(h5resolution, n_bins, h5opts)
        # indexes (just initialize bin1offets)
        grp = h5resolution.create_group('indexes')
        write_chrom_offset(grp, chr_offsets, h5opts)
        # pixels (written later)
        grp = h5resolution.create_group('pixels')
        initialize_pixels(grp, n_bins, h5opts)
        # w.r.t. metadata, each resolution is considered a full file
        info = {}
        info['nchroms'] = len(chr_names)
        info['nbins'] = n_bins
        info['bin-type'] = 'fixed'
        info['bin-size'] = resolution
        info['format'] = MAGIC
        info['format-url'] = URL
        info['generated-by'] = 'hic2cool-' + __version__
        info['genome-assembly'] = genome
        info['creation-date'] = datetime.now().isoformat()
        h5resolution.attrs.update(info)
    return by_chr_offset_map


def write_chroms(grp, chrs, h5opts):
    """
    Write the chroms table, which includes chromosome names and length
    """
    # format chr_names and chr_lengths np arrays
    chr_names = np.array(
        [val[1] for val in chrs.values() if val[1] != 'all'],
        dtype=CHROM_DTYPE)
    chr_lengths = np.array(
        [val[2] for val in chrs.values() if val[1] != 'all'])
    grp.create_dataset(
        'name',
        shape=(len(chr_names),),
        dtype=CHROM_DTYPE,
        data=chr_names,
        **h5opts)
    grp.create_dataset(
        'length',
        shape=(len(chr_lengths),),
        dtype=CHROMSIZE_DTYPE,
        data=chr_lengths,
        **h5opts)
    return chr_names


def create_bins(chrs, binsize):
    """
    Cooler requires a bin file with each line having <chrId, binStart, binEnd>
    where chrId is an enum and binStart and binEnd are bp locations.
    Use the chromosome info from the header along with given binsize to build
    this table in numpy. Returns the table (with dimensions #bins x 3) and chr
    offset information with regard to bins
    """
    bins_array = []
    offsets = [0]
    by_chr_offsets = {}
    chr_n_bins = {}
    for c_idx in chrs:
        chr_name = chrs[c_idx][1]
        if chr_name == 'all':
            continue
        chr_size = chrs[c_idx][2]
        bin_start = 0
        chr_bin_count = 0
        while bin_start < chr_size:
            bin_end = bin_start + binsize
            if bin_end > chr_size:
                bin_end = chr_size
            bins_array.append([chr_name, bin_start, bin_end])
            chr_bin_count += 1
            bin_start = bin_end
        # offsets are the number of bins of all chromosomes prior to this one
        # by_chr_offsets are bin offsets for each chromosome, indexed by chr_idx
        by_chr_offsets[chrs[c_idx][0]] = offsets[-1]
        chr_n_bins[c_idx] = chr_bin_count
        if len(offsets) < len(chrs):
            offsets.append(int(math.ceil(chr_size/binsize) + offsets[-1]))
    return np.array(bins_array), np.array(offsets), chr_n_bins, by_chr_offsets


def write_bins(grp, req, chroms, bins, chr_n_bins, chr_footer_info, h5opts):
    """
    Write the bins table, which has columns: chrom, start (in bp), end (in bp),
    and one column for each normalization type, named for the norm type.
    'weight' is a reserved column for `cooler balance`.
    Chrom is an index which corresponds to the chroms table.
    Also write a dataset for each normalization vector within the hic file.
    """
    n_chroms = len(chroms)
    n_bins = len(bins)
    idmap = dict(zip(chroms, range(n_chroms)))
    # Convert chrom names to enum
    chrom_ids = [idmap[chrom.encode('UTF-8')] for chrom in bins[:, 0]]
    starts = [int(val) for val in bins[:, 1]]
    ends = [int(val) for val in bins[:, 2]]
    enum_dtype = h5py.special_dtype(enum=(CHROMID_DTYPE, idmap))
    # Store bins
    grp.create_dataset('chrom',
                       shape=(n_bins,),
                       dtype=enum_dtype,
                       data=chrom_ids,
                       **h5opts)
    grp.create_dataset('start',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=starts,
                       **h5opts)
    grp.create_dataset('end',
                       shape=(n_bins,),
                       dtype=COORD_DTYPE,
                       data=ends,
                       **h5opts)
    # write columns for normalization vectors
    for norm in NORMS:
        norm_data = []
        for chr_idx in chr_n_bins:
            if chr_idx not in chr_footer_info[norm]:
                error_str = (
                    'ERROR. Normalization vector %s does not exist for chr idx %s.'
                    '\nThis chromosome may not exist in the hic file. Try running'
                    ' with the "-e" option to exclude MT.' % (norm, chr_idx))
                force_exit(error_str, req)
            norm_vector = read_normalization_vector(req, chr_footer_info[norm][chr_idx])
            chr_bin_end = chr_n_bins[chr_idx]
            # possible issue to look into
            # norm_vector returned by read_normalization_vector has an extra
            # entry at the end (0.0) that is never used and causes the
            # norm vectors to be longer than n_bins for a given chr
            # restrict the length of norm_vector to chr_bin_end for now.
            norm_data.extend(list(map(norm_convert, norm_vector[:chr_bin_end])))
        if len(norm_data) != n_bins:
            error_str = (
                'ERROR. Length of normalization vector %s does not match the'
                ' number of bins.\nThis is likely a problem with the hic file' % (norm))
            force_exit(error_str, req)
        grp.create_dataset(
            norm,
            shape=(len(norm_data),),
            dtype=NORM_DTYPE,
            data=np.array(norm_data, dtype=NORM_DTYPE),
            **h5opts)


def norm_convert(val):
    """
    Convert between hic and cool normalization values
    """
    if val != 0.0:
        return 1 / val
    else:
        return np.inf


def prepare_chunks_dataset(grp, n_bins, h5opts):
    """
    Temporary dataset that is deleted before cooler file is finalized.
    Stores the chunks that are iteratively generated from the hic file.
    See generate_counts_chunk function for chunk details.
    """
    # max_size = n_bins * (n_bins - 1) // 2 + n_bins
    max_size = None
    init_size = 0
    grp.create_dataset(
        'chunks',
        shape=(init_size,),
        maxshape=(max_size,),
        chunks=(CHUNK_SIZE,),
        dtype=CHUNK_DTYPE,
        **h5opts)


def write_chrom_offset(grp, chr_offsets, h5opts):
    """
    Write the chrom offset information (bin offset information written later)
    """
    grp.create_dataset(
        'chrom_offset',
        shape=(len(chr_offsets),),
        dtype=CHROMOFFSET_DTYPE,
        data=chr_offsets,
        **h5opts)


def initialize_pixels(grp, n_bins, h5opts):
    """
    Initialize the pixels datasets with a given max_size (expression from
    cooler). These are resizable datasets that are chunked with size of
    CHUNK_SIZE and initialized with lenth of init_size.
    """
    # max_size = n_bins * (n_bins - 1) // 2 + n_bins
    max_size = None
    init_size = 0
    grp.create_dataset('bin1_id',
                       dtype=BIN_DTYPE,
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    grp.create_dataset('bin2_id',
                       dtype=BIN_DTYPE,
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)
    grp.create_dataset('count',
                       dtype=COUNT_DTYPE,
                       shape=(init_size,),
                       maxshape=(max_size,),
                       **h5opts)


def write_bin1_offset(grp, bin1_offsets, h5opts):
    """
    Write the bin1_offset information
    """
    grp.create_dataset(
        'bin1_offset',
        shape=(len(bin1_offsets),),
        dtype=BIN1OFFSET_DTYPE,
        data=bin1_offsets,
        **h5opts)


def write_pixels_chunk(outfile, resolution, chunks, low_mem):
    """
    Converts the chunks dumped from the hic file into a 2D pixels array,
    with columns bin1 (int idx), bin2 (int idx), and count (int).
    An offsets array is also generated, which gives the number of bins preceding
    any given bin (with regard to bin1 position). Thus, bin1[0] will always have
    an offset of 0 and bin1[1] will have an offset equal to the number of times
    bin1=0.
    The pixel array is sorted by bin1 first, then bin2. The following could be
    a 2D pixel array (where the columns refer to bin1, bin2, and count
    respectively)
    0   0   3
    0   1   1
    0   2   0
    0   3   1
    1   2   1
    1   3   2
    2   3   1

    After pixels group has been written and bin1_offsets written within
    the indexes group, the temporary chunks dataset is deleted.
    """
    h5opts = dict(compression='gzip', compression_opts=6)
    with h5py.File(outfile, "r+") as h5file:
        h5resolution = h5file[str(resolution)] if multi_res else h5file
        n_bins = h5resolution.attrs['nbins']
        # do I have to load this into memory?
        if low_mem:
            chunks = h5resolution['chunks'].value
        chunks.sort(order=['bin1_id', 'bin2_id'])
        nnz = chunks.shape[0]
        # write ordered chunk to pixels
        grp = h5resolution['pixels']
        for set_name in ['bin1_id', 'bin2_id', 'count']:
            dataset = grp[set_name]
            prev_size = dataset.shape[0]
            dataset.resize(prev_size + nnz, axis=0)
            dataset[prev_size:] = chunks[set_name]
        # reset temporary chunks dataset
        if 'chunks' in h5resolution.keys():
            del h5resolution['chunks']
            prepare_chunks_dataset(h5resolution, n_bins, h5opts)


def finalize_resolution_cool(outfile, resolution):
    """
    Finally, remove the chunks dataset and write bin1_offsets and nnz attr
    for one resolution in the cool file
    """
    h5opts = dict(compression='gzip', compression_opts=6)
    with h5py.File(outfile, "r+") as h5file:
        h5resolution = h5file[str(resolution)] if multi_res else h5file
        if 'chunks' in h5resolution.keys():
            del h5resolution['chunks']
        nnz = h5resolution['pixels']['count'].shape[0]
        n_bins = h5resolution.attrs['nbins']
        # now we need bin offsets
        # adapated from _writer.py in cooler
        offsets = np.zeros(n_bins + 1, dtype=BIN1OFFSET_DTYPE)
        curr_val = 0
        for start, length, value in zip(*rlencode(h5resolution['pixels']['bin1_id'], 1000000)):
            offsets[curr_val:value + 1] = start
            curr_val = value + 1
        offsets[curr_val:] = nnz
        write_bin1_offset(h5resolution['indexes'], offsets, h5opts)
        # update attributes with nnz
        info = {'nnz': nnz}
        h5resolution.attrs.update(info)


def rlencode(array, chunksize=None):
    """
    Run length encoding.
    Based on http://stackoverflow.com/a/32681075, which is based on the rle
    function from R.

    TAKEN FROM COOLER

    Parameters
    ----------
    x : 1D array_like
        Input array to encode
    dropna: bool, optional
        Drop all runs of NaNs.

    Returns
    -------
    start positions, run lengths, run values

    """
    where = np.flatnonzero
    if not isinstance(array, h5py.Dataset):
        array = np.asarray(array)
    n = len(array)
    if n == 0:
        return (np.array([], dtype=int),
                np.array([], dtype=int),
                np.array([], dtype=array.dtype))

    if chunksize is None:
        chunksize = n

    starts, values = [], []
    last_val = np.nan
    for i in range(0, n, chunksize):
        x = array[i:i+chunksize]
        locs = where(x[1:] != x[:-1]) + 1
        if x[0] != last_val:
            locs = np.r_[0, locs]
        starts.append(i + locs)
        values.append(x[locs])
        last_val = x[-1]
    starts = np.concatenate(starts)
    lengths = np.diff(np.r_[starts, n])
    values = np.concatenate(values)

    return starts, lengths, values


def write_zooms_for_higlass(h5res):
    """
    NOT CURRENTLY USED

    Add max-zoom column needed for higlass, but only if the resolutions are
    successively divisible by two. Otherwise, warn that this file will not be
    higlass compatible
    """
    resolutions = sorted(f['resolutions'].keys(), key=int, reverse=True)

    # test if resolutions are higlass compatible
    higlass_compat = True
    for i in range(len(resolutions)):
        if i == len(resolutions) - 1:
            break
        if resolutions[i] * 2 != resolutions[i+1]:
            higlass_compat = False
            break
    if not higlass_compat:
        print('WARNING: This hic file is not higlass compatible! Will not add [max-zoom] attribute.')
        return

    print('INFO: This hic file is higlass compatible! Adding [max-zoom] attribute.')
    max_zoom = len(resolutions) - 1

    # Assign max-zoom attribute
    h5res.attrs['max-zoom'] = max_zoom
    print('max-zoom: {}'.format(max_zoom))

    # Make links to zoom levels
    for i, res in enumerate(resolutions):
        print('zoom {}: {}'.format(i, res))
        h5res[str(i)] = h5py.SoftLink('/resolutions/{}'.format(res))


def hic2cool_convert(infile, outfile, resolution=0, exclude_MT=False, low_mem=False, command_line=False):
    """
    Main function that coordinates the reading of header and footer from infile
    and uses that information to parse the hic matrix.
    Opens outfile and writes in form of .cool file
    Params:
    <infile> str .hic filename
    <outfile> str .cool output filename
    <resolution> int bp bin size. If 0, use all. Defaults to 0.
                Final .cool structure will change depending on this param (see README)
    <exclude_MT> bool. If True, ignore MT contacts. Defaults to False.
    <low_mem> bool. If True, attempt to reduce memory used by writing pixel chunks
            to an intermediate dataset in the output cool file.
    <command_line> bool. True if executing from run_hic.py. Prompts hic headers
                be printed to stdout.
    """
    unit = 'BP'  # only using base pair unit for now
    resolution = int(resolution)
    req, used_chrs, resolutions, masteridx, genome = read_header(infile)
    if command_line:  # print hic header info for command line usage
        chr_names = [used_chrs[key][1] for key in used_chrs.keys()]
        print('################')
        print('### hic2cool ###')
        print('################')
        print('hic file header info:')
        print('Chromosomes: ', chr_names)
        print('Resolutions: ', resolutions)
        print('Genome: ', genome)
    if exclude_MT:  # remove chr25, which is MT, if this flag is set
        used_chrs.pop(25, None)
    # ensure user input binsize is a resolution supported by the hic file
    if resolution != 0 and resolution not in resolutions:
        error_str = (
            'ERROR. Given binsize (in bp) is not a supported resolution in '
            'this file.\nPlease use 0 (all resolutions) or use one of: ' +
            str(resolutions))
        force_exit(error_str, req)
    use_resolutions = resolutions if resolution == 0 else [resolution]
    global multi_res
    multi_res = resolution == 0
    for binsize in use_resolutions:
        t_start = time.time()
        req, pair_footer_info, chr_footer_info = read_footer(
            req, masteridx, unit, binsize)
        # initialize cooler file. return per resolution bin offset maps
        chr_offset_map = initialize_cool(outfile, req, used_chrs, genome,
                                        binsize, chr_footer_info, low_mem)
        covered_chr_pairs = []
        for chr_a in used_chrs:
            total_chunk = np.zeros(shape=0, dtype=CHUNK_DTYPE) if not low_mem else None
            if used_chrs[chr_a][1].lower() == 'all':
                continue
            for chr_b in used_chrs:
                if used_chrs[chr_b][1].lower() == 'all':
                    continue
                c1 = min(chr_a, chr_b)
                c2 = max(chr_a, chr_b)
                # ensure this is true
                # since matrices are upper triangular, no need to cover c1-c2
                # and c2-c1 reciprocally
                if str(c1) + "_" + str(c2) in covered_chr_pairs:
                    continue
                tmp_chunk = parse_hic(req, outfile, used_chrs[c1], used_chrs[c2],
                          unit, binsize, covered_chr_pairs, pair_footer_info,
                          chr_offset_map, low_mem)
                if not low_mem:
                    total_chunk = build_chunk(total_chunk, tmp_chunk)
                    del tmp_chunk
            # write at the end of every chr_a
            write_pixels_chunk(outfile, binsize, total_chunk, low_mem)
        # finalize to remove chunks and write a bit of metadata
        finalize_resolution_cool(outfile, binsize)
        t_parse = time.time()
        elapsed_parse = t_parse - t_start
        print('Resolution %s took: %s seconds.' % (binsize, elapsed_parse))
    req.close()


def force_exit(message, req):
    """
    Exit the program due to some error. Print out message and close the given
    input files.
    """
    req.close()
    print(message, file=sys.stderr)
    sys.exit()
