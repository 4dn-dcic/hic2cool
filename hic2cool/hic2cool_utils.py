#### hic2cool ####

# CLI for a converter between .hic files (from juicer) and .cool files (for
# cooler).  Both hic and cool files describe Hi-C contact matrices.
# Intended to be lightweight, this is a stand-alone Python file written by
# Carl Vitzthum of the HMS DBMI Park lab. Originally published 1/26/17.

# This code was originally based off of the straw project by Neva C. Durand and
# Yue Wu (https://github.com/theaidenlab/straw).
# The cooler file writing was based off of much of the CLI code contained in
# this repo: https://github.com/mirnylab/cooler.

# The most import method is hic2cool_convert, which can be imported if the package is pip installed.

# Usage on the command line is: python hic2cool.py infile outfile resolution normalization include_MT
# If an invalid resolution is given (i.e. one that does not correspond to a hic resolution,
# the program will terminate and prompt you to enter a valid one)
# See the README or main() function in this file for more usage information.

from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import struct
import zlib
import numpy as np
import h5py
import math
from collections import OrderedDict
from ._version import __version__


#read function
def readcstr(f):
    # buf = bytearray()
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            # return str(buf,encoding="utf-8", errors="strict")
            return buf.decode("utf-8", errors="ignore")
        else:
            buf += b
            # buf.append(b)


def read_header(infile):
    """
    Takes in a .hic file and returns a dictionary containing information about
    the chromosome. Keys are chromosome index numbers (0 through # of chroms contained
    in file) and values are [chr idx (int), chr name (str), chrom length (str)].
    Returns the masterindex used by the file as well as the open file object.
    """
    req=open(infile, 'rb')
    chrs = {}
    resolutions = []
    magic_string = struct.unpack('<3s', req.read(3))[0]
    req.read(1)
    if (magic_string != b"HIC"):
        print('This does not appear to be a HiC file; magic string is incorrect')
        sys.exit()
    global version
    version = struct.unpack('<i',req.read(4))[0]
    masterindex = struct.unpack('<q',req.read(8))[0]
    genome = b""
    c=req.read(1)
    while (c != b'\0'):
        genome += c
        c=req.read(1)
    genome = genome.decode('ascii')
    nattributes = struct.unpack('<i',req.read(4))[0]
    for x in range(nattributes):
        key = readcstr(req)
        value = readcstr(req)
    nChrs = struct.unpack('<i',req.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(req)
        length = struct.unpack('<i',req.read(4))[0]
        if name and length:
            formatted_name = 'chr' + name if ('all' not in name.lower() and 'chr' not in name.lower()) else name
            formatted_name = 'chrM' if formatted_name == 'chrMT' else formatted_name
            chrs[i] = [i, formatted_name, length]
    nBpRes = struct.unpack('<i',req.read(4))[0]
    # find bp delimited resolutions supported by the hic file
    for x in range(0, nBpRes):
      res = struct.unpack('<i',req.read(4))[0]
      resolutions.append(res)
    return req, chrs, resolutions, masterindex, genome


def read_footer(req, master, norm, unit, resolution):
    """
    Takes in an open hic file and generates two dictionaries. pair_footer_info
    contains the file position of info for any given chromosome pair (formatted
    as a string). Chr_footer_info gives chromosome-level size and position info
    relative to the file. This way, this function only has to run once
    All of the unused read() code is used to find the correct place in the file,
    supposedly. This is code from straw.
    """
    pair_footer_info={}
    chr_footer_info={}
    req.seek(master)
    nBytes = struct.unpack('<i', req.read(4))[0]
    nEntries = struct.unpack('<i', req.read(4))[0]
    found = False
    for i in range(nEntries):
        stri = readcstr(req)
        fpos = struct.unpack('<q', req.read(8))[0]
        sizeinbytes = struct.unpack('<i', req.read(4))[0]
        pair_footer_info[stri] = fpos
    nExpectedValues = struct.unpack('<i',req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        binSize = struct.unpack('<i',req.read(4))[0]
        nValues = struct.unpack('<i',req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack('<d',req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i',req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack('<i',req.read(4))[0]
            v = struct.unpack('<d',req.read(8))[0]
    nExpectedValues = struct.unpack('<i',req.read(4))[0]
    for i in range(nExpectedValues):
        str_ = readcstr(req)
        str_ = readcstr(req)
        binSize = struct.unpack('<i',req.read(4))[0]
        nValues = struct.unpack('<i',req.read(4))[0]
        for j in range(nValues):
            v = struct.unpack('<d',req.read(8))[0]
        nNormalizationFactors = struct.unpack('<i',req.read(4))[0]
        for j in range(nNormalizationFactors):
            chrIdx = struct.unpack('<i',req.read(4))[0]
            v = struct.unpack('<d',req.read(8))[0]
    nEntries = struct.unpack('<i',req.read(4))[0]
    for i in range(nEntries):
        normtype = readcstr(req)
        chrIdx = struct.unpack('<i',req.read(4))[0]
        unit1 = readcstr(req)
        resolution1 = struct.unpack('<i',req.read(4))[0]
        filePosition = struct.unpack('<q',req.read(8))[0]
        sizeInBytes = struct.unpack('<i',req.read(4))[0]
        if (normtype==norm and unit1==unit and resolution1==resolution):
            chr_footer_info[chrIdx] = {'position':filePosition, 'size':sizeInBytes}

    return req, pair_footer_info, chr_footer_info


#FUN(fin, unit, resolution, blockMap) Return [storeBlockData, myBlockBinCount, myBlockColumnCount]
def read_matrix_zoom_data(req, myunit, mybinsize, blockMap):
    unit = readcstr(req)
    temp = struct.unpack('<i',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    temp2 = struct.unpack('<f',req.read(4))[0]
    binSize = struct.unpack('<i',req.read(4))[0]
    blockBinCount = struct.unpack('<i',req.read(4))[0]
    blockColumnCount = struct.unpack('<i',req.read(4))[0]
    storeBlockData = False
    #for the initial
    myBlockBinCount = -1
    myBlockColumnCount = -1
    if (myunit==unit and mybinsize==binSize):
        myBlockBinCount=blockBinCount
        myBlockColumnCount=blockColumnCount
        storeBlockData=True
    nBlocks = struct.unpack('<i',req.read(4))[0]
    for b in range(nBlocks):
        blockNumber = struct.unpack('<i',req.read(4))[0]
        filePosition = struct.unpack('<q',req.read(8))[0]
        blockSizeInBytes = struct.unpack('<i',req.read(4))[0]
        entry=dict()
        entry['size'] = blockSizeInBytes
        entry['position'] = filePosition
        if (storeBlockData):
            blockMap[blockNumber] = entry
    return [storeBlockData, myBlockBinCount, myBlockColumnCount]


#FUN(fin, myFilePos, unit, binsize, blockMap) Return [blockBinCount, blockColumnCount]
def read_matrix(req, myFilePos, unit, binsize, blockMap):
    req.seek(myFilePos)
    c1 = struct.unpack('<i',req.read(4))[0]
    c2 = struct.unpack('<i',req.read(4))[0]
    nRes = struct.unpack('<i',req.read(4))[0]
    i = 0
    found = False
    blockBinCount = -1
    blockColumnCount = -1
    while (i<nRes and (not found)):
        list1 = read_matrix_zoom_data(req, unit, binsize, blockMap)
        found = list1[0]
        if(list1[1]!=-1 and list1[2]!=-1):
            blockBinCount = list1[1]
            blockColumnCount = list1[2]
        i=i+1
    if (not found):
        print("Error finding block data\n")
        return -1
    return [blockBinCount, blockColumnCount]


#FUN(regionIndices, blockBinCount, blockColumnCount, intra) Return blocksSet
def get_block_numbers_for_region_from_bin_position(regionIndices, blockBinCount, blockColumnCount, intra):
    col1=int(regionIndices[0]/blockBinCount)
    col2=int((regionIndices[1]+1)/blockBinCount)
    row1=int(regionIndices[2]/blockBinCount)
    row2=int((regionIndices[3]+1)/blockBinCount)
    blocksSet=set()
    for r in range(row1, row2+1):
        for c in range(col1, col2+1):
            blockNumber=r*blockColumnCount+c
            blocksSet.add(blockNumber)
    if (intra):
        for r in range(col1, col2+1):
            for c in range(row1, row2+1):
                blockNumber=r*blockColumnCount+c
                blocksSet.add(blockNumber)
    return blocksSet


#FUN(fin, blockNumber, blockMap)
def read_block(req, blockNumber, blockMap):
    idx=dict()
    if(blockNumber in blockMap):
        idx=blockMap[blockNumber]
    else:
        idx['size']=0
        idx['position']=0
    if (idx['size']==0):
        return []
    req.seek(idx['position'])
    compressedBytes = req.read(idx['size'])
    uncompressedBytes = zlib.decompress(compressedBytes)
    nRecords = struct.unpack('<i',uncompressedBytes[0:4])[0]
    v=[]
    global version
    if (version < 7):
        for i in range(nRecords):
            binX = struct.unpack('<i',uncompressedBytes[(12*i):(12*i+4)])[0]
            binY = struct.unpack('<i',uncompressedBytes[(12*i+4):(12*i+8)])[0]
            counts = struct.unpack('<f',uncompressedBytes[(12*i+8):(12*i+12)])[0]
            record = dict()
            record['binX'] = binX
            record['binY'] = binY
            record['counts'] = counts
            v.append(record)
    else:
        binXOffset = struct.unpack('<i',uncompressedBytes[4:8])[0]
        binYOffset = struct.unpack('<i',uncompressedBytes[8:12])[0]
        useShort = struct.unpack('<b',uncompressedBytes[12:13])[0]
        type_ = struct.unpack('<b',uncompressedBytes[13:14])[0]
        index=0
        if (type_==1):
            rowCount = struct.unpack('<h',uncompressedBytes[14:16])[0]
            temp = 16
            for i in range(rowCount):
                y = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                temp=temp+2
                binY = y + binYOffset
                colCount = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                temp=temp+2
                for j in range(colCount):
                    x = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                    temp=temp+2
                    binX = binXOffset + x
                    if (useShort==0):
                        c = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                        temp=temp+2
                        counts = c
                    else:
                        counts = struct.unpack('<f',uncompressedBytes[temp:(temp+4)])[0]
                        temp=temp+4
                    record = dict()
                    record['binX'] = binX
                    record['binY'] = binY
                    record['counts'] = counts
                    v.append(record)
                    index = index + 1
        elif (type_== 2):
            temp=14
            nPts = struct.unpack('<i',uncompressedBytes[temp:(temp+4)])[0]
            temp=temp+4
            w = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
            temp=temp+2
            for i in range(nPts):
                row=int(i/w)
                col=i-row*w
                bin1=int(binXOffset+col)
                bin2=int(binYOffset+row)
                if (useShort==0):
                    c = struct.unpack('<h',uncompressedBytes[temp:(temp+2)])[0]
                    temp=temp+2
                    if (c != -32768):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = c
                        v.append(record)
                        index = index + 1
                else:
                    counts = struct.unpack('<f',uncompressedBytes[temp:(temp+4)])[0]
                    temp=temp+4
                    if (countsnot != 0x7fc00000):
                        record = dict()
                        record['binX'] = bin1
                        record['binY'] = bin2
                        record['counts'] = counts
                        v.append(record)
                        index = index + 1
    return v


#FUN(fin, entry) Return Norm
def read_normalization_vector(req, entry):
    req.seek(entry['position'])
    nValues = struct.unpack('<i',req.read(4))[0]
    value = []
    for i in range(nValues):
        d = struct.unpack('<d',req.read(8))[0]
        value.append(d)
    return value


def parse_hic(norm, req, h5file, chr1, chr2, unit, binsize, covered_chr_pairs, pair_footer_info, chr_footer_info, bin_map, count_map):
    """
    Adapted from the straw() function in the original straw package.
    Mainly, since all chroms are iterated over, the read_header and read_footer
    functions were placed outside of straw() and made to be reusable across
    any chromosome pair.
    Main function is to build a bin_map, which contains normalization values
    for every bin, and a count_map, which is a nested dictionary which contains
    the contact count for any two bins.
    """
    blockMap = {}
    magic_string = ""
    if (not (norm=="NONE" or norm=="VC" or norm=="VC_SQRT" or norm=="KR")):
        print("Norm specified incorrectly, must be one of <NONE/VC/VC_SQRT/KR>")
        force_exit(warn_string, req, h5file)
    if (not (unit=="BP" or unit=="FRAG")):
        print("Unit specified incorrectly, must be one of <BP/FRAG>")
        force_exit(warn_string, req, h5file)
    chr1ind=chr1[0]
    chr2ind=chr2[0]
    c1pos1=0
    c1pos2=int(chr1[2])
    c2pos1=0
    c2pos2=int(chr2[2])
    c1=min(chr1ind,chr2ind)
    c2=max(chr1ind,chr2ind)
    chr_key = str(c1) + "_" + str(c2)
    origRegionIndices=[]
    regionIndices=[]
    if (chr1ind > chr2ind):
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        regionIndices.append(int(c2pos1/binsize))
        regionIndices.append(int(c2pos2/binsize))
        regionIndices.append(int(c1pos1/binsize))
        regionIndices.append(int(c1pos2/binsize))
    else:
        origRegionIndices.append(c1pos1)
        origRegionIndices.append(c1pos2)
        origRegionIndices.append(c2pos1)
        origRegionIndices.append(c2pos2)
        regionIndices.append(int(c1pos1/binsize))
        regionIndices.append(int(c1pos2/binsize))
        regionIndices.append(int(c2pos1/binsize))
        regionIndices.append(int(c2pos2/binsize))
    try:
        pair_footer_info[chr_key]
    except KeyError:
        warn_string = ('ERROR. There is a discrepancy between the chrs declared in the infile header and the actual information it contains.\nThe intersection between ' + chr1[1] + ' and ' + chr2[1] + ' could not be found in the file.')
        force_exit(warn_string, req, h5file)
    myFilePos=pair_footer_info[chr_key]
    if (norm != "NONE"):
        c1Norm = read_normalization_vector(req, chr_footer_info[c1])
        c2Norm = read_normalization_vector(req, chr_footer_info[c2])
    list1 = read_matrix(req, myFilePos, unit, binsize, blockMap)
    blockBinCount=list1[0]
    blockColumnCount=list1[1]
    blockNumbers = get_block_numbers_for_region_from_bin_position(regionIndices, blockBinCount, blockColumnCount, c1==c2)
    for i_set in (blockNumbers):
        records=read_block(req, i_set, blockMap)
        for j in range(len(records)):
            rec=records[j]
            # x and y are bin indices
            x=rec['binX']
            y=rec['binY']
            c=rec['counts']
            if (norm != "NONE"):
                normBinX = c1Norm[x]
                normBinY = c2Norm[y]
                if normBinX != 0.0:
                    nX = 1/normBinX
                else:
                    nX = 'inf'
                if normBinY != 0.0:
                    nY = 1/normBinY
                else:
                    nY = 'inf'
            else:
                nX = 1.0
                nY = 1.0
            if ((x>=origRegionIndices[0] and x<=origRegionIndices[1] and y>=origRegionIndices[2] and y<=origRegionIndices[3]) or ((c1==c2) and y>=origRegionIndices[0] and y<=origRegionIndices[1] and x>= origRegionIndices[2] and x<=origRegionIndices[3])):
                c1_key = str(c1) + ":" + str(x)
                c2_key = str(c2) + ":" + str(y)
                bin_map[c1_key] = nX
                bin_map[c2_key] = nY
                if c1_key not in count_map:
                    count_map[c1_key] = {}
                if c2_key in count_map[c1_key]:
                    print('\nWARNING: multiple count entries found for the following pixel\n', c1_key, ' and ', c2_key, '   (in form chrom_idx:bin)\n')
                    print('c1:  ', str(c1)+':'+str(x*binsize)+'|'+str(c2)+':'+str(y*binsize),'. Conflicting count found: ', c, '. Will use: ', count_map[c1_key][c2_key])
                else:
                    count_map[c1_key][c2_key] = c
    covered_chr_pairs.append(chr_key)


def write_cool(h5res, chr_info, binsize, bin_map, count_map, norm, genome):
    """
    Use various information to write a cooler file in HDF5 format. Tables
    included are chroms, bins, pixels, and indexes.
    For an in-depth explanation of the data structure, please see:
    http://cooler.readthedocs.io/en/latest/datamodel.html
    """
    grp = h5res.create_group('chroms')
    h5opts = dict(compression='gzip', compression_opts=6)
    chr_names = write_chroms(grp, chr_info, h5opts)
    # create the 'bins' table required by cooler_root
    bin_table, chr_offsets, by_chr_offset_map = create_bins(chr_info, binsize, bin_map)
    grp = h5res.create_group('bins')
    write_bins(grp, chr_names, bin_table, h5opts, norm)
    pixels_table, bin1_offsets = create_pixels(count_map, by_chr_offset_map, len(bin_table))
    grp = h5res.create_group('pixels')
    write_pixels(grp, pixels_table, h5opts)
    grp = h5res.create_group('indexes')
    write_indexes(grp, chr_offsets, bin1_offsets, h5opts)
    # update info attributes of the h5res
    info = {}
    info['nchroms'] = len(chr_names)
    info['nbins'] = len(bin_table)
    info['nnz'] = len(pixels_table)
    info['bin-type'] = 'fixed'
    info['bin-size'] = binsize
    info['format'] = 'HDF5::Cooler'
    info['generated-by'] = 'hic2cool-' + __version__
    info['genome-assembly'] = genome
    h5res.attrs.update(info)



def write_chroms(grp, chrs, h5opts):
    """
    Write the chroms table, which includes chromosome names and length
    """
    # format chr_names and chr_lengths np arrays
    chr_names = np.array([val[1] for val in chrs.values() if val[1] != 'All'], dtype=np.dtype('S32'))
    chr_lengths = np.array([val[2] for val in chrs.values() if val[1] != 'All'])
    grp.create_dataset('name', shape=(len(chr_names),), dtype=np.dtype('S32'), data=chr_names, **h5opts)
    grp.create_dataset('length', shape=(len(chr_names),), dtype=np.int32, data=chr_lengths, **h5opts)
    return chr_names


def create_bins(chrs, binsize, bin_map):
    """
    Cooler requires a bin file with each line having <chrId, binStart, binEnd, weight>
    where binStart and binEnd are bp locations and weight is a float.
    Use the chromosome info from the header along with given binsize to build
    this table in numpy
    Returns the table (with dimensions #bins x 4) and chr offset information
    with regard to bins
    """
    bins_array = []
    offsets = [0]
    by_chr_offsets = {}
    for c_idx in chrs:
        chr_name = chrs[c_idx][1]
        if chr_name == 'All':
            continue
        chr_size = chrs[c_idx][2]
        bin_start = 0
        while bin_start < chr_size:
            bin_end = bin_start + binsize
            if bin_end > chr_size:
                bin_end = chr_size
            weight_key = str(chrs[c_idx][0]) + ':' + str(int(math.ceil(bin_start/binsize)))
            weight = bin_map[weight_key] if weight_key in bin_map else 1.0
            bins_array.append([chr_name, bin_start, bin_end, weight])
            bin_start = bin_end
        # offsets are the number of bins of all chromosomes prior to this one
        # by_chr_offsets are bin offsets for each chromosome, indexed by chr_idx
        by_chr_offsets[chrs[c_idx][0]] = offsets[-1]
        if len(offsets) < len(chrs):
            offsets.append(math.ceil(chr_size/binsize) + offsets[-1])
    return np.array(bins_array), np.array(offsets), by_chr_offsets


def write_bins(grp, chroms, bins, h5opts, norm):
    """
    Write the bins table, which has columns: chrom, start (in bp), end (in bp),
    and weight (float). Chrom is an index which corresponds to the chroms table.
    Only writes weight column if norm != "NONE"
    """
    n_chroms = len(chroms)
    n_bins = len(bins)
    idmap = dict(zip(chroms, range(n_chroms)))
    # Convert chrom names to enum
    chrom_ids = [idmap[chrom.encode('UTF-8')] for chrom in bins[:,0]]
    starts = [int(val) for val in bins[:,1]]
    ends = [int(val) for val in bins[:,2]]
    weights = [float(val) for val in bins[:,3]]
    enum_dtype = h5py.special_dtype(enum=(np.int32, idmap))
    grp.create_dataset('chrom', shape=(n_bins,), dtype=enum_dtype, data=chrom_ids, **h5opts)
    grp.create_dataset('start',  shape=(n_bins,), dtype=np.int32, data=starts, **h5opts)
    grp.create_dataset('end', shape=(n_bins,), dtype=np.int32, data=ends, **h5opts)
    if norm != "NONE":
        grp.create_dataset('weight', shape=(n_bins,), dtype=np.float64, data=weights, **h5opts)


def create_pixels(count_map, by_chr_offset_map, n_bins):
    """
    Converts the countmap dumped by from the hic file into a 2D pixels array,
    with columns bin1 (int idx), bin2 (int idx), and count (int).
    An offset array is also generated, which gives the number of bins preceding
    any given bin (with regard to bin1 position). Thus, bin1[0] will always have
    an offset of 0 and bin1[1] will have an offset equal to the number of times
    bin1=0.
    The pixel array is sorted by bin1 first, then bin2. The following could be
    a 2D pixel array (where the columns refer to bin1, bin2, and count respectively)
    0   0   3
    0   1   1
    0   2   0
    0   3   1
    1   2   1
    1   3   2
    2   3   1
    """
    pixel_2d_array = []
    bin1_offsets = []
    # order count_map by bin2 then bin1
    for key in count_map:
        # turn inner bin2 keys from chr_key to bin number
        converted_bin2 = {chr_key_to_bin(k, by_chr_offset_map): v for k,v in count_map[key].items()}
        # order them by bin number
        count_map[key] = OrderedDict(sorted(converted_bin2.items(), key=lambda t: t[0]))
    # do the same for bin1 (chr_key to bin number)
    converted_bin1 = {chr_key_to_bin(k, by_chr_offset_map): v for k,v in count_map.items()}
    ordered_counts_by_bins = OrderedDict(sorted(converted_bin1.items(), key=lambda t: t[0]))
    # build the 2d array. Columns are [bin1, bin2, count]. Sorted by bin1 then bin2 (ascending order)
    # loop through all possible bin indices, adding them to the pixel table if they are a valid bin1_idx
    # every index is added
    for bin1_idx in range(0, n_bins):
        # append the row offset for this bin1 in the pixels table
        bin1_offsets.append(len(pixel_2d_array))
        # loop through bin2s and counts
        if bin1_idx in ordered_counts_by_bins:
            for bin2, count in ordered_counts_by_bins[bin1_idx].items():
                pixel_2d_array.append([bin1_idx, bin2, count])
    return np.array(pixel_2d_array), np.array(bin1_offsets)


def chr_key_to_bin(key, offset_map):
    """
    Function made to covert chr_keys (in form <chridx>:<bin_idx>) to an offset
    bin idx using the offset_map (which uses chr_idxs as keys)
    """
    split_key = key.strip().split(':')
    chr_idx = int(split_key[0])
    offset = offset_map[chr_idx]
    return offset + int(split_key[1])


def write_pixels(grp, pixels, h5opts):
    """
    Write the pixel table. Columns are bin1ID, bin2ID, count
    """
    n_pixels = len(pixels)
    grp.create_dataset('bin1_id', shape=(n_pixels,), dtype=np.int64, data=pixels[:,0], **h5opts)
    grp.create_dataset('bin2_id',  shape=(n_pixels,), dtype=np.int64, data=pixels[:,1], **h5opts)
    grp.create_dataset('count', shape=(n_pixels,), dtype=np.int32, data=pixels[:,2], **h5opts)


def write_indexes(grp, chr_offsets, bin1_offsets, h5opts):
    """
    Write the chrom offset information and bin offset information
    """
    grp.create_dataset('chrom_offset', shape=(len(chr_offsets),), dtype=np.int32, data=chr_offsets, **h5opts)
    grp.create_dataset('bin1_offset',  shape=(len(bin1_offsets),), dtype=np.int32, data=bin1_offsets, **h5opts)


def hic2cool_convert(infile, outfile, resolution=0, norm='KR', exclude_MT=False, command_line=False):
    """
    Main function that coordinates the reading of header and footer from infile
    and uses that information to parse the hic matrix.
    Opens outfile and writes in form of .cool file
    Params:
    <infile> str .hic filename
    <outfile> str .cool output filename
    <resolution> int bp bin size. If 0, use all. Defaults to 0.
                Final .cool structure will change depending on this param (see README)
    <norm> str normalization type. Defaults to KR, optionally NONE, VC, or VC_SQRT
    <exclude_MT> bool. If True, ignore MT contacts. Defaults to False.
    <command_line> bool. True if executing from run_hic.py. Prompts hic headers
                be printed to stdout.
    """
    unit='BP' # only using base pair unit for now
    resolution = int(resolution)
    req, used_chrs, resolutions, masteridx, genome = read_header(infile)
    if command_line: # print hic header info for command line usage
        chr_names = [used_chrs[key][1] for key in used_chrs.keys()]
        print('################')
        print('### hic2cool ###')
        print('################')
        print('hic file header info:')
        print('Chromosomes: ', chr_names)
        print('Resolutions: ', resolutions)
        print('Genome: ', genome)
    if exclude_MT: # remove chr25, which is MT, if this flag is set
        used_chrs.pop(25, None)
    if resolution == 0: # use multi-res formatting
        h5file = h5py.File(outfile, 'w')
        h5resolutions = h5file.create_group('resolutions')
    else:
        h5file = h5py.File(outfile, 'w')
        h5res = h5file
    # ensure user input binsize is a resolution supported by the hic file
    if resolution != 0 and resolution not in resolutions:
        error_str = ('ERROR. Given binsize (in bp) is not a supported resolution in this file.\nPlease use 0 (all resolutions) or use one of: ' + resolutions)
        force_exit(error_str, req, h5file)
    use_resolutions = resolutions if resolution == 0 else [resolution]
    for binsize in use_resolutions:
        bin_map = {}
        count_map = {}
        req, pair_footer_info, chr_footer_info = read_footer(req, masteridx, norm, unit, binsize)
        covered_chr_pairs = []
        for chr_x in used_chrs:
            if used_chrs[chr_x][1].lower() == 'all':
                continue
            for chr_y in used_chrs:
                if used_chrs[chr_y][1].lower() == 'all':
                    continue
                c1=min(chr_x,chr_y)
                c2=max(chr_x,chr_y)
                # ensure this is true
                # since matrices are upper triangular, no need to cover c1-c2 and c2-c1 reciprocally
                if str(c1) + "_" + str(c2) in covered_chr_pairs:
                    continue
                parse_hic(norm, req, h5file, used_chrs[c1], used_chrs[c2], unit, binsize, covered_chr_pairs, pair_footer_info, chr_footer_info, bin_map, count_map)
        if resolution == 0:
            h5res = h5resolutions.create_group(str(binsize))
        write_cool(h5res, used_chrs, binsize, bin_map, count_map, norm, genome)
    req.close()
    h5file.close()


def force_exit(message, req, h5file):
    """
    Exit the program due to some error. Print out message and close the given
    input files.
    """
    req.close()
    h5file.close()
    print(message, file=sys.stderr)
    sys.exit()
