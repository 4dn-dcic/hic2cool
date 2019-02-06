from __future__ import absolute_import, division, print_function, unicode_literals
import sys
import struct
import zlib
import numpy as np
import h5py
import math
import cooler
import pandas as pd
from collections import OrderedDict

"""
This script uses the opposite approach of hic2cool. Whereas hic2cool starts
with a hic files and generates a cool or mcool file from scratch with the hic
normalization vectors included, this script requires an exisiting cool/mcool
file and adds the hic normalization vectors to it manually.
"""

version = 'dummy'


NORMS = ['VC', 'VC_SQRT', 'KR']


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
            formatted_name = ('chr' + name if ('all' not in name.lower() and
                              'chr' not in name.lower()) else name)
            formatted_name = ('chrM' if formatted_name.lower() == 'chrmt' else
                              formatted_name)
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
    possibleNorms = req.read(4)
    if not possibleNorms:
        return req, pair_footer_info, chr_footer_info
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


#FUN(fin, entry) Return Norm
def read_normalization_vector(req, entry):
    req.seek(entry['position'])
    nValues = struct.unpack('<i',req.read(4))[0]
    value = []
    for i in range(nValues):
        d = struct.unpack('<d',req.read(8))[0]
        value.append(d)
    return value


def parse_norm(norm, req, chr1, chr2, unit, binsize, covered_chr_pairs,
               pair_footer_info, chr_footer_info, chrom_map):
    """
    Adapted from the straw() function in the original straw package.
    Mainly, since all chroms are iterated over, the read_header and read_footer
    functions were placed outside of straw() and made to be reusable across
    any chromosome pair.
    Main function is to build a bin_map, which contains normalization values
    for every bin, and a count_map, which is a nested dictionary which contains
    the contact count for any two bins.

    """
    magic_string = ""
    if (not (norm=="VC" or norm=="VC_SQRT" or norm=="KR")):
        print("Norm specified incorrectly, must be one of <NONE/VC/VC_SQRT/KR>")
        force_exit(warn_string, req)
    if (not (unit=="BP" or unit=="FRAG")):
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

    try:
        pair_footer_info[chr_key]
    except KeyError:
        warn_string = (
            'ERROR. There is a discrepancy between the chrs declared in the ' +
            'infile header and the actual information it contains.\nThe '
            'intersection between ' + chr1[1] + ' and ' + chr2[1] +
            ' could not be found in the file.')
        force_exit(warn_string, req)
    myFilePos = pair_footer_info[chr_key]

    # chr_footer_info will be an empty dictionary if norm not in file
    if norm != "NONE" and chr_footer_info:
        c1Norm = read_normalization_vector(req, chr_footer_info[c1])
        c2Norm = read_normalization_vector(req, chr_footer_info[c2])
        chrom_map[chr1[1]] = c1Norm
        chrom_map[chr2[1]] = c2Norm

    covered_chr_pairs.append(chr_key)


def hic2cool_extractnorms(infile, outfile, resolution=0,
                          exclude_MT=False, command_line=False):
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
    from collections import OrderedDict
    import cooler

    unit = 'BP' # only using base pair unit for now
    resolution = int(resolution)

    req, used_chrs, resolutions, masteridx, genome = read_header(infile)

    chromosomes = [used_chrs[i][1] for i in range(1, len(used_chrs))]
    lengths = [used_chrs[i][2] for i in range(1, len(used_chrs))]
    chromsizes = pd.Series(index=chromosomes, data=lengths)

    chr_names = [used_chrs[key][1] for key in used_chrs.keys()]
    if command_line: # print hic header info for command line usage
        print('################')
        print('### hic2cool ###')
        print('################')
        print('hic file header info:')
        print('Chromosomes: ', chr_names)
        print('Resolutions: ', resolutions)
        print('Genome: ', genome)

    if exclude_MT: # remove mitchondrial chr by name if this flag is set
        # find index of chr with formatted name == chrM
        found_idxs = [idx for idx, fv in used_chrs.items() if fv[1].lower() == 'chrm']
        if len(found_idxs) == 1:
            used_chrs.pop(found_idxs[0], None)
        elif len(found_idxs) > 1:
            error_str = (
                'ERROR. More than one chromosome with name chrM was found when '
                ' attempting to exclude MT. Found chromosomes: %s' % chr_names
            )
            force_exit(error_str, req)
        else:
            if command_line: print('No chrM found for exlcude MT. Found chromosomes: %s' % chr_names)

    # ensure user input binsize is a resolution supported by the hic file
    if resolution != 0 and resolution not in resolutions:
        error_str = (
            'ERROR. Given binsize (in bp) is not a supported resolution in ' +
            'this file.\nPlease use 0 (all resolutions) or use one of: ' +
            resolutions)
        force_exit(error_str, req)

    use_resolutions = resolutions if resolution == 0 else [resolution]

    cooler_groups = {}
    for path in cooler.io.ls(outfile):
        binsize = cooler.Cooler(outfile + '::' + path).info['bin-size']
        cooler_groups[binsize] = path
    print('MCOOL contents:')
    print(cooler_groups)

    for norm in NORMS:
        print('Norm:', norm)

        for binsize in use_resolutions:
            chrom_map = {}
            bins = cooler.binnify(chromsizes, binsize)
            req, pair_footer_info, chr_footer_info = read_footer(
                req, masteridx, norm, unit, binsize)

            covered_chr_pairs = []
            for chr_x in used_chrs:
                if used_chrs[chr_x][1].lower() == 'all':
                    continue
                for chr_y in used_chrs:
                    if used_chrs[chr_y][1].lower() == 'all':
                        continue
                    c1 = min(chr_x, chr_y)
                    c2 = max(chr_x, chr_y)
                    # ensure this is true
                    # since matrices are upper triangular, no need to cover
                    # c1-c2 and c2-c1 reciprocally
                    if str(c1) + "_" + str(c2) in covered_chr_pairs:
                        continue
                    parse_norm(
                        norm,
                        req,
                        used_chrs[c1],
                        used_chrs[c2],
                        unit,
                        binsize,
                        covered_chr_pairs,
                        pair_footer_info,
                        chr_footer_info,
                        chrom_map
                    )

            # case where normalization info not present
            if not chrom_map:
                print('Normalization not present in hic file at resolution %s. Continuing...' % binsize)
                continue

            lengths_in_bins = bins.groupby('chrom').size()
            # hic normalization vector lengths have inconsistent lengths...
            # truncate appropriately
            vector = np.concatenate([
                chrom_map[chrom][:lengths_in_bins.loc[chrom]]
                for chrom in chromosomes
            ])
            bins[norm] = vector
            print('Resolution:', binsize)
            print(bins.head())
            print('Writing to cool file...')
            group_path = cooler_groups[binsize]
            cooler.io.append(
                outfile + '::' + group_path,
                'bins',
                {norm: bins[norm].values},
                force=True)
    req.close()


def force_exit(message, req):
    """
    Exit the program due to some error. Print out message and close the given
    input files.
    """
    req.close()
    print(message, file=sys.stderr)
    sys.exit()


if __name__ == '__main__':

    import argparse

    def main():
        """
        Execute the program from the command line
        Args are:
        python hic2hic2cool_extractnorms.py <infile (.hic)> <outfile (.cool)>
        <resolutions desired (defaults to all, optionally bp int)>
        <normalization type (defaults to 'KR', optionally 'NONE', 'VC', or 'VC_SQRT')>
        <exclude MT (default False)>

        """
        parser = argparse.ArgumentParser()
        parser.add_argument("infile", help=".hic input file")
        parser.add_argument("outfile", help=".cool output file")
        parser.add_argument("-r", "--resolution",
            help="integer bp resolution desired in cooler file. "
                 "Setting to 0 (default) will use all resolutions. ",
                 type=int, default=0)
        parser.add_argument("-e", "--exclude_MT",
            help="if used, exclude the mitochondria (MT) from the output",
            action="store_true")
        args = parser.parse_args()
        # these parameters adapted from theaidenlab/straw
        # KR is default normalization type and BP is the unit for binsize
        hic2cool_extractnorms(
            args.infile,
            args.outfile,
            args.resolution,
            #args.normalization,
            args.exclude_MT,
            True)
    main()
