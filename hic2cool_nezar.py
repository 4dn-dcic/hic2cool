from __future__ import division, print_function
import struct
import mmap
import zlib

import numpy as np
import pandas as pd
import cooler


def readcstr(f):
    """
    Read byte-by-byte until encountering a null char. Returns decoded string.

    """
    buf = b""
    while True:
        b = f.read(1)
        if b is None or b == b"\0":
            return buf.decode("utf-8", errors="ignore")
        else:
            buf += b


def read_header(f):
    magic_string = struct.unpack('<3s', f.read(3))[0]
    f.read(1)
    if (magic_string != b"HIC"):
        raise ValueError('This does not appear to be a HiC file; '
              'magic string is incorrect')

    # version and master index location
    version = struct.unpack('<i', f.read(4))[0]
    masterindex = struct.unpack('<q', f.read(8))[0]

    # genome assembly name
    genome = b""
    c = f.read(1)
    while (c != b'\0'):
        genome += c
        c = f.read(1)
    genome = genome.decode('ascii')

    # metadata
    metadata = {}
    nattributes = struct.unpack('<i', f.read(4))[0]
    for _ in range(nattributes):
        key = readcstr(f)
        value = readcstr(f)
        metadata[key] = value

    # chromosomes
    chrs = {}
    nChrs = struct.unpack('<i', f.read(4))[0]
    for i in range(0, nChrs):
        name = readcstr(f)
        length = struct.unpack('<i', f.read(4))[0]
        if name and length:
            chrs[i] = [i, name, length]

    # bp delimited resolutions
    resolutions = []
    nBpRes = struct.unpack('<i', f.read(4))[0]
    for _ in range(0, nBpRes):
        res = struct.unpack('<i', f.read(4))[0]
        resolutions.append(res)

    return version, genome, metadata, masterindex, chrs, resolutions


def read_footer(f, buf, masterindex):
    f.seek(masterindex)

    cpair_info = {}
    nBytes = struct.unpack('<i', f.read(4))[0]
    nEntries = struct.unpack('<i', f.read(4))[0]
    for _ in range(nEntries):
        key = readcstr(f)
        fpos = struct.unpack('<q', f.read(8))[0]
        sizeinbytes = struct.unpack('<i', f.read(4))[0]
        cpair_info[key] = fpos

    expected = {}
    factors = {}
    # raw
    nExpectedValues = struct.unpack('<i', f.read(4))[0]
    for _ in range(nExpectedValues):
        unit = readcstr(f)
        binsize = struct.unpack('<i', f.read(4))[0]

        nValues = struct.unpack('<i', f.read(4))[0]
        expected['RAW', unit, binsize] = np.frombuffer(
            buf,
            dtype=np.dtype('<d'),
            count=nValues,
            offset=f.tell())
        f.seek(nValues * 8, 1)

        nNormalizationFactors = struct.unpack('<i', f.read(4))[0]
        factors['RAW', unit, binsize] = np.frombuffer(
            buf,
            dtype=np.dtype([('chrom', '<i'), ('factor', '<d')]),
            count=nNormalizationFactors,
            offset=f.tell())
        f.seek(nNormalizationFactors * 12, 1)
    # normalized
    nExpectedValues = struct.unpack('<i', f.read(4))[0]
    for _ in range(nExpectedValues):
        normtype = readcstr(f)
        unit = readcstr(f)
        binsize = struct.unpack('<i', f.read(4))[0]

        nValues = struct.unpack('<i', f.read(4))[0]
        expected[normtype, unit, binsize] = np.frombuffer(
            buf,
            dtype='<d',
            count=nValues,
            offset=f.tell())
        f.seek(nValues * 8, 1)

        nNormalizationFactors = struct.unpack('<i', f.read(4))[0]
        factors[normtype, unit, binsize] = np.frombuffer(
            buf,
            dtype=np.dtype([('chrom', '<i'), ('factor', '<d')]),
            count=nNormalizationFactors,
            offset=f.tell())
        f.seek(nNormalizationFactors * 12, 1)

    norm_info = {}
    nEntries = struct.unpack('<i', f.read(4))[0]
    for _ in range(nEntries):
        normtype = readcstr(f)
        chrIdx = struct.unpack('<i', f.read(4))[0]
        unit = readcstr(f)
        resolution = struct.unpack('<i', f.read(4))[0]
        filePosition = struct.unpack('<q', f.read(8))[0]
        sizeInBytes = struct.unpack('<i', f.read(4))[0]
        norm_info[normtype, unit, resolution, chrIdx] = {
            'filepos': filePosition,
            'size': sizeInBytes
        }

    return cpair_info, expected, factors, norm_info


def read_blockinfo(f, buf, filepos, unit, binsize):
    """
    Read the locations of the data blocks of a cpair

    """
    f.seek(filepos)

    c1 = struct.unpack('<i', f.read(4))[0]
    c2 = struct.unpack('<i', f.read(4))[0]
    nRes = struct.unpack('<i', f.read(4))[0]
    for _ in range(nRes):
        unit_ = readcstr(f)
        f.seek(5 * 4, 1)
        binsize_ = struct.unpack('<i', f.read(4))[0]
        blockBinCount = struct.unpack('<i', f.read(4))[0]
        blockColCount = struct.unpack('<i', f.read(4))[0]
        nBlocks = struct.unpack('<i', f.read(4))[0]
        if unit == unit_ and binsize == binsize_:
            break
        f.seek(nBlocks * 16, 1)

    dt = np.dtype([('block_id', '<i'), ('filepos', '<q'), ('size', '<i')])
    data = np.frombuffer(buf, dt, count=nBlocks, offset=f.tell())
    return data, blockBinCount, blockColCount


def read_normvector(f, buf, filepos):
    f.seek(filepos)
    nValues = struct.unpack('<i', f.read(4))[0]
    return np.frombuffer(buf, dtype=np.dtype('<d'), count=nValues, offset=filepos+4)


def read_block(f, filepos, size):
    f.seek(filepos)
    bytes_ = zlib.decompress(f.read(size))

    nRecords = struct.unpack('<i', bytes_[0:4])[0]
    binXOffset = struct.unpack('<i', bytes_[4:8])[0]
    binYOffset = struct.unpack('<i', bytes_[8:12])[0]
    useShort = (struct.unpack('<b', bytes_[12:13])[0] == 0)
    type_ = struct.unpack('<b', bytes_[13:14])[0]

    block = np.zeros(
        nRecords,
        dtype=[
            ('binX', '<i'),
            ('binY', '<i'),
            ('counts', '<i') # '<h' if useShort else '<i')
        ]
    )
    binX = block['binX']
    binY = block['binY']
    counts = block['counts']

    if type_ == 1:
        k = 0
        ptr = 14
        rowCount = struct.unpack('<h', bytes_[ptr:ptr+2])[0]; ptr += 2
        for i in range(rowCount):
            y = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
            colCount = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
            y += binYOffset
            for j in range(colCount):
                x = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
                x += binXOffset
                if useShort:
                    c = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
                else:
                    c = struct.unpack('<f', bytes_[ptr:(ptr+4)])[0]; ptr += 4
                binX[k] = x
                binY[k] = y
                counts[k] = c
                if c < 0 or counts[k] < 0:
                    print(c, counts[k])
                k += 1
    elif type_ ==  2:
        k = 0
        ptr = 14
        nPts = struct.unpack('<i', bytes_[ptr:(ptr+4)])[0]; ptr += 4
        w = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
        for i in range(nPts):
            row = int(i / w)
            col = i - row * w
            x = int(binXOffset + col)
            y = int(binYOffset + row)
            if useShort:
                c = struct.unpack('<h', bytes_[ptr:(ptr+2)])[0]; ptr += 2
                if (c != -32768):
                    binX[k] = x
                    binY[k] = y
                    counts[k] = c
                    k += 1
            else:
                c = struct.unpack('<f',bytes_[ptr:(ptr+4)])[0]; ptr += 4
                if (c != 0x7fc00000):
                    binX[k] = x
                    binY[k] = y
                    counts[k] = c
                    k += 1
    else:
        raise ValueError("Unknown block layout type")

    del bytes_
    return block


def fetch_region(f, blockstats,
                 chromX_offset,
                 chromY_offset,
                 binX0, binX1,
                 binY0, binY1):

    blkBinCount = blockstats['binCount']
    blkColCount = blockstats['colCount']
    blkinfo = blockstats['info']
    block_ids = set(blkinfo['block_id'])
    i0 = binX0 // blkBinCount
    i1 = binX1 // blkBinCount + 1
    j0 = binY0 // blkBinCount
    j1 = binY1 // blkBinCount + 1

    query_block_ids = set()
    for r in range(i0, i1):
        for c in range(j0, j1):
            query_block_ids.add(r + blkColCount * c)

    if chromX_offset == chromY_offset:
        for r in range(j0, j1):
            for c in range(i0, i1):
                query_block_ids.add(r + blkColCount * c)

    result = []
    for blk_id, filepos, size in blkinfo:
        if blk_id in query_block_ids:
            data = read_block(f, int(filepos), int(size))
            mask = ((data['binX'] >= binX0) &
                    (data['binX'] < binX1) &
                    (data['binY'] >= binY0) &
                    (data['binY'] < binY1))
            data = data[mask]
            data['binX'] += chromX_offset
            data['binY'] += chromY_offset
            result.append(data)
    return result


def make_bins(f, buf, chrs, binsize, norm_info):
    unit = 'BP'
    chrs = [chrs[i] for i in sorted(chrs.keys()) if chrs[i][1].upper() != 'ALL']
    chromsizes = pd.Series(
        index=[c[1] for c in chrs],
        data=[c[2] for c in chrs])
    chrom_ids = [c[0] for c in chrs]

    bins = cooler.binnify(chromsizes, binsize)
    grouped = bins.groupby('chrom')
    nbins_per_chrom = {c[0]: len(grouped.get_group(c[1]))
                           for c in chrs}
    for norm in ["VC", "VC_SQRT", "KR"]:
        vector = np.concatenate([
            read_normvector(
                f,
                buf,
                norm_info[norm, unit, binsize, cid]['filepos'])[:nbins_per_chrom[cid]]
            for cid in chrom_ids
        ])
        bins[norm] = vector
    return bins


def iter_pixels(f, buf, chrs, binsize, cpair_info, bins):
    unit = 'BP'
    chrs = [chrs[i] for i in sorted(chrs.keys()) if chrs[i][1].upper() != 'ALL']
    chroms = [c[1] for c in chrs]
    chromsizes = pd.Series(
        index=[c[1] for c in chrs],
        data=[c[2] for c in chrs])
    idmap = {c[1]: c[0] for c in chrs}
    chrom_sizes_in_bins = bins.groupby('chrom', sort=False).size()
    x = bins['chrom'].values
    chrom_offsets = np.r_[0, np.flatnonzero(x[1:] != x[:-1]) + 1]

    blockstats = {}
    for i, chromX in enumerate(chroms):
        for chromY in chroms[i:]:
            chr_key = '{}_{}'.format(idmap[chromX], idmap[chromY])
            blkinfo, bc, cc = read_blockinfo(
                f, buf, cpair_info[chr_key], unit, binsize)
            blockstats[chromX, chromY] = {
                'binCount': bc,
                'colCount': cc,
                'info': blkinfo,
            }

    for i in range(len(chroms)):
        chromX = chroms[i]
        chromX_offset = chrom_offsets[i]

        step = max(blockstats[chromX, y]['binCount'] for y in chroms[i:])
        n_rows = chrom_sizes_in_bins[i] // step + 1
        print(chromX, '#rows of blocks = {}'.format(n_rows))
        print('step size (bins)', step)

        for r in range(n_rows):
            row_of_blocks = []

            for j in range(i, len(chroms)):
                chromY = chroms[j]
                chromY_offset = chrom_offsets[j]

                row_of_blocks.extend(
                    fetch_region(
                        f,
                        blockstats[chromX, chromY],
                        chromX_offset,
                        chromY_offset,
                        r * step,
                        (r+1) * step,
                        0,
                        chrom_sizes_in_bins[j])
                )

            if len(row_of_blocks):
                out = np.concatenate(row_of_blocks, axis=0)
                out.sort(order=['binX', 'binY'])
                print("Finished row {}".format(r))
                yield {
                    'bin1_id': out['binX'],
                    'bin2_id': out['binY'],
                    'count': out['counts'],
                }
                del row_of_blocks
                del out



if __name__ == '__main__':

    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "infile",
        help=".hic input file")
    parser.add_argument(
        "outfile",
        help=".cool output file")
    parser.add_argument(
        "-r", "--resolution",
        type=int,
        default=5000)
    args = vars(parser.parse_args())

    f = open(args['infile'], 'rb')
    buf = mmap.mmap(f.fileno(), 0, access=mmap.ACCESS_READ)
    version, genome, metadata, masterindex, chrs, resolutions = read_header(f)
    cpair_info, expected, factors, norm_info = read_footer(f, buf, masterindex)

    binsize = args['resolution']
    bins = make_bins(f, buf, chrs, binsize, norm_info)
    pixels = iter_pixels(f, buf, chrs, binsize, cpair_info, bins)
    cooler.io.create(
        args['outfile'],
        bins,
        pixels,
        metadata=metadata,
        assembly=genome)
