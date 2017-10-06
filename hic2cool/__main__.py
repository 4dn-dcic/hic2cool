from __future__ import absolute_import
from . import hic2cool_convert
import argparse


def main():
    """
    Execute the program from the command line
    Args are:
    python -m hic2cool <infile (.hic)> <outfile (.cool)>
    <resolutions desired (defaults to all, optionally bp int)>
    <exclude MT (default False)>

    """
    parser = argparse.ArgumentParser('hic2cool')
    parser.add_argument(
        "infile",
        help=".hic input file")
    parser.add_argument(
        "outfile",
        help=".cool output file")
    parser.add_argument(
        "-r", "--resolution",
        help="integer bp resolution desired in cooler file. Setting to 0 "
             "(default) will use all resolutions. If all resolutions are used, "
             "a multi-res .cool file will be created, which has a different "
             "hdf5 structure. See the README for more info",
        type=int,
        default=0)
    parser.add_argument(
        "-e", "--exclude_missing",
        help="if used, allow for missing chromosome contacts",
        action="store_true")
    parser.add_argument(
        "-l", "--low_mem",
        help="if used, use less RAM at the expense of runtime",
        action="store_true")
    args = parser.parse_args()

    # these parameters adapted from theaidenlab/straw
    # KR is default normalization type and BP is the unit for binsize
    hic2cool_convert(args.infile, args.outfile, args.resolution,
                     args.exclude_missing, args.low_mem, True)


if __name__ == '__main__':
    main()
