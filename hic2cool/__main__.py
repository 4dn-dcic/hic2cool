from __future__ import absolute_import
from . import hic2cool_convert
import argparse


def main():
    """
    Execute the program from the command line
    Args are:
    python -m hic2cool <infile (.hic)> <outfile (.cool)> 
    <resolutions desired (defaults to all, optionally bp int)> 
    <normalization type (defaults to 'KR', optionally 'NONE', 'VC', or 'VC_SQRT')> 
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
        "-n", "--normalization", 
        help="string normalization type. Defaults to KR, optionally NONE, VC, "
             "or VC_SQRT",
        choices=['KR', 'NONE', 'VC', 'VC_SQRT'], 
        default='KR')
    parser.add_argument(
        "-e", "--exclude_MT", 
        help="if used, exclude the mitochondria (MT) from the output", 
        action="store_true")
    args = parser.parse_args()

    # these parameters adapted from theaidenlab/straw
    # KR is default normalization type and BP is the unit for binsize
    hic2cool_convert(args.infile, args.outfile, args.resolution, 
                     args.normalization, args.exclude_MT, True)


if __name__ == '__main__':
    main()
