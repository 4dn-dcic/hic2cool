from __future__ import absolute_import
from . import hic2cool_convert, hic2cool_update, __version__
import argparse
import sys


def main():
    """
    Execute the program from the command line
    """
    # the primary parser is used for hic2cool -v or -h
    primary_parser = argparse.ArgumentParser(prog='hic2cool', add_help=False)
    primary_parser.add_argument('-v', '--version', action='version',
                                version='%(prog)s ' + __version__)
    # the secondary parser is used for the specific run mode
    secondary_parser = argparse.ArgumentParser(prog='hic2cool', parents=[primary_parser])
    # the subparsers collect the args used to run the hic2cool mode
    subparsers = secondary_parser.add_subparsers(
        title='program modes',
        description='choose one of the following modes to run hic2cool:',
        dest='mode',
        metavar='mode: {convert, update}'
    )
    subparsers.required = True
    # add a subparser for the 'convert' command
    convert_help = 'convert a hic file to a cooler file'
    convert_subparser = subparsers.add_parser('convert', help=convert_help, description=convert_help)
    convert_subparser.add_argument("infile", help="hic input file path")
    convert_subparser.add_argument("outfile", help="cooler output file path")
    convert_subparser.add_argument(
        "-r", "--resolution",
        help="integer bp resolution desired in cooler file. Setting to 0 "
             "(default) will use all resolutions. If all resolutions are used, "
             "a multi-res .cool file will be created, which has a different "
             "hdf5 structure. See the README for more info",
        type=int,
        default=0
    )
    # add a subparser for the 'update' command
    update_help = 'update a cooler file produced by hic2cool'
    update_subparser = subparsers.add_parser('update', help=update_help, description=update_help)
    update_subparser.add_argument("infile", help="cooler input file path")
    update_subparser.add_argument("-o", "--outfile", help="optional new output file path", default='')
    # arguments shared by all subparsers
    for sp in [convert_subparser, update_subparser]:
        sp.add_argument(
            "-w", "--warnings",
            help="if used, print out non-critical WARNING messages, which are "
                 "hidden by default.",
            action="store_true"
        )

    # two step argument parsing
    # first check for top level -v or -h (i.e. `hic2cool -v`)
    (primary_namespace, remaining) = primary_parser.parse_known_args()
    # lastly, get the mode and specific args
    args = secondary_parser.parse_args(args=remaining, namespace=primary_namespace)
    if args.mode == 'convert':
        hic2cool_convert(args.infile, args.outfile, args.resolution, args.warnings, True)
    elif args.mode == 'update':
        hic2cool_update(args.infile, args.outfile, args.warnings, True)


if __name__ == '__main__':
    main()
