import sys
import argparse
from taxsbp import __version__


def main(arguments: str = None):

    if arguments is not None:
        sys.argv = arguments

    parser = argparse.ArgumentParser(
        prog="taxsbp", conflict_handler="resolve", add_help=True
    )
    parser.add_argument(
        "-i",
        "--input-file",
        required=True,
        help="Input file: unique_id <tab> weigth <tab> node",
    )
    parser.add_argument(
        "-o",
        "--output-file",
        help="Output file: id <tab> weigth <tab> node <tab> bin. Default: STDOUT",
    )
    parser.add_argument(
        "-t", "--taxonomy-file", help="Taxonomy file: node <tab> parent <tab> rank"
    )
    parser.add_argument(
        "-l",
        "--bin-len",
        type=int,
        help="Maximum bin length (in bp). Use this parameter insted of -b to define the number of bins. Default: length of the biggest group [Mutually exclusive -b]",
    )
    parser.add_argument(
        "-b",
        "--bins",
        type=int,
        help="Approximate number of bins (estimated by total length/bin number). [Mutually exclusive -l]",
    )
    parser.add_argument(
        "-p",
        "--pre-cluster",
        type=str,
        default="",
        help="",
    )
    parser.add_argument(
        "-e",
        "--bin-exclusive",
        type=str,
        default="",
        help="",
    )
    parser.add_argument(
        "-s",
        "--silent",
        default=False,
        action="store_true",
        help="Ignore warning",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="version: %(prog)s " + __version__,
        help="Show program's version number and exit.",
    )

    if len(sys.argv) == 1:  # Print help calling script without parameters
        parser.print_help()
        return False

    args = parser.parse_args()  # read sys.argv[1:] by default

    print(args)


if __name__ == "__main__":
    sys.exit(0 if main() else 1)
