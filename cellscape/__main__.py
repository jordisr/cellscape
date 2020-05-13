import argparse, sys, os

def make_cartoon(args):
    print("test")

def main():
    # set up argument parser
    parser = argparse.ArgumentParser(description='CellScape: Vector graphics for molecular visualization')
    subparsers = parser.add_subparsers(dest="command")
    subparsers.required=True

    # cartoon (formerly pdb2svg.py)
    parser_cartoon = subparsers.add_parser('cartoon', help='', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser_cartoon.set_defaults(func=make_cartoon)
    parser_cartoon.add_argument('--test', help='')

    # parse arguments and call corresponding command
    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    main()
