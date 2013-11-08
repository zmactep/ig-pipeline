__author__ = 'mactep'

import argparse
from pair import report_generator


def main():
    parser = argparse.ArgumentParser(description="IG pairwise (VH-VL) report generator")
    parser.add_argument('hdir', action='store', help='clusterization result for heavy')
    parser.add_argument('ldir', action='store', help='clusterization result for light')
    parser.add_argument('--fix-suffix', action='store_const', dest='fix', default=False, const=True,
                        help="cut last '_' suffix")
    parser.add_argument('out', action='store', help='output json')

    args = parser.parse_args()
    report_generator.run(args.hdir, args.ldir, args.out, args.fix)


if __name__ == "__main__":
    main()