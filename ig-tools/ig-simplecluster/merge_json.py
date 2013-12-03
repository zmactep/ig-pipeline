__author__ = 'pavel'

import sys
import os
import json
import logging
import argparse


def merge_json(target, groups):
    try:
        t = json.loads(target)
    except:
        logging.error("Target JSON file corrupted.")
        raise
    try:
        g = json.loads(groups)['groups']
    except:
        logging.error("Pre-level JSON file corrupted.")
        raise

    try:
        for group in t['groups']:
            d = {}
            for name in t['groups'][group]:
                d[name] = g[name]
            t['groups'][group] = d
    except:
        logging.error("JSON structure is wrong.")
        raise

    return json.dumps(t)


def main():
    parser = argparse.ArgumentParser(description="Merge info json")
    parser.add_argument('--next', action='store', help='current json')
    parser.add_argument('--prev', action='store', help='pred json with info about groups')
    parser.add_argument('--out', action='store', help='output json')

    args, unknown = parser.parse_known_args()

    if not os.path.isfile(args.next) or not args.next:
        logging.error("Next file does not exist.")
        raise FileExistsError
    if not os.path.isfile(args.prev) or not args.prev:
        logging.error("Prev file does not exist.")
        raise FileExistsError

    if not args.out:
        logging.error("Output file is not specified.")
        raise FileExistsError

    logging.info("Started with right input parameters.")
    logging.info("Next: %s, Prev: %s, Out: %s." % (args.next, args.prev, args.out))

    with open(args.next, "rt") as fd_next, open(args.prev, "rt") as fd_prev:
        result = merge_json(fd_next.readline(), fd_prev.readline())

        with open(args.out, "wt") as fd:
            fd.write(result)


if __name__ == "__main__":
    FORMAT = "[%(levelname)s] %(asctime)s | %(message)s"
    logging.basicConfig(format=FORMAT, level=logging.INFO, stream=sys.stderr)
    try:
        main()
    except:
        logging.critical("Terminated.")