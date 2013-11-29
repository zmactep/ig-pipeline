__author__ = 'pavel'

import json
import argparse


def merge_json(target, groups):
    t = json.loads(target)
    g = json.loads(groups)['groups']

    for group in t['groups']:
        d = {}
        for name in t['groups'][group]:
            d[name] = g[name]
        t['groups'][group] = d

    return json.dumps(t)


def main():
    parser = argparse.ArgumentParser(description="Merge info json")
    parser.add_argument('--next', action='store', help='current json')
    parser.add_argument('--prev', action='store', help='pred json with info about groups')
    parser.add_argument('--out', action='store', help='output json')

    args, unknown = parser.parse_known_args()
    with open(args.out, "wt") as fd:
        fd.write(merge_json(open(args.next, "rt").readline(), open(args.prev, "rt").readline()))


if __name__ == "__main__":
    main()