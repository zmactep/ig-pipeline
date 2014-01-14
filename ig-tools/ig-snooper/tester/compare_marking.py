from __future__ import division

import argparse
import json


def parse_marking_line(line):
    elems = line.split()[1:]
    bounds = [(int(elems[2 * i]) - 1, int(elems[2 * i + 1])) for i in range(7)]
    return bounds


def get_marking_dict(marking):
    with open(marking, "rt") as fd:
        marking_strings = fd.readlines()
    return {line.split()[0]: parse_marking_line(line) for line in marking_strings}


def compare_marking(mdict1, mdict2):
    total_error = 0
    total = 0
    seq_error = 0
    misfit_error = 0
    for name in mdict1:
        if name not in mdict2:
            continue
        total += 1
        if mdict1[name] != mdict2[name]:
            total_error += 1
            for (i, j) in zip(mdict1[name], mdict2[name]):
                seq_error += int(i[0] != j[0]) + int(i[1] != j[0])
                misfit_error += abs((i[1] - i[0]) - (j[1] - j[0]))
    return {'wrong_rate': (total_error / total), 'error_rate': (seq_error / total),
            'misfit_rate': (misfit_error / total),
            'error_wrong_rate': (seq_error / total_error), 'misfit_wrong_rate': (misfit_error / total_error)}


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('marking1', help='input marking [fasta]')
    parser.add_argument('marking2', help='input marking [igblast]')
    parser.add_argument('--json', help='json dump output')
    return parser.parse_args()


def main():
    args = parse_args()
    d = compare_marking(get_marking_dict(args.marking1), get_marking_dict(args.marking2))
    print("Wrong sequence annotation rate: %.4f" % d['wrong_rate'])
    print("Average error rate per sequence: %.4f" % d['error_rate'])
    print("Average region misfit per sequence: %.4f" % d['misfit_rate'])
    print("Average error rate per wrong sequence: %.4f" % d['error_wrong_rate'])
    print("Average region misfit per wrong sequence: %.4f" % d['misfit_wrong_rate'])
    if args.json:
        with open(args.json, "wt") as fd:
            fd.write(json.dumps(d))


if __name__ == "__main__":
    main()
