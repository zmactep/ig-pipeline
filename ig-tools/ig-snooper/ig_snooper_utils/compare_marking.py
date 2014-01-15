from __future__ import division

__author__ = 'pavel'

import argparse
import json


def compare_marking(marking1, marking2, part, chothia):
    hb = open(marking1, "rt")
    hp = open(marking2, "rt")

    up = 10 if part else 14

    total = 0
    total_error = 0
    seq_error = 0
    misfit_error = 0
    for i, j in zip(hb, hp):
        total += 1
        if i == j:
            continue
        line_b = [int(k) for k in i.strip().split('\t')[1:]]
        line_p = [int(k) for k in j.strip().split('\t')[1:]]
        if chothia == 1:
            line_b[1] += 5
            line_b[2] += 5
        elif chothia == 2:
            line_p[1] += 5
            line_p[2] += 5
        if line_b[:up] != line_p[:up]:
            total_error += 1
            for o, (k, l) in enumerate(list(zip(line_p, line_b))[:up]):
                seq_error += int(k != l)
                misfit_error += k - l if o % 2 == 0 else l - k
    return total_error, seq_error, misfit_error, total


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('marking1', help='input marking [igblast]')
    parser.add_argument('marking2', help='input marking [igblast]')
    parser.add_argument('--part', action="store_const", default=False, const=True,
                        help='compare just till cdr3 (default: false)')
    parser.add_argument('--json', help='json dump output')
    parser.add_argument('--chothia', metavar='N', type=int, help='use chothia for markingN')
    return parser.parse_args()


def main():
    args = parse_args()
    total_error, seq_error,  misfit_error, total = compare_marking(args.marking1, args.marking2, args.part, args.chothia)
    d = {
        'wrong_rate': (total_error / total) if total else 0,
        'error_rate': (seq_error / total) if total else 0,
        'misfit_rate': (misfit_error / total) if total else 0,
        'error_wrong_rate': (seq_error / total_error) if total_error else 0,
        'misfit_wrong_rate': (misfit_error / total_error) if total_error else 0
    }
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
