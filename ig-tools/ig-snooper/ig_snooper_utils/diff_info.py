from __future__ import division

__author__ = 'pavel'

import argparse


def get_diff(marking1, marking2, part, chothia):
    up = 10 if part else 14

    diff = [0 for _ in range(up)]
    dsgn = [0 for _ in range(up)]
    dusg = [0 for _ in range(up)]

    hb = open(marking1, "rt")
    hp = open(marking2, "rt")

    error_list = []

    total = 0.
    for b, p in zip(hb, hp):
        total += 1
        if p == b:
            continue
        line_b = [int(k) for k in b.strip().split('\t')[1:up + 1]]
        line_p = [int(k) for k in p.strip().split('\t')[1:up + 1]]

        if chothia == 1:
            line_b[1] += 5
            line_b[2] += 5
        elif chothia == 2:
            line_p[1] += 5
            line_p[2] += 5

        f = False
        for i, j in zip(line_b, line_p):
            if abs(i - j) > 20:
                f = True
                break
        if f:
            error_list.append(b.strip().split('\t')[0])
            continue

        indicator = [int(i != j) for i, j in zip(line_b, line_p)]
        sndicator = [i - j for i, j in zip(line_b, line_p)]
        undicator = [abs(i - j) for i, j in zip(line_b, line_p)]
        for k, (v, s, u) in enumerate(zip(indicator, sndicator, undicator)):
            diff[k] += v
            dsgn[k] += s
            dusg[k] += u
    hb.close()
    hp.close()

    x = ["%s%d %s" % (r, n, se)
         for n in range(1, 5)
         for r in ["FR", "CDR"]
         for se in ["start", "end"]][:-6 if part else -2]

    stot = map(lambda i: i / total, dsgn)
    utot = map(lambda i: i / total, dusg)

    serr = map(lambda i: i[0] / i[1] if i[1] else 0, zip(dsgn, diff))
    uerr = map(lambda i: i[0] / i[1] if i[1] else 0, zip(dusg, diff))
	
    errr = map(lambda i: i / total, diff)

    return zip(x, zip(diff, stot, utot, serr, uerr, errr)), error_list


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('marking1', help='input marking [fasta]')
    parser.add_argument('marking2', help='input marking [igblast]')
    parser.add_argument('--chothia', metavar='N', type=int, help='use chothia for markingN')
    parser.add_argument('--part', action="store_const", default=False, const=True,
                        help='compare just till cdr3 (default: false)')
    return parser.parse_args()


def main():
    args = parse_args()
    d, error_list = get_diff(args.marking1, args.marking2, args.part, args.chothia)
    print("REGIONS:\terrors\t(signed\t/ unsigned)\t(signed\t/ unsigned)\t error rate")
    for name, (v, s, u, se, ue, er) in d:
        print("%s:\t%d\t(%.3f\t/ %.3f)\t(%.3f\t/ %.3f)\t%.4f" % (name, v, s, u, se, ue, er))
    print("Double-sequence list:")
    for e in error_list:
        print(e)


if __name__ == "__main__":
    main()
