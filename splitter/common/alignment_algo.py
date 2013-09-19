import numpy as np
import align
from Bio.Seq import Seq

from common.enums import Direction

S = np.array([[1 if i == j else -1
               for i in range(256)] for j in range(256)], dtype=np.short)


def align_shortlist(seq, shortlist):
    result = []

    nseq = list(align.string_to_alignment(seq))
    for sh in shortlist:
        score, nr, shr = align.align(nseq, list(align.string_to_alignment(sh)),
                                     -1, -1, S, True, True)
        result.append((align.alignment_to_string(nr),
                       align.alignment_to_string(shr)))
    return result


def cut_local_ends(seq, ends):
    ends_ng = [end.replace("-", "") for end in ends]
    s = seq.find(ends_ng[0]) + len(ends_ng[0])
    e = seq.find(ends_ng[1])
    return seq[s:e]


def cut_record(rec, ends):
    seq = cut_local_ends(str(rec.seq), ends)
    r = rec
    r.seq = Seq(seq)
    return r


def cut_record_by_dir(rec, direction, c):
    seq = str(rec.seq)

    if direction == Direction.forward:
        rfs, _ = zip(*align_shortlist(seq, c.f))
        return cut_record(rec, rfs)
    elif direction == Direction.reversed:
        rrs, _ = zip(*align_shortlist(seq, c.r))
        recv = cut_record(rec, rrs).reverse_complement()
        recv.id = rec.id + " (RC)"
        recv.description = rec.description
        return recv
    assert "Wrong direction"


def local_check(*zipped):
    for x, y in list(*zipped):
        nx = list(align.string_to_alignment(x.replace("-", "")))
        ny = list(align.string_to_alignment(y))
        score, xr, yr = align.align(nx, ny, -1, -1, S, True, True)
        xr = align.alignment_to_string(xr)
        yr = align.alignment_to_string(yr)
        # diff = abs(len(y.replace("-", "")) - len(xr.replace("-", "")))
        # if diff > 1:
        #     return False
        if hamming_dist(y, xr) > 1:
            # print(y, xr, yr, hamming_dist(y, xr))
            return False
    return True


def hamming_dist(s, t):
    if len(s) != len(t):
        return float("inf")
    return sum(s[i] != t[i] for i in range(len(s)))
