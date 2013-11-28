from IPython.parallel import Client


buckets_cache = {}

last_rec_id = None

_ms = 12

mids = {
    'l-vh': ('TATGATACGC', 'TCACTCATAC'),
    'l-vk': ('TCATCGAGTC', 'TCGAGCTCTC'),
    'l-vl': ('TCGCAGACAC', 'TCTGTCTCGC'),
    'm-vk': ('CGACACTATC', 'CGAGACGCGC'),
    'm-vkh': ('TAGCTCTATC', 'TATAGACATC'),
    'm-vl': ('CGTATGCGAC', 'CGTCGATCTC'),
    'm-vlh': ('TAGACTGCAC', 'TAGCGCGCGC')
}

view = (Client()).load_balanced_view()


def align_test(s, m):
    import numpy as np
    import align
    S = np.array([[1 if i == j else -1 for i in range(256)] for j in range(256)], dtype=np.short)
    score, p, _ = align.align(list(align.string_to_alignment(s)),
                              list(align.string_to_alignment(m)),
                              -1, -1, S, True, True)
    p = align.alignment_to_string(p).replace('-', '')
    return score, s.find(p)


def test(seq, f, r, t_max):
    fscore, fpos = align_test(seq, f)
    rscore, rpos = align_test(seq, r)
    if fpos < rpos and fscore + rscore >= _ms and fscore + rscore >= t_max:
        return fscore + rscore


def get_rc(rec):
    from Bio import SeqIO
    return SeqIO.SeqRecord(rec.seq.reverse_complement(),
                           rec.id + " (RC)", rec.name + " (RC)",
                           rec.description,
                           letter_annotations={"phred_quality": rec.letter_annotations["phred_quality"][::-1]})

@view.parallel()
def process_record(rec):
    from Bio.Seq import Seq
    import numpy as np

    buckets = {"l-vl": [], "l-vk": [], "l-vh": [], "m-vk": [], "m-vl": [], "m-vlh": [], "m-vkh": [],
               "other": [], "trash": []}
    if len(rec) not in range(250, 551) or np.average(rec.letter_annotations["phred_quality"]) < 21:
        buckets["trash"].append(rec)
        return buckets
    seq = str(rec.seq).upper()
    rseq = str(rec.seq.reverse_complement).upper()
    t_pax = float("-inf")
    t_max = float("-inf")
    t_mid = None
    t_dir = True  # True - forward, False - reverse
    for mid in mids:
        f, r = mids[mid][0], str(Seq(mids[mid][1]).reverse_complement())
        # Forward test
        mx = test(seq, f, r, t_max)
        if mx:
            t_pax = t_max
            t_max, t_mid, t_dir = mx, mid, True
            # Reverse test
        mx = test(rseq, f, r, t_max)
        if mx:
            t_pax = t_max
            t_max, t_mid, t_dir = mx, mid, False
    if t_mid and t_max != t_pax:
        buckets[t_mid].append(rec if t_dir else get_rc(rec))
    else:
        buckets['other'].append(rec)

    return buckets


def merge_dicts(ds):
    m = {}

    for d in ds:
        for k, v in d.items():
            if k in m:
                m[k].extend(v)
            else:
                m[k] = v
    return m


def split_dataset(recs, ms=12):
    _ms = ms
    result = process_record.map(recs)

    result.wait_interactive()

    return merge_dicts(result.result)


    # todo log progress
    #    if i % 1000 == 0:
    #    print({mid: len(buckets[mid]) for mid in buckets})
    #    buckets_cache = buckets
    #    last_rec_id = rec.id


