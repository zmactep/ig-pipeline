import numpy as np
import align

from Bio.Seq import Seq
from Bio import SeqIO


buckets_cache = {}

last_rec_id = None


mids = {
    'l-vh': ('TATGATACGC', 'TCACTCATAC'),
    'l-vk': ('TCATCGAGTC', 'TCGAGCTCTC'),
    'l-vl': ('TCGCAGACAC', 'TCTGTCTCGC'),
    'm-vk': ('CGACACTATC', 'CGAGACGCGC'),
    'm-vkh': ('TAGCTCTATC', 'TATAGACATC'),
    'm-vl': ('CGTATGCGAC', 'CGTCGATCTC'),
    'm-vlh': ('TAGACTGCAC', 'TAGCGCGCGC')
}


primers = {
    'l-vh': ('SARSTKCAGCTGSWGSAGTCTGG', 'CCTGARGAGACGGTGACCAKKGTYCC'),
    'l-vk': ('GAYRTYSWGWTGACYCAGWCTCC', 'TTTRATHTCCASYYKKGTCCC'),
    'l-vl': ('CWGBYTGTGYTGACKCARCC', 'CTGKARRACGGTSASCTKGGTCCC'),
    'm-vh': ('GCCTACGGCAGCCGCTGGATTGTTATTAC', 'CACAGACGGGCCTTTTGTAGAC'),
    'm-vk': ('GCTGCTGCTGGTCTGCTGCTCCTCGCTG', 'GGCGGGAAAATAAAAACAGACGG'),
    'm-vl': ('GCTGCTGCTGGTCTGCTGCTCCTCGCTG', 'GGCGGGAAAATAAAAACAGACGG')
}


def get_rc(rec):
    return SeqIO.SeqRecord(rec.seq.reverse_complement(), 
			   rec.id + " (RC)", rec.name + " (RC)", 
			   rec.description, 
			   letter_annotations={"phred_quality": rec.letter_annotations["phred_quality"][::-1]})


def align_test(s, m):
    S = np.array([[1 if i == j else -1 for i in range(256)]
                  for j in range(256)], dtype=np.short)
    score, p, _ = align.align(list(align.string_to_alignment(s)),
                              list(align.string_to_alignment(m)),
                              -1, -1, S, True, True)
    p = align.alignment_to_string(p).replace('-', '')
    return score, s.find(p)


def split_dataset(recs, mids, ms=12):
    def test(seq, f, r, t_max):
        fscore, fpos = align_test(seq, f)
        rscore, rpos = align_test(seq, r)
        if fpos < rpos and fscore + rscore >= ms and fscore + rscore >= t_max:
            return fscore + rscore
  
    buckets = {"l-vl": [], "l-vk": [], "l-vh": [],
               "m-vk": [], "m-vl": [], "m-vlh": [], "m-vkh": [],
               "other": [], "trash": []}
    for i, rec in enumerate(recs):
        if len(rec) > 550 or \
           len(rec) < 250 or \
           np.average(rec.letter_annotations["phred_quality"]) < 21:
            buckets["trash"].append(rec)
            continue
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
        if i % 1000 == 0:
            print({mid:len(buckets[mid]) for mid in buckets})
            buckets_cache = buckets
            last_rec_id = rec.id
    return buckets
