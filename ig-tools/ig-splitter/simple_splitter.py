import json
import tempfile
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import align

def align_test(s, m):
    S = np.array([[1 if i == j else -1 for i in range(256)] for j in range(256)], dtype=np.short)
    score, p, _ = align.align(list(align.string_to_alignment(s)),
                              list(align.string_to_alignment(m)),
                              -1, -1, S, True, True)
    p = align.alignment_to_string(p).replace('-', '')
    return score, s.find(p)


def get_rc(rec):
    return SeqIO.SeqRecord(rec.seq.reverse_complement(),
                           rec.id + " (RC)", rec.name + " (RC)",
                           rec.description,
                           letter_annotations={"phred_quality": rec.letter_annotations["phred_quality"][::-1]})


def test(seq, f, r, t_max, ms=12):
        fscore, fpos = align_test(seq, f)
        rscore, rpos = align_test(seq, r)
        if fpos < rpos and fscore + rscore >= ms and fscore + rscore >= t_max:
            return fscore + rscore


def load_mids():
    try:
        with open('config/simple_config.json', 'r') as conf:
            return json.load(conf)['mids']
    except (OSError, ValueError) as e:
        raise RuntimeError("Error reading configuration") from e


def split_dataset(recs):
    class Split:
        def __init__(self, mids):
            self.buckets = {k: [] for k in list(mids) + ['trash', 'other']}

        def __enter__(self):
            self.fasta = tempfile.NamedTemporaryFile(mode='w', suffix='.fasta')

        def __exit__(self, type, value, tb):
            self.fasta.close()

        def accept(self, key, val):
            self.buckets[key].append(val)
            SeqIO.write(val, self.fasta, 'fasta')

    mids = load_mids()

    with Split(mids) as split:
        for rec in recs:
            if len(rec) not in range(250, 551) or np.average(rec.letter_annotations["phred_quality"]) < 21:
                split.accept("trash", rec)
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
                split.accept(t_mid, rec if t_dir else get_rc(rec))
            else:
                split.accept('other', rec)

        return split.buckets