import json
import os
from Bio import SeqIO
from Bio.Seq import Seq
import numpy as np
import align
from IPython.config import Application

# keeps list of dump file names from previous split runs, to be removed by a call to cleanup()
_dump_files = []

class Split:
    def get_dump_file(self, key):
        name = 'dump-{0}-{1}.fasta'.format(key, os.getpid())
        return open(name, 'w')

    def __init__(self, mids):
        keys = list(mids) + ['trash', 'other']
        self.__buckets = {k: [] for k in keys}
        self.__dumps = {k: self.get_dump_file(k) for k in keys}

    def __enter__(self):
        return self

    def __exit__(self, type, value, tb):
        for dump in self.__dumps.values():
            dump.close()

    def accept(self, key, val):
        self.__buckets[key].append(val)
        SeqIO.write(val, self.__dumps[key], 'fasta')

    def complete(self):
        return self.__buckets, [lambda name=f.name: os.remove(name) for f in self.__dumps.values()]


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


def load_mids():
    try:
        with open('config/simple_config.json', 'r') as conf:
            return json.load(conf)['mids']
    except (OSError, ValueError) as e:
        raise RuntimeError("Error reading configuration") from e


def split_dataset(recs):
    def test(seq, f, r, t_max, ms=12):
        fscore, fpos = align_test(seq, f)
        rscore, rpos = align_test(seq, r)
        if fpos < rpos and fscore + rscore >= ms and fscore + rscore >= t_max:
            return fscore + rscore

    mids = load_mids()

    with Split(mids) as split:
        print(split)

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

        buckets, dumps = split.complete()
        _dump_files.extend(dumps)
        return buckets


def cleanup():
    if _dump_files:
        for remover in _dump_files:
            remover()
