import os
import re

from Bio.Seq import Seq
from Bio import SeqIO

from itertools import chain

from common.enums import DomainType, Direction
from splitter.configuration import SplitterConfiguration
from splitter.svm_split import svm_create, svm_split
from config import config as cfg
from common.alignment_algo import align_shortlist, local_check, cut_record
from common.small_algo import fasta_head


class DomainMarkers(object):
    def __init__(self, middict):
        mm, mp = middict['mid'], middict['primer']
        self.fm = (mm[Direction.forward],
                   str(Seq(mm[Direction.reversed]).reverse_complement()))
        self.rm = (mm[Direction.reversed],
                   str(Seq(mm[Direction.forward]).reverse_complement()))
        self.fstart = self.fm[0] + mp[Direction.forward]
        self.rstart = self.rm[0] + mp[Direction.reversed]
        self.f = (self.fstart, str(Seq(self.rstart).reverse_complement()))
        self.r = (self.rstart, str(Seq(self.fstart).reverse_complement()))


def load_configuration():
    sc = SplitterConfiguration()
    with open(cfg.getConf(cfg.splitter_configuration), "rt") as fd:
        lines = "".join(chain(i.strip() for i in fd.readlines()))
        sc.loads(lines)
    if not sc.mids:
        return
    sdict = {mid.d_type: {'mid': {Direction.forward: None,
                                  Direction.reversed: None},
                          'primer': {Direction.forward: None,
                                     Direction.reversed: None}}
             for mid in sc.mids}
    for mid in sc.mids:
        sdict[mid.d_type]['mid'][mid.direction] = mid.mid
        sdict[mid.d_type]['primer'][mid.direction] = mid.primer
    return {d_type: DomainMarkers(sdict[d_type]) for d_type in sdict}


def detect_chain(rec, c, additional_check=False):
    seq = str(rec.seq)
    cost_func = lambda X, Y: sum(len(x) - len(y.replace("-", ""))
                                 for x, y in zip(X, Y))

    rfs, rfl = zip(*align_shortlist(seq, c.f))
    if cost_func(c.f, rfl) < cfg.splitter_max_errors:
        if not additional_check or local_check(zip(rfs, c.fm)):
            return True, cut_record(rec, rfs)

    rrs, rrl = zip(*align_shortlist(seq, c.r))
    if cost_func(c.r, rrl) < cfg.splitter_max_errors:
        if not additional_check or local_check(zip(rrs, c.rm)):
            recv = cut_record(rec, rrs).reverse_complement()
            recv.id = rec.id + " (RC)"
            recv.description = rec.description
            return True, recv

    return False, rec


def mid_split(recs, config):
    srecs = {d_type: [] for d_type in config}
    unsplitted = []
    bad_reads = []
    for i, rec in enumerate(recs):
        if len(rec) < cfg.splitter_least_len:
            bad_reads.append(rec)
        else:
            for d_type in config:
                result, nrec = detect_chain(rec, config[d_type], True)
                if result:
                    srecs[d_type].append(nrec)
                    break
            else:
                unsplitted.append(rec)
        if not ((i + 1) % 100):
            tpl = tuple(map(len, srecs.values()))
            print("[%i %s / %i / %i]" % (sum(tpl), str(tpl),
                                         len(unsplitted),
                                         len(bad_reads)))
            total = sum(tpl) + len(unsplitted) + len(bad_reads)
            print(i + 1, total)
    return srecs, unsplitted, bad_reads


def remove_duplicates(v):
    rv = {}
    for d_type in v:
        rv[d_type] = list(set(v[d_type]))
    return rv


def merge_splits(v1, v2, rdup=True):
    v = {d_type: v1[d_type] for d_type in v1}
    for d_type in v2:
        if d_type not in v:
            v[d_type] = v2[d_type]
        else:
            for i in v2[d_type]:
                v[d_type].append(i)
    return remove_duplicates(v) if rdup else v


def pop_vh(v):
    vh = {DomainType.VH: v[DomainType.VH]}
    del v[DomainType.VH]
    return v, vh


def vh_h_split(vh):
    vh_h = {DomainType.VH: [], DomainType.VHH: []}
    pattern = re.compile(cfg.splitter_vhh_pattern)
    for d_type in vh:
        if d_type != DomainType.VH:
            continue
        for rec in vh[d_type]:
            seq = str(rec.seq)
            for i in range(2):
                pseq = str(Seq(seq[i:]).translate())
                if re.match(pattern, pseq):
                    vh_h[DomainType.VHH].append(rec)
                    break
            else:
                vh_h[DomainType.VH].append(rec)
    return vh_h


def length_filter(v):
    vn = {}
    for d_type in v:
        vn[d_type] = list(filter(lambda x: len(x) > cfg.splitter_least_len,
                            v[d_type]))
    return vn


def dump_fasta(results_dir, v):
    for d_type in v:
        SeqIO.write(v[d_type],
                    os.path.join(results_dir,
                                 cfg.splitter_outfasta.format(d_type.name)),
                    "fasta")


def print_info(v):
    tpl = tuple(map(len, v.values()))
    print("[%i %s]" % (sum(tpl), str(tpl)))


def run(fasta_source, results_dir):
    config = load_configuration()
    v, us, b = mid_split(SeqIO.parse(fasta_source, "fasta"), config)
    # v, us, b = mid_split(fasta_head(fasta_source, 2000), config)
    print_info(v)
    cl = svm_create(v)
    v2 = svm_split(cl, config, us)
    print_info(v2)
    v = merge_splits(v, v2)
    print_info(v)
    if cfg.splitter_vhh_pattern:
        v, vh = pop_vh(v)
        vh_h = vh_h_split(vh)
        v = merge_splits(v, vh_h, False)
    v = length_filter(v)
    print_info(v)
    dump_fasta(results_dir, v)
