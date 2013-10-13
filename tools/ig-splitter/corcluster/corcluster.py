from confing import config as cfg
from common.small_algo import kmer_gen


def create_kstat(v_list, k):
    kstat = {}
    for i, rec in enumerate(v_list):
        for kmer in kmer_gen(str(rec.seq), k):
            kstat[kmer].append(i)
    return kstat


def create_tree(v_list, kstat):
    pass


def reduce_tree(cl):
    pass


def run(v_list):
    kst = create_kstat(v_list, cfg.cluster_kmer_length)
    cl = create_tree(v_list, kst)
    cl = reduce_tree(cl)
