__author__ = 'mactep'

import os
import json
from itertools import chain

from Bio import SeqIO, AlignIO

import common


def fix_suffix_a(cls):
    names = list(chain(cls.values()))
    st = 1
    while len(list(filter(lambda name: names[0][-st:] == name[-st:], names))) == len(names):
        st += 1
    for key in cls.keys():
        cls[key] = list(map(lambda x: x[:-st+1], cls[key]))
    return cls


def fix_suffix(cls):
    for key in cls.keys():
        cls[key] = list(map(lambda x: x[:x.rfind('_')], cls[key]))
    return cls


def load_clusters(abs_src):
    clusters_dir = os.path.join(abs_src, common.CLUSTER_DIR)

    clusters = map(lambda x: os.path.join(clusters_dir, x),
                   filter(lambda x: x.endswith(".fasta"), os.listdir(clusters_dir)))
    trash = os.path.join(abs_src, common.TRASH_FILE)

    cls = {}
    for cluster in clusters:
        name = os.path.basename(cluster)
        name = name[:name.rfind('.')]
        cls[name] = list(map(lambda rec: rec.id,
                        SeqIO.parse(cluster, "fasta")))
    if os.path.isfile(trash):
        cls[common.TRASH_NAME] = list(map(lambda rec: rec.id,
                            SeqIO.parse(trash, "fasta")))
    else:
        cls[common.TRASH_NAME] = []

    return cls


def load_consensus(abs_src):
    return SeqIO.to_dict(AlignIO.read(os.path.join(abs_src, common.CONSENSUS_FILE), "fasta"))


def load_info(abs_src):
    info = os.path.join(abs_src, common.INFO_FILE)
    if os.path.isfile(info):
        return json.loads(open(info, "rt").readline())
    else:
        return {}