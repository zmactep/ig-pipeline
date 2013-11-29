__author__ = 'mactep'

import os
import shutil
import json

from Bio import SeqIO, AlignIO
from Bio.Align.Applications import ClustalOmegaCommandline

import common


def run_clustalo(src, align_file):
    cline = ClustalOmegaCommandline(infile=os.path.abspath(src),
                                    outfmt="fasta", outfile=align_file,
                                    threads=4, force=True)
    stdout, stderr = cline()
    return stdout, stderr


def get_similars(align_file, m, skipf):
    alignment_dict = SeqIO.to_dict(AlignIO.read(align_file, "fasta"))
    sim_dict = {}
    for key in alignment_dict.keys():
        s = str(alignment_dict[key].seq)
        s = s[m:int(skipf * len(s))].replace('-', '')
        if s not in sim_dict:
            sim_dict[s] = []
        sim_dict[s].append(key)
    return alignment_dict, sim_dict


def save_clusters(abs_out, alignment_dict, sim_dict, mcp):
    cluster_path = os.path.join(abs_out, common.CLUSTER_DIR)
    if os.path.isdir(cluster_path):
        shutil.rmtree(cluster_path)
    os.mkdir(cluster_path)
    for cluster_id, seq in enumerate(sim_dict.keys()):
        to_write = []
        for key in sim_dict[seq]:
            to_write.append(alignment_dict[key])
        SeqIO.write(to_write, os.path.join(cluster_path, common.CLUSTER_MASK % (mcp, cluster_id + 1) + ".fasta"), "fasta")


def write_consensus(abs_out, alignment_dict, sim_dict, is_shortest, mcp):
    cons = []
    for cluster_id, seq in enumerate(sim_dict.keys()):
        val = None
        for key in sim_dict[seq]:
            if not val or \
               (is_shortest and len(val) > len(alignment_dict[key])) or \
               (not is_shortest and len(val) < len(alignment_dict[key])):
                val = alignment_dict[key].seq
        cons.append(SeqIO.SeqRecord(val, common.CLUSTER_MASK % (mcp, cluster_id + 1), name="", description=""))
    SeqIO.write(cons, os.path.join(abs_out, common.CONSENSUS_FILE), "fasta")


def write_info(abs_out, minlen, m, skipn, sim_dict, mcp, is_shortest):
    trash_size = len(SeqIO.to_dict(SeqIO.parse(os.path.join(abs_out, common.TRASH_FILE), "fasta")))
    with open(os.path.join(abs_out, common.INFO_FILE), "wt") as fd:
        fd.write(json.dumps({'minlen': minlen,
                             'head-skip': m,
                             'tail-prct': skipn,
                             'use-short': is_shortest,
                             "trash": trash_size,
                             "groups": {common.CLUSTER_MASK % (mcp, i + 1): sim_dict[d]
                                        for i, d in enumerate(sim_dict.keys())}}))


def run(src, out, minlen, m, skipn, is_shortest):
    abs_out = os.path.abspath(out)
    trash_path = os.path.join(abs_out, common.TRASH_FILE)
    copy_path = os.path.join(abs_out, os.path.basename(src))
    align_file = os.path.join(abs_out, common.ALIGNMENT_FILE)

    if not m:
        m = 0
    if not skipn:
        skipn = 100

    mcp = common.split_and_save(src, minlen, copy_path, trash_path)
    run_clustalo(copy_path, align_file)
    alignment_dict, sim_dict = get_similars(align_file, m, skipn / 100.0)
    save_clusters(abs_out, alignment_dict, sim_dict, mcp)
    write_consensus(abs_out, alignment_dict, sim_dict, is_shortest, mcp)
    write_info(abs_out, minlen, m, skipn, sim_dict, mcp, is_shortest)