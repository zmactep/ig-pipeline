__author__ = 'mactep'

import os
import logging
from Bio import SeqIO
from Bio.Align.Applications import ClustalOmegaCommandline


TRASH_FILE = "others.fasta"
DATA_FILE = "data.fasta"
TREE_FILE = "cluster_tree.dnd"
ALIGNMENT_FILE = "alignment.fasta"
CONSENSUS_FILE = "consensus.fasta"
DENDOGRAM_FILE = "dendogram.png"
INFO_FILE = "info.json"

CLUSTER_SPLT = "G"
CLUSTER_MASK = "%s{}%i".format(CLUSTER_SPLT)
UTREE_MASK = "%s_unrooted.pdf"

TRASH_NAME = "other"

CLUSTER_DIR = "clusters"


def run_clustalo(src, tree_path, align_file):
    cline = ClustalOmegaCommandline(infile=os.path.abspath(src),
                                    guidetree_out=tree_path,
                                    outfmt="fasta", outfile=align_file,
                                    threads=4, force=True)
    logging.info("Aligner started")
    try:
        stdout, stderr = cline()
    except:
        logging.error("Aligner call error")
        raise

    if stdout:
        logging.error("Aligner internal error: %s" % stderr)
        raise RuntimeError
    return stdout, stderr


def split_and_save(src, minlen, copy_path, trash_path):
    logging.info("Splitting and filtering dataset.")
    try:
        recs = SeqIO.parse(src, "fasta")
    except:
        logging.error("Source file was not found in splitter action.")
        raise
    maxcomprefix = None

    c, t = [], []
    for rec in recs:
        if not maxcomprefix:
            maxcomprefix = rec.id
        else:
            pos = 0
            while pos < len(rec.id) and pos < len(maxcomprefix) and maxcomprefix[pos] == rec.id[pos]:
                pos += 1
            maxcomprefix = maxcomprefix[:pos]
        if minlen and len(str(rec.seq).replace('-', '')) < minlen:
            t.append(rec)
        else:
            c.append(rec)
    try:
        SeqIO.write(c, copy_path, "fasta")
        SeqIO.write(t, trash_path, "fasta")
    except:
        logging.error("Cannot write split step result.")
        raise
    if maxcomprefix.endswith(CLUSTER_SPLT):
         maxcomprefix = maxcomprefix[:-1]
    return maxcomprefix