__author__ = 'mactep'

from Bio import SeqIO


TRASH_FILE = "others.fasta"
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


def split_and_save(src, minlen, copy_path, trash_path):
    recs = SeqIO.parse(src, "fasta")
    maxcomprefix = None

    if not minlen:
        SeqIO.write(recs, copy_path, "fasta")
    else:
        c, t = [], []
        for rec in recs:
            if not maxcomprefix:
                maxcomprefix = rec.id
            else:
                pos = 0
                while pos < len(rec.id) and pos < len(maxcomprefix) and maxcomprefix[pos] == rec.id[pos]:
                    pos += 1
                maxcomprefix = maxcomprefix[:pos]
            if len(str(rec.seq).replace('-', '')) < minlen:
                t.append(rec)
            else:
                c.append(rec)
        SeqIO.write(c, copy_path, "fasta")
        SeqIO.write(t, trash_path, "fasta")
    return maxcomprefix