__author__ = 'mactep'

import os
import shutil
import itertools
import difflib
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

from Bio import Phylo, SeqIO, AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment
from Bio.Align.Applications import ClustalOmegaCommandline

import common


def distance(record1, record2):
    return 1 - difflib.SequenceMatcher(None, str(record1.seq),
                                       str(record2.seq)).ratio()


# Fails in current version of Biopython
def save_utree(tree_path):
    prefix_name = tree_path[:tree_path.rfind('.')]
    unrooted_tree = Phylo.read(tree_path, 'newick')
    unrooted_tree.ladderize()
    Phylo.draw_graphviz(unrooted_tree)
    plt.savefig(common.UTREE_MASK % prefix_name)


def make_distance_matrix(align_file):
    alignment_dict = SeqIO.to_dict(AlignIO.read(align_file, "fasta"))
    ids = dict(enumerate(alignment_dict.keys()))
    distance_matrix = np.zeros([len(ids)] * 2)
    for i, j in itertools.combinations(range(len(ids)), r=2):
        distance_matrix[i][j] = distance_matrix[j][i] = \
            distance(alignment_dict[ids[i]], alignment_dict[ids[j]])
    return distance_matrix, alignment_dict


def save_dendrogram(abs_out, distance_matrix, Y, cutoff):
    # Compute and plot dendrogram
    fig = plt.figure()
    axdendro = fig.add_axes([0.09, 0.1, 0.2, 0.8])
    Z = dendrogram(Y, orientation="right", color_threshold=cutoff)
    axdendro.set_yticks([])

    # Plot distance matrix
    axmatrix = fig.add_axes([0.3, 0.1, 0.6, 0.8])
    index = Z["leaves"]
    distance_matrix = distance_matrix[index, :]
    distance_matrix = distance_matrix[:, index]
    im = axmatrix.matshow(distance_matrix, aspect="auto", origin="lower")
    axmatrix.set_xticks([])
    axmatrix.set_yticks([])

    # Plot colorbar
    axcolor = fig.add_axes([0.91, 0.1, 0.02, 0.8])
    plt.colorbar(im, cax=axcolor)

    # Display and save figure
    dendogram_path = os.path.join(abs_out, common.DENDOGRAM_FILE)
    fig.savefig(dendogram_path)


def make_clusters(distance_matrix):
    Y = linkage(distance_matrix, method="centroid")
    cutoff = 0.5 * max(Y[:, 2])
    clusters = fcluster(Y, cutoff, "distance")

    return clusters, Y, cutoff


def make_fasta_clusters(clusters, alignment_dict):
    ids = dict(enumerate(alignment_dict.keys()))
    fasta_clusters = defaultdict(list)
    for i, cluster in enumerate(clusters):
        fasta_id = ids[i]
        fasta_clusters[cluster].append(alignment_dict[fasta_id])
    return fasta_clusters


def save_clusters(abs_out, fasta_clusters, mcp):
    cluster_path = os.path.join(abs_out, common.CLUSTER_DIR)
    if os.path.isdir(cluster_path):
        shutil.rmtree(cluster_path)
    os.mkdir(cluster_path)
    for cluster_id, cluster in fasta_clusters.items():
        SeqIO.write(cluster, os.path.join(cluster_path, common.CLUSTER_MASK % (mcp, cluster_id) + ".fasta"), 'fasta')


def write_consensus(abs_out, fasta_clusters, mcp):
    cons = list()
    for cluster_id, cluster in fasta_clusters.items():
        consensus = AlignInfo.SummaryInfo(MultipleSeqAlignment(cluster)).dumb_consensus()
        cons.append(SeqIO.SeqRecord(consensus, common.CLUSTER_MASK % (mcp, cluster_id), name="", description=""))
    SeqIO.write(cons, os.path.join(abs_out, common.CONSENSUS_FILE), "fasta")


def run_clustalo(src, tree_path, align_file):
    cline = ClustalOmegaCommandline(infile=os.path.abspath(src),
                                    guidetree_out=tree_path,
                                    outfmt="fasta", outfile=align_file,
                                    threads=4, force=True)
    stdout, stderr = cline()
    return stdout, stderr


def run(src, out, minlen=None):
    abs_out = os.path.abspath(out)
    trash_path = os.path.join(abs_out, common.TRASH_FILE)
    copy_path = os.path.join(abs_out, os.path.basename(src))
    tree_path = os.path.join(abs_out, common.TREE_FILE)
    align_file = os.path.join(abs_out, common.ALIGNMENT_FILE)

    mcp = common.split_and_save(src, minlen, copy_path, trash_path)
    run_clustalo(copy_path, tree_path, align_file)
    #save_utree(tree_path)
    distance_matrix, alignment_dict = make_distance_matrix(align_file)
    clusters, Y, cutoff = make_clusters(distance_matrix)
    save_dendrogram(abs_out, distance_matrix, Y, cutoff)
    fasta_clusters = make_fasta_clusters(clusters, alignment_dict)
    save_clusters(abs_out, fasta_clusters, mcp)
    write_consensus(abs_out, fasta_clusters, mcp)