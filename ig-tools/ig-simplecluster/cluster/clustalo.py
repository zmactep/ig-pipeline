__author__ = 'mactep'

import os
import json
import shutil
import logging
import itertools
import difflib
from collections import defaultdict

import numpy as np
import matplotlib.pyplot as plt
from scipy.cluster.hierarchy import linkage, dendrogram, fcluster

from Bio import Phylo, SeqIO, AlignIO
from Bio.Align import AlignInfo, MultipleSeqAlignment

import common


def distance(record1, record2):
    return 1 - difflib.SequenceMatcher(None, str(record1.seq),
                                       str(record2.seq)).ratio()


# Fails in current version of Biopython
def save_utree(tree_path):
    logging.info("Saving unrooted tree.")
    prefix_name = tree_path[:tree_path.rfind('.')]
    unrooted_tree = Phylo.read(tree_path, 'newick')
    unrooted_tree.ladderize()
    Phylo.draw_graphviz(unrooted_tree)
    try:
        plt.savefig(common.UTREE_MASK % prefix_name)
    except:
        logging.error("Error while saving unrooted tree.")
        raise


def make_distance_matrix(align_file):
    logging.info("Copmputing distances matrix.")
    alignment_dict = SeqIO.to_dict(AlignIO.read(align_file, "fasta"))
    ids = dict(enumerate(alignment_dict.keys()))
    distance_matrix = np.zeros([len(ids)] * 2)
    for i, j in itertools.combinations(range(len(ids)), r=2):
        distance_matrix[i][j] = distance_matrix[j][i] = \
            distance(alignment_dict[ids[i]], alignment_dict[ids[j]])
    return distance_matrix, alignment_dict


def save_dendrogram(abs_out, distance_matrix, Y, cutoff):
    logging.info("Generating dendrogram.")
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
    try:
        fig.savefig(dendogram_path)
    except:
        logging.error("Error while saving dendrogram.")
        raise


def make_clusters(distance_matrix):
    logging.info("Computing clusters.")
    Y = linkage(distance_matrix, method="centroid")
    cutoff = 0.5 * max(Y[:, 2])
    clusters = fcluster(Y, cutoff, "distance")

    return clusters, Y, cutoff


def make_fasta_clusters(clusters, alignment_dict):
    logging.info("Extracing fasta clusters.")
    ids = dict(enumerate(alignment_dict.keys()))
    fasta_clusters = defaultdict(list)
    for i, cluster in enumerate(clusters):
        fasta_id = ids[i]
        fasta_clusters[cluster].append(alignment_dict[fasta_id])
    return fasta_clusters


def save_clusters(abs_out, fasta_clusters, mcp):
    logging.info("Saving clusters.")
    cluster_path = os.path.join(abs_out, common.CLUSTER_DIR)
    try:
        if os.path.isdir(cluster_path):
            shutil.rmtree(cluster_path)
        os.mkdir(cluster_path)
    except:
        logging.error("Cannot create clusters directory.")
        raise

    for cluster_id, cluster in fasta_clusters.items():
        try:
            SeqIO.write(cluster, os.path.join(cluster_path, common.CLUSTER_MASK % (mcp, cluster_id) + ".fasta"), 'fasta')
        except:
            logging.error("Error while saving cluster %s." % common.CLUSTER_MASK % (mcp, cluster_id))
            raise


def write_consensus(abs_out, fasta_clusters, mcp):
    logging.info("Writing consensus sequences.")
    cons = list()
    for cluster_id, cluster in fasta_clusters.items():
        consensus = AlignInfo.SummaryInfo(MultipleSeqAlignment(cluster)).dumb_consensus()
        cons.append(SeqIO.SeqRecord(consensus, common.CLUSTER_MASK % (mcp, cluster_id), name="", description=""))
    try:
        SeqIO.write(cons, os.path.join(abs_out, common.CONSENSUS_FILE), "fasta")
    except:
        logging.error("Error while saving consensus.")
        raise


def write_info(abs_out, fasta_clusters, mcp):
    logging.info("Writing JSON info")
    trash_size = len(SeqIO.to_dict(SeqIO.parse(os.path.join(abs_out, common.TRASH_FILE), "fasta")))
    groups = {common.CLUSTER_MASK % (mcp, cluster_id): [rec.id for rec in cluster]
              for cluster_id, cluster in fasta_clusters.items()}
    with open(os.path.join(abs_out, common.INFO_FILE), "wt") as fd:
        fd.write(json.dumps({"trash": trash_size,
                             "groups": groups}))


def run(src, out, minlen=None):
    abs_out = os.path.abspath(out)
    trash_path = os.path.join(abs_out, common.TRASH_FILE)
    copy_path = os.path.join(abs_out, os.path.basename(src))
    tree_path = os.path.join(abs_out, common.TREE_FILE)
    align_file = os.path.join(abs_out, common.ALIGNMENT_FILE)

    mcp = common.split_and_save(src, minlen, copy_path, trash_path)
    common.run_clustalo(copy_path, tree_path, align_file)
    #save_utree(tree_path)
    distance_matrix, alignment_dict = make_distance_matrix(align_file)
    clusters, Y, cutoff = make_clusters(distance_matrix)
    save_dendrogram(abs_out, distance_matrix, Y, cutoff)
    fasta_clusters = make_fasta_clusters(clusters, alignment_dict)
    save_clusters(abs_out, fasta_clusters, mcp)
    write_consensus(abs_out, fasta_clusters, mcp)
    write_info(abs_out, fasta_clusters, mcp)