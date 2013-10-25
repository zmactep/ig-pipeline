import os

splitter_configuration = "splitter.conf.json"
splitter_least_len = 300
splitter_outfasta = "{}-splitted.fasta"
splitter_max_errors = 8

splitter_vhh_pattern = "ER[A-Z*]?[GSAF]"

svm_train_size = 30
svm_train_radius = 10

cluster_kmer_length = 33


def getConf(conf_variable):
    return os.path.join(os.path.abspath(os.path.dirname(__file__)),
                        conf_variable)
