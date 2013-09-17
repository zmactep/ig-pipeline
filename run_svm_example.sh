#!/bin/sh

python ./data/split_fasta_to_learn_and_train_datasets.py ./data/germline/human/VJK_combinations.fasta 70
./svm_find_regions.sh ./train.fasta ./data/nomenclature/human/VJK_combinations.kabat test.fasta 11 20 10 0
