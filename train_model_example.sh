#!/bin/sh

python ./data/split_fasta_to_learn_and_train_datasets.py ./data/germline/human/VJK_combinations.fasta 70
./train_model.sh ./train.fasta ./data/nomenclature/human/VJK_combinations.kabat 13
./predict.sh ./data/nomenclature/human/VJK_combinations.kabat test.fasta 13 1 7
