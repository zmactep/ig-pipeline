#!/bin/sh

mkdir ./tmp/
python ./split_fasta_to_learn_and_train_datasets.py ../../data/germline/human/VJK_combinations.fasta ./tmp/ 70
./train_model.sh ./tmp/train.fasta ../../data/nomenclature/human/VJK_combinations.kabat 13 ./tmp/ ../
./predict.sh ../../data/nomenclature/human/VJK_combinations.kabat ./tmp/test.fasta 13 1 7 ./tmp/ ../ ./tmp/model.model
