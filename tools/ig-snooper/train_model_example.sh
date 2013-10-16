#!/bin/sh

mkdir ./tmp2/
python ./split_fasta_to_learn_and_train_datasets.py ../../data/germline/human/VJK_combinations.fasta ./tmp2/ 70
./train_model.sh ./tmp2/train.fasta ../../data/nomenclature/human/VJK_combinations.kabat 13 ./tmp2/ ../
./predict.sh ../../data/nomenclature/human/VJK_combinations.kabat ./tmp2/test.fasta 13 1 7 ./tmp2/ ../ ./tmp2/model.model
