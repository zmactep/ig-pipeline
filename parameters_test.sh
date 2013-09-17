#!/bin/sh

python ./data/split_fasta_to_learn_and_train_datasets.py ./data/germline/human/VJK_combinations.fasta 70

for SVM_WINDOW_SIZE in 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31 
do
  for AVG_WINDOW_SIZE in 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31
  do
    for MERGE_THRESHOLD in 3 5 7 9 11 13 15 17 19 21 23 25 27 29 31
    do 
      for BORDERS in 0 1
      do
        ./svm_find_regions.sh ./train.fasta ./data/nomenclature/human/VJK_combinations.kabat test.fasta ${SVM_WINDOW_SIZE} ${AVG_WINDOW_SIZE} ${MERGE_THRESHOLD} ${BORDERS} > run_${SVM_WINDOW_SIZE}_${AVG_WINDOW_SIZE}_${MERGE_THRESHOLD}_${BORDERS}.log
        echo "SVM_WINDOW_SIZE: $SVM_WINDOW_SIZE AVG_WINDOW_SIZE: ${AVG_WINDOW_SIZE} MERGE_THRESHOLD: ${MERGE_THRESHOLD} BORDERS: ${BORDERS}"
        cat run_${SVM_WINDOW_SIZE}_${AVG_WINDOW_SIZE}_${MERGE_THRESHOLD}_${BORDERS}.log | grep 'Absolute score'
      done
    done
  done
done