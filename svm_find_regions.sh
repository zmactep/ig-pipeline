#!/bin/sh

EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: svm_find_regions.sh input_train.fasta input_kabat input_predict.fasta window_size"
  exit
fi


#echo "Generating datafiles..."
#cd ./svm_data_generator/bin
#./svm_data_generator train ../../$1 ../../$2 $4 > traindata.svm
#./svm_data_generator predict ../../$3 $4 > predictdata.svm
#echo "Done."

#cd ../../
#echo "traing SVM..."


#cd ./libsvm-3.17/
#./svm-train -h 0 ../svm_data_generator/bin/traindata.svm 
#./svm-predict ../svm_data_generator/bin/predictdata.svm traindata.svm.model output.txt
#cd ..
#echo "Done"

awk '{print $NF}' ./svm_data_generator/bin/predictdata.svm > read_names.txt
paste read_names.txt ./libsvm-3.17/output.txt > data.txt
python ./parse_svm_output.py > result.txt

#rm data.txt ./libsvm-3.17/traindata.svm.model ./libsvm-3.17/output.txt read_names.txt
#echo "Done. Result is in result.txt"