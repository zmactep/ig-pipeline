#!/bin/sh

EXPECTED_ARGS=5

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: svm_find_regions.sh input_train.fasta input_kabat input_predict.fasta svm_window_size avg_window_size"
  exit
fi


echo "Generating datafiles..."
cd ./svm_data_generator/bin
./svm_data_generator train ../../$1 ../../$2 $4 > traindata.svm
./svm_data_generator predict ../../$3 $4 > predictdata.svm
echo "Done."

cd ../../
echo "training SVM..."

cd ./common_lib/third_party/libsvm-string-3.17/
./svm-train -h 0 -t 5 ../../../svm_data_generator/bin/traindata.svm 
./svm-predict ../../../svm_data_generator/bin/predictdata.svm traindata.svm.model output.txt
cd ../../..
echo "Done"

paste ./svm_data_generator/bin/read_names.txt ./common_lib/third_party/libsvm-string-3.17/output.txt > data.txt
python ./parse_svm_output.py data.txt $5 $4

rm data.txt ./common_lib/third_party/libsvm-string-3.17/traindata.svm.model ./common_lib/third_party/libsvm-string-3.17/output.txt ./svm_data_generator/bin/read_names.txt
echo "Done. Result is in results.txt. Debug output is in debug_prediction.txt and debug_prediction_avg.txt"
