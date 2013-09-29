#!/bin/sh

EXPECTED_ARGS=6

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: predict.sh input_kabat input_predict.fasta svm_window_size sliding_window_size merge_threshold working_dir"
  exit
fi

echo "Start date: " `date` 
echo "Generating predict data in libsvm format..."
if [ ! -f ./svm_data_generator/bin/svm_data_generator ] 
then 
  echo "svm_data_generator not found. Abort."
  exit
fi

./svm_data_generator/bin/svm_data_generator predict $2 $3 ${6}> /dev/null 2> /dev/null

if [ ! -f ${6}predict.libsvm ] 
then 
  echo "Error in svm_data_generator: no output found. Abort."
  exit
fi

echo "Done. Applying NumericToNominal conversion..."
export CLASSPATH=./common_lib/third_party/weka-3.6.10/weka.jar
java -Xmx4096M weka.filters.unsupervised.attribute.NumericToNominal -i ${6}predict.libsvm -o ${6}predict_nominal.arff 2> /dev/null
if [ ! -f ${6}predict_nominal.arff ] 
then 
  echo "Error in NumericToNominal conversion: no output found. Abort."
  exit
fi

if [ ! -f ${6}model.model ] 
then 
  echo "Error: model.model not found. Abort."
  exit
fi

SKIPLINES=$(cat ${6}header.txt | wc -l)
cat ${6}header.txt > ${6}predict_nominal_tmp.arff
tail -n +$((SKIPLINES + 1)) ${6}predict_nominal.arff >> ${6}predict_nominal_tmp.arff

echo "Predict.."
java -Xmx4096M weka.classifiers.trees.RandomForest -no-cv -p 0 -l ${6}model.model -T ${6}predict_nominal_tmp.arff > ${6}prediction.txt
if [ ! -f ${6}prediction.txt ] 
then 
  echo "Error in prediction: ${6}prediction.txt not found. Abort."
  exit
fi

tail -n +6 ${6}prediction.txt | awk '{print $3}' | cut -d: -f2 > ${6}prediction.filtered.txt 
echo "Done."
echo "Current time: " `date`

paste ${6}read_names.txt ${6}prediction.filtered.txt > ${6}data.txt
python ./parse_svm_output.py --input_file ${6}data.txt --sliding_window_size $4 --merge_threshold $5 --output $6
python ./data/compare_kabat.py --ref $1 --input ${6}results.txt --output $6

echo "Done. Result is in results.txt. Debug output is in debug_prediction.txt and debug_prediction_avg.txt. Comparison is in comparison.kabat. All in $6"
echo "End date: " `date` 
