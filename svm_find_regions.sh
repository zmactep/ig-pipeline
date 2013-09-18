#!/bin/sh

EXPECTED_ARGS=7

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: svm_find_regions.sh input_train.fasta input_kabat input_predict.fasta svm_window_size sliding_window_size merge_threshold no_region_borders_in_train"
  exit
fi

echo "Start date: " `date` 
echo "Generating datafiles..."
cd ./svm_data_generator/bin
./svm_data_generator train ../../$1 ../../$2 $4 $7 > /dev/null 2> /dev/null
./svm_data_generator predict ../../$3 $4 > /dev/null 2> /dev/null
echo "Done."

cd ../../

echo "Prepare data..."
export CLASSPATH=./common_lib/third_party/weka-3.6.10/weka.jar
java -Xmx2048M weka.filters.unsupervised.attribute.NumericToNominal -i ./svm_data_generator/bin/train.libsvm -o train_nominal.arff 2> /dev/null
java -Xmx2048M weka.filters.unsupervised.attribute.NumericToNominal -i ./svm_data_generator/bin/predict.libsvm -o predict_nominal.arff 2> /dev/null
CLASS_NAMES=$(cat ./train_nominal.arff | grep class) 
#fix headers to make it look the same
if [ -e /bin/sed ]; then
	SED=sed
else
	SED=gsed
fi
$SED -r -i 's/([0-9]*) \{.*\}/\1 {65,66,67,71,78,84}/g' ./train_nominal.arff
$SED -r -i 's/([0-9]*) \{.*\}/\1 {65,66,67,71,78,84}/g' ./predict_nominal.arff
$SED -i "s/@attribute class {65,66,67,71,78,84}/$CLASS_NAMES/g" ./train_nominal.arff 
$SED -i "s/@attribute class {65,66,67,71,78,84}/$CLASS_NAMES/g" ./predict_nominal.arff 

echo "Train..."
java -Xmx2048M weka.classifiers.trees.RandomForest -I 10 -K 0 -S 1 -no-cv -p 0 -t ./train_nominal.arff -d ./model.model > /dev/null
echo "Current time: " `date`
echo "Predict.."
java -Xmx2048M weka.classifiers.trees.RandomForest -no-cv -p 0 -l ./model.model -T ./predict_nominal.arff > prediction.txt
tail -n +6 prediction.txt | awk '{print $3}' | cut -d: -f2 > prediction.filtered.txt 
echo "Done."
echo "Current time: " `date`

paste ./svm_data_generator/bin/read_names.txt ./prediction.filtered.txt > data.txt
python3.3 ./parse_svm_output.py --input_file data.txt --sliding_window_size $5 --merge_threshold $6

python3.3 ./data/compare_kabat.py --ref $2 --input results.txt
#rm data.txt ./svm_data_generator/bin/read_names.txt train_nominal.arff model.model predict_nominal.arff prediction.filtered.txt prediction.txt 
echo "Done. Result is in results.txt. Debug output is in debug_prediction.txt and debug_prediction_avg.txt. Comparison is in comparison.kabat"
echo "End date: " `date` 

