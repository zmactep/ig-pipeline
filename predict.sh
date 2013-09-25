#!/bin/sh

EXPECTED_ARGS=5

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: predict.sh input_kabat input_predict.fasta svm_window_size sliding_window_size merge_threshold"
  exit
fi

echo "Start date: " `date` 
echo "Generating predict data in libsvm..."
if [ ! -f ./svm_data_generator/bin/svm_data_generator ] 
then 
  echo "svm_data_generator not found. Abort."
  exit
fi

./svm_data_generator/bin/svm_data_generator predict $2 $3 > /dev/null 2> /dev/null

if [ ! -f predict.libsvm ] 
then 
  echo "Error in svm_data_generator: no output found. Abort."
  exit
fi
mv ./predict.libsvm ./tmp/predict.libsvm
mv read_names.txt ./tmp/read_names.txt

echo "Done. Applying NumericToNominal conversion..."
export CLASSPATH=./common_lib/third_party/weka-3.6.10/weka.jar
java -Xmx4096M weka.filters.unsupervised.attribute.NumericToNominal -i ./tmp/predict.libsvm -o ./tmp/predict_nominal.arff 2> /dev/null
if [ ! -f ./tmp/predict_nominal.arff ] 
then 
  echo "Error in NumericToNominal conversion: no output found. Abort."
  exit
fi

if [ ! -f ./tmp/model.model ] 
then 
  echo "Error: model.model not found. Abort."
  exit
fi

SKIPLINES=$(cat ./tmp/header.txt | wc -l)
cat ./tmp/header.txt > ./tmp/predict_nominal_tmp.arff
tail -n +$((SKIPLINES + 1)) ./tmp/predict_nominal.arff >> ./tmp/predict_nominal_tmp.arff

echo "Predict.."
java -Xmx4096M weka.classifiers.trees.RandomForest -no-cv -p 0 -l ./tmp/model.model -T ./tmp/predict_nominal_tmp.arff > ./tmp/prediction.txt
if [ ! -f ./tmp/prediction.txt ] 
then 
  echo "Error in prediction: ./tmp/prediction.txt not found. Abort."
  exit
fi

tail -n +6 ./tmp/prediction.txt | awk '{print $3}' | cut -d: -f2 > ./tmp/prediction.filtered.txt 
echo "Done."
echo "Current time: " `date`

paste ./tmp/read_names.txt ./tmp/prediction.filtered.txt > ./tmp/data.txt
python ./parse_svm_output.py --input_file ./tmp/data.txt --sliding_window_size $4 --merge_threshold $5
python ./data/compare_kabat.py --ref $1 --input ./results.txt

mv ./results.txt ./tmp/results.txt
mv ./debug_prediction.txt ./tmp/debug_prediction.txt 
mv ./debug_prediction_avg.txt ./tmp/debug_prediction_avg.txt 
mv ./comparison.kabat ./tmp/comparison.kabat

echo "Done. Result is in results.txt. Debug output is in debug_prediction.txt and debug_prediction_avg.txt. Comparison is in comparison.kabat. All in ./tmp"
echo "End date: " `date` 
