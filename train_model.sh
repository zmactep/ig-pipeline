#!/bin/sh

EXPECTED_ARGS=3

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: train_model.sh input_train.fasta input_kabat svm_window_size"
  exit
fi

mkdir tmp
echo "Start date: " `date` 
echo "Generating train data in libsvm..."
if [ ! -f ./svm_data_generator/bin/svm_data_generator ] 
then 
  echo "svm_data_generator not found. Abort."
  exit
fi

./svm_data_generator/bin/svm_data_generator train $1 $2 $3 1 > /dev/null 2> /dev/null
if [ ! -f train.libsvm ] 
then 
  echo "Error in svm_data_generator: no output found. Abort."
  exit
fi

mv ./train.libsvm ./tmp/train.libsvm
echo "Done. Applying NumericToNominal conversion..."
export CLASSPATH=./common_lib/third_party/weka-3.6.10/weka.jar
java -Xmx4096M weka.filters.unsupervised.attribute.NumericToNominal -i ./tmp/train.libsvm -o ./tmp/train_nominal.arff 2> /dev/null
if [ ! -f ./tmp/train_nominal.arff ] 
then 
  echo "Error in NumericToNominal conversion: no output found. Abort."
  exit
fi

#fix headers to make it look the same
python3.3 ./fix_weka_header.py ./tmp/train_nominal.arff
rm ./tmp/train_nominal.arff
mv ./result.arff ./tmp/train_nominal.arff
mv ./header.txt ./tmp/

echo "Train..."
java -Xmx4096M weka.classifiers.trees.RandomForest -I 10 -K 0 -S 1 -no-cv -p 0 -t ./tmp/train_nominal.arff -d ./tmp/model.model > /dev/null
echo "Done. All files in ./tmp"
echo "End date: " `date` 
