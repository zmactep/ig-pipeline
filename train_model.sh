#!/bin/sh

EXPECTED_ARGS=4

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: train_model.sh input_train.fasta input_kabat svm_window_size working_dir"
  exit
fi

echo "Start date: " `date` 
echo "Generating train data in libsvm format..."
if [ ! -f ./svm_data_generator/bin/svm_data_generator ] 
then 
  echo "svm_data_generator not found. Abort."
  exit
fi

./svm_data_generator/bin/svm_data_generator train $1 $2 $3 1 $4> /dev/null 2> /dev/null
if [ ! -f ${4}/train.libsvm ] 
then 
  echo "Error in svm_data_generator: no output found. Abort."
  exit
fi

echo "Done. Applying NumericToNominal conversion..."
export CLASSPATH=./common_lib/third_party/weka-3.6.10/weka.jar
java -Xmx4096M weka.filters.unsupervised.attribute.NumericToNominal -i ${4}train.libsvm -o ${4}train_nominal.arff 2> /dev/null
if [ ! -f ${4}train_nominal.arff ] 
then 
  echo "Error in NumericToNominal conversion: no output found. Abort."
  exit
fi

#fix headers to make it look the same
python ./fix_weka_header.py ${4}train_nominal.arff $4
mv ${4}result.arff ${4}train_nominal.arff

echo "Train..."
java -Xmx4096M weka.classifiers.trees.RandomForest -I 10 -K 0 -S 1 -no-cv -p 0 -t ${4}train_nominal.arff -d ${4}model.model > /dev/null
echo "Done. All files in $4"
echo "End date: " `date` 
