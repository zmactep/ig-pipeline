#!/bin/sh

PREDICTOR=BaggingREPTree
TRAIN_FASTA=/opt/ig-pipeline/data/test/germline/test-13/train-10/vh-train.fasta
TRAIN_KABAT=/opt/ig-pipeline/data/test/germline/test-13/train-10/vh-train.kabat
TEST_FASTA=/opt/ig-pipeline/data/test/germline/test-13/test-data/vh-test.fasta
TEST_KABAT=/opt/ig-pipeline/data/test/germline/test-13/test-data/vh-test.kabat
ML_WINDOW_SIZE=13
AVG_WINDOW_SIZE=3
MERGE_THRESHOLD=3
TOOLS_ROOT=/opt/ig-pipeline/ig-tools
OUTDIR=/tmp/tests

export PYTHONPATH=$PYTHONPATH:/opt/ig-pipeline/ig-tools/common_lib/python:/opt/ig-pipeline/ig-tools/ig-snooper:/opt/ig-pipeline/ig-tools/ig-snooper/ig_snooper_utils

rm -rf $OUTDIR
mkdir -p $OUTDIR

echo "Train..."
time python $TOOLS_ROOT/ig-snooper/train.py --predictor=$PREDICTOR --fasta=$TRAIN_FASTA --kabat=$TRAIN_KABAT --ml_window_size=$ML_WINDOW_SIZE --tools_root=$TOOLS_ROOT --outdir=$OUTDIR
echo "Done."

echo "Predict..."
time python $TOOLS_ROOT/ig-snooper/predict.py --predictor=$PREDICTOR --fasta=$TEST_FASTA --kabat=$TEST_KABAT --ml_window_size=$ML_WINDOW_SIZE --avg_window_size=$AVG_WINDOW_SIZE --merge_threshold=$MERGE_THRESHOLD --tools_root=$TOOLS_ROOT --outdir=$OUTDIR --model_path=$OUTDIR
echo "Done."

echo "Scores..."
python $TOOLS_ROOT/ig-snooper/ig_snooper_utils/compare_marking.py $OUTDIR/results.kabat $TEST_KABAT > $OUTDIR/compare_info.txt
python $TOOLS_ROOT/ig-snooper/ig_snooper_utils/diff_info.py $OUTDIR/results.kabat $TEST_KABAT > $OUTDIR/diff_info.txt
echo "Done."
