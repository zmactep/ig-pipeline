#!/bin/sh

EXPECTED_ARGS=5

if [ $# -ne $EXPECTED_ARGS ]
then
  echo "Usage: run.sh working_dir input.fasta V.fasta D.fasta J.fasta"
  echo "If you do not have any of V D J fasta files, just write path to non-existing file, utility will handle it"
  exit
fi

WORKDIR=$1
INPUT=$2
VGENE=$3
DGENE=$4
JGENE=$5

MPWD=`pwd`

GENE_ALIGNER="gene_aligner"
GENE_DETECTOR="gene_detector"

function run_aligner {
    cd ./$GENE_ALIGNER/bin
    ./$GENE_ALIGNER
    cd $MPWD
    echo "Done"
}

echo "Working directory: $WORKDIR"
echo "Input fasta file: $INPUT"
echo "V genes file: $VGENE"
echo "D genes file: $DGENE"
echo "J genes file: $JGENE"

echo "Creating config for V genes"
cat ./$GENE_ALIGNER/bin/config/config.template > ./$GENE_ALIGNER/bin/config/config.ini
echo "reference_file=$VGENE" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "input_file=$INPUT" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_align=${WORKDIR}/v_align" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_regions=${WORKDIR}/v_regions" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "Done"

echo "Run aligner for V gene"
run_aligner


echo "Creating config for D genes"
cat ./$GENE_ALIGNER/bin/config/config.template > ./$GENE_ALIGNER/bin/config/config.ini
echo "reference_file=$DGENE" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "input_file=$INPUT" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_align=${WORKDIR}/d_align" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_regions=${WORKDIR}/d_regions" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "Done"

echo "Run aligner for D gene"
run_aligner

echo "Creating config for J genes"
cat ./$GENE_ALIGNER/bin/config/config.template > ./$GENE_ALIGNER/bin/config/config.ini
echo "reference_file=$JGENE" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "input_file=$INPUT" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_align=${WORKDIR}/j_align" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "output_regions=${WORKDIR}/j_regions" >> ./$GENE_ALIGNER/bin/config/config.ini
echo "Done"

echo "Run aligner for J gene"
run_aligner

echo "Prepare to run detector for collected data"
echo "[data]" > ./$GENE_DETECTOR/bin/config/config.ini 
echo "output_file=${WORKDIR}/result.txt" >> ./$GENE_DETECTOR/bin/config/config.ini 
echo "v_data_file=${WORKDIR}/v_regions" >> ./$GENE_DETECTOR/bin/config/config.ini 
echo "d_data_file=${WORKDIR}/d_regions" >> ./$GENE_DETECTOR/bin/config/config.ini 
echo "j_data_file=${WORKDIR}/j_regions" >> ./$GENE_DETECTOR/bin/config/config.ini 
cd ./$GENE_DETECTOR/bin

echo "Start detector"
./$GENE_DETECTOR
echo "Done. You're results are saved to ${WORKDIR}/result.txt"


