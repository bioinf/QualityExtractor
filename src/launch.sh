#!/bin/sh
VERSION="0.53"
java -jar $DIRNAME/QualityExtractor.jar $TSP_FILEPATH_BAM "/results/plugins/QualityExtractor/cftr.bed" $TSP_FILEPATH_GENOME_FASTA $RESULTS_DIR
