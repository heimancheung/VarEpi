#!/bin/bash

Genotype=$1
GWAS=$2
qv=$3
dir=$4

awk 'NR>1{print $1}' variable_ano_res.txt  |sed 's/-/\t/g' |intersectBed -a ${Genotype} -b - -wa -wb > genotype_overlap_peak.txt

Rscript ${dir}/SNP_associated_with_epi.R DESEQ2normalized_merge_matrix.txt genotype_overlap_peak.txt ${qv} ${GWAS} cluster.txt

