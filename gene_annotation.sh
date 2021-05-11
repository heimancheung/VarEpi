#!/bin/bash

gff3=$1
cloneGene=$2


awk -v OFS="\t" -F "=|;|\t" '$3=="gene"{print $1,$4,$5,$10}' ${gff3} > all_gene.bed

awk 'NR>1{print $1"\t"$3}' non_variable_ano_res.txt|sed 's/-/\t/g'|intersectBed -a - -b all_gene.bed  -wa -wb |awk 'BEGIN{print "peak_ID\tF\tgene_ID"}{print $1"-"$2"-"$3"\t"$4"\t"$8}' > non_variable_ano_res.annotation
awk 'NR>1{print $1"\t"$3}' variable_ano_res.txt|sed 's/-/\t/g'|intersectBed -a - -b all_gene.bed  -wa -wb |awk 'BEGIN{print "peak_ID\tF\tgene_ID"}{print $1"-"$2"-"$3"\t"$4"\t"$8}' > variable_ano_res.annotation

awk 'NR>1{print $1"\t"$3}' variable_ano_res.txt|sed 's/-/\t/g'|intersectBed -a - -b all_gene.bed  -wa -wb |awk '{print $1"-"$2"-"$3"\t"$4"\t"$8}'|awk 'BEGIN{print "peak_ID\tF\tgene_ID\tclone_gene_Symbol"} NR==FNR{a[$1]=$2;next} $3 in a{print $0"\t"a[$3]} !($3 in a){print $0"\t-"}' ${cloneGene} - > variable_ano_res.clone.annotation
