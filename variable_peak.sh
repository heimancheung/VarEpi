#!/bin/bash

ds=$1
Fvalue=$2
CF=$3
qv=$4
group=$5
dir=$6


Rscript ${dir}/variable_peak_ancestry_variable_peak.R ${ds} ${Fvalue} ${CF} ${qv} ${group}
