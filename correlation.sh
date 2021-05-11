#!/bin/bash

dir=$1
method=$2

multiBamSummary bins --bamfiles ${dir}/*.bam --minMappingQuality 30  -out sample.npz
plotCorrelation -in sample.npz --corMethod ${method} --skipZeros --plotTitle "spearman Correlation of sample" --whatToPlot heatmap -o  heatmap.pdf   --plotNumbers --outFileCorMatrix panicle.tab --plotFileFormat pdf