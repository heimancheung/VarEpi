#!/bin/bash

module load FastQC/0.11.5-Java-1.8.0_92
module load Trimmomatic/0.32-Java-1.8.0_92
module load BWA/0.7.17
module load SAMtools/1.9
module load picard/2.1.1-Java-1.8.0_92

index=$1
genome_length=$2
control_bam=$3

date=$(date +%Y%m%d)

mkdir -p ${date}_statDir


#########1.预处理
#paired end
	

	sample=$5


    mkdir -p fastqc_output

    java -jar $EBROOTTRIMMOMATIC/trimmomatic-0.32.jar PE -threads 8 ${sample}_R1.fq.gz ${sample}_R2.fq.gz ${sample}_R1_trim.fq.gz ${sample}_R1_trim_unpaired.fq.gz  ${sample}_R2_trim.fq.gz  ${sample}_R2_trim_unpaired.fq.gz  ILLUMINACLIP:/public/home/xwzhang/software/trimmomatic-0.38/adapters/TruSeq3-PE.fa:2:30:10:8:True  SLIDINGWINDOW:4:15 MINLEN:50 HEADCROP:10 LEADING:5 TRAILING:5 2>> ./${sample}_trimStat.log && echo "###${sample}" &&  fastqc ${sample}*trim.f*q.gz -o fastqc_output 


#########.比对
##################mapping&samtoolsflagstat#######

###ChIP-seq

        fq_R1_file=./${sample}_R1_trim.fq.gz
        fq_R2_file=./${sample}_R2_trim.fq.gz
        mkdir -p ./${date}_sam
        sam_file=./${date}_sam/${sample}.sam
       log=./${date}_sam/${sample}.log

      bwa mem -t 8  ${index}  ${fq_R1_file}  ${fq_R2_file}  > ${sam_file} 2> ${log} && samtools flagstat ./${date}_sam/${sample}.sam >> ./${date}_statDir/flagstat.log && echo "###${sample}" >> ./${date}_statDir/flagstat.log



#########3.sam转bam（过滤Q30和DUPLICATE)

######samtools .sam->.bam
##########check tools samtools



    mkdir -p ${date}_bam
    sam_file=./${date}_sam/${sample}.sam
    sort_file=./${date}_sam/${sample}_sort.bam
    temp_file=./${date}_sam/${sample}_temp.bam
    unique_bam=./${date}_bam/${sample}.bam

    samtools view -q 30 -h -b -S   ${sam_file} | samtools sort -T ${temp_file} -o ${sort_file} && java -jar ${EBROOTPICARD}/picard.jar  MarkDuplicates REMOVE_DUPLICATES=true  I=${sort_file} O=${unique_bam} M=${sample}.metrics > ./${date}_sam/${sample}_picard.log 2>&1


module unload FastQC/0.11.5-Java-1.8.0_92
module unload Trimmomatic/0.32-Java-1.8.0_92
module unload BWA/0.7.17
module load deepTools/2.5.3
module load MACS2/2.1.1
module load BEDTools/2.27
module load R/3.5.2



		cd ${date}_bam
        mkdir -p ../${date}_statDir/qual
		
		
#########4.1 visualized file		
       samtools index  ${sample}.bam && bamCoverage -b ${sample}.bam --normalizeUsingRPKM -o ${sample}.bw
#########4.2 RSC&NSC	   
       Rscript /public/home/xwzhang/phantompeakqualtools/run_spp.R -c=./${sample}.bam -savp -out=../${date}_statDir/qual/${sample}.qual > ../${date}_statDir/qual/${sample}.Rout

	   
	   
#########5.call peak	   
        treat_bam=${sample}.bam
        macs2_name=${sample}
        out_dir=./${sample}_callpeaks


		if [ $4 -eq 2 ]
		then
			###########################broad peak
			mkdir -p ${sample}_callpeaks_broad
			macs2 callpeak -t ${treat_bam} -n ${macs2_name}  -c ${control_bam}  --outdir ${out_dir}_broad -f BAMPE -g ${genome_length} -B --broad    > ${out_dir}_broad/${sample}.log  2>&1 && bamToBed -i ${sample}.bam > ${sample}.bed && intersectBed -a ./${sample}.bed -b ${out_dir}_broad/${sample}_peaks.broadPeak -c -f 0.20 > ./${sample}_broad.intersectBed && rm ./${sample}.bed
		elif [ $4 -eq 1 ]
		then
			mkdir -p ${sample}_callpeaks
			##########################narrow peak   
			macs2 callpeak -t ${treat_bam} -n ${macs2_name}  -c ${control_bam}  --outdir ${out_dir} -f BAMPE -B -q 0.05 -g ${genome_length} > ./${sample}_callpeaks/${sample}.log  2>&1 && bamToBed -i ${sample}.bam > ${sample}.bed && intersectBed -a ./${sample}.bed -b ${sample}_callpeaks/${sample}_peaks.narrowPeak -c -f 0.20 > ./${sample}.intersectBed && rm ./${sample}.bed

		elif [ $4 -eq 3 ]
		then
			#########################DNase & FAIRE & ATAC
			macs2 callpeak -t ${treat_bam} --nomodel --shift -100 --extsize 200  -n ${macs2_name} -B -q 0.05 -g ${genome_length} --outdir ${out_dir} > ./${sample}_callpeaks/${sample}.log 2>&1 && bamToBed -i ${sample}.bam > ${sample}.bed && intersectBed -a ./${sample}.bed -b ${sample}_callpeaks/${sample}_peaks.narrowPeak -c -f 0.20 > ./${sample}.intersectBed && rm ./${sample}.bed
		
		fi

		

		


# bash chip-seq.sh /public/home/xwzhang/reference/RIGW/RS2/MH63/lowerclass/BWA/MH63RS2_low.fasta 3.6e8 /public/home/xwzhang/GMR/20200816_bam/seedlings_WGC094053-HPM-066_combined.bam 1 AN20200817-RMCS-1085
