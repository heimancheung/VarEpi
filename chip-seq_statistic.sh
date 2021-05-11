#######6.statistics
date=$1

##raw reads count
for i in *_trimStat.log ;do  echo -n ${i%_*}$'\t' >> ./${date}_statDir/raw_reads_count.txt; grep -o -e "Input Read Pairs: [0-9]* Both Surviving: [0-9]*"  $i >> ./${date}_statDir/raw_reads_count.txt;done && sed -i -e 's/Input Read Pairs: /\t/g' -e 's/Both Surviving: /\t/g' ./${date}_statDir/raw_reads_count.txt


cd ${date}_bam
##Frip
for i in *.intersectBed
do
NAME=${i%%.*}
reads=`cat ${i}|awk '{if($7>0) {print}}'|wc -l`
totalreads=`cat ${i}|wc -l`
FRiP=`bc <<< "scale=2;100*$reads/$totalreads"`
echo "${NAME}" >>  ../${date}_statDir/frip.txt
echo "reads: ${reads}" >>  ../${date}_statDir/frip.txt
echo "totalreads: ${totalreads}" >>  ../${date}_statDir/frip.txt
echo '==> FRiP value:' >>  ../${date}_statDir/frip.txt
echo "${FRiP}"'%' >>  ../${date}_statDir/frip.txt
echo " " >>  ../${date}_statDir/frip.txt
done		

##peak number
for i in *_callpeaks* 
do 
	if [ -e ${i}/*broadPeak ] ;then
		wc -l ${i}/*broadPeak >> ../${date}_statDir/peak_number.txt  
	elif [ -e ${i}/*narrowPeak ] ; then
		wc -l ${i}/*narrowPeak >> ../${date}_statDir/peak_number.txt
	else
		echo "None vaild peak file"
	fi
done

##Q30	
for i in *.bam ;do samtools view -c ${i} >> ../${date}_statDir/Q30.txt && echo ${i%.*} >> ../${date}_statDir/Q30.txt ;done 


cd ../${date}_sam
##duplicate reads
for i in *_sort.bam ;do samtools view -c ${i} >> ../${date}_statDir/duplicate.txt && echo ${i%.*} >> ../${date}_statDir/duplicate.txt ;done 


	
cd ../${date}_statDir/
###mapping rate
grep -E -o -e '[0-9]+ \+ [0-9]+ mapped \([0-9]*.[0-9]*% :' -e "###[[:print:]]+" flagstat.log |sed -r -e 's/mapped \(//' -e 's/ ://' -e 's/ \+ [0-9]+/\t/' -e 's/###//'|awk 'NR%2{rate=$0} !(NR%2){print $0"\t"rate}' >> mapping_rate.txt

##NSC&RSC
echo "sample	NSC	RSC" > ./qual/total.file
for i in ./qual/*.qual ;do cat ${i}|awk '{print $1 "\t" $9 "\t"$10 "\t"$11}' >> ./qual/total.file ;done


#######print result
echo "#raw reads count"
cat raw_reads_count.txt
echo ""

echo "#mapping rate"
cat mapping_rate.txt
echo ""

echo "#duplicate reads"
awk 'NR%2{count=$0} !(NR%2){print $0"\t"count}'  duplicate.txt
echo ""

echo "#Q30"
awk 'NR%2{count=$0} !(NR%2){print $0"\t"count}' Q30.txt
echo ""

echo "#peak number"
awk -F " |_" '{print $1"\t"$2}' peak_number.txt
echo ""

echo "#Frip"
awk 'NR%6==1{name=$0} NR%6==5{print name"\t"$0} !(NR%6==1)&&!(NR%6==5){next}' frip.txt
echo ""

echo "#NSC&RSC"
cat ./qual/total.file
echo ""



# bash chip-seq_statistic.sh 20201025
