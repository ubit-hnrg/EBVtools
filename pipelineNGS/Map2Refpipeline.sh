#!/bin/bash
set -e

sample=$1;
s=$sample #alias
inputpath=$2
outpath=$3
referenceEBV=$4
mask=$5
FilterBinaryCode='1548'

####################################################
# create output path if do not exist
# create variables for writing outputs
outp=$outpath/$sample
mkdir -p $outp
outtrimmed=$outp/'trimmed'/$s

refdir=$outp/refGenomes
maskedReference=$refdir/maskedReference.fa

mkdir -p $refdir
#####################################################


###############################################################################
## Get zero-based bedfile for masking in the right way. 

if [ "$mask" != "None" ];then
zeroMaskFile=$outp/maskZB.bed

	cat $mask | while read chr start end 
	do
		echo -e $chr'\t'$(($start-1))'\t'$(($end-1)) 
	done  > $zeroMaskFile

	# mask original reference
	bedtools maskfasta -fi $referenceEBV -bed $zeroMaskFile -fo $maskedReference
fi


if [ "$mask" != "None" ];then
	interval=$outp/keeping_intervalZB.bed
	#Get complement of masing file and name it $interval
	len=$(awk '/^>/{if (l!="") ;; l=0; next}{l+=length($0)}END{print l}' $referenceEBV)
	refid=$(head -n1 $mask|cut -f1)
	echo -e $refid'\t'$len > $outp/referenceLength.tsv

	## amplio uno a izquierda y uno a derecha porque el complement de bedtools no permite sacar los extremos
	maskforcomplement=$outp/maskforcomplement.tmp
		cat $zeroMaskFile | while read chr start end 
		do
			echo -e $chr'\t'$(($start-1))'\t'$(($end+1))   
		done  > $maskforcomplement

	bedtools complement -i $maskforcomplement -g $outp/referenceLength.tsv > $interval
fi

################################################################################


#cat $samplefile |while read s;
#	do

################## PREPROCESSING #################
########### ONLY PERFORM THIS STEP IF PREPROCESSED FILE DO NOT EXIST ###########
echo 'processing sample id: '$s
##Generamos una nueva carpeta dentro del directorio Trimmed que contendrá las lecturas resultantes del preprocesamiento de las muestras
mkdir -p $outtrimmed  

if [ -f "$outtrimmed/$s.good.trimmed_1.fastq.gz" ]; then
	echo 'skiping preprocessing step'
else
	##Agrupamos los archivos descargados en dos grupos:
	zcat $inputpath/*_1.fastq.gz |gzip > $outtrimmed/$s.R1.fq.gz
	zcat $inputpath/*_2.fastq.gz |gzip > $outtrimmed/$s.R2.fq.gz

	##Trimmed y eliminacion de lecturas con baja calidad
	fastp --disable_adapter_trimming -i $outtrimmed/$s.R1.fq.gz -I $outtrimmed/$s.R2.fq.gz -o $outtrimmed/$s.R1.trimmed.fq.gz \
        -O $outtrimmed/$s.R2.trimmed.fq.gz -h \
        $outtrimmed/$s.report_fastp.html -j $outtrimmed/$s.report_fastp.json \
        --trim_front1=15 --trim_tail1=5 -e 28 --length_required 50;


	##Descomprimir las lecturas trimeadas : requisito necesario para entrar a prinseq
	gzip -c -d $outtrimmed/$s.R1.trimmed.fq.gz > $outtrimmed/$s.R1.trimmed.fq
	gzip -c -d $outtrimmed/$s.R2.trimmed.fq.gz > $outtrimmed/$s.R2.trimmed.fq
	
    ##Remove duplicates
	prinseq-lite -fastq $outtrimmed/$s.R1.trimmed.fq -fastq2 $outtrimmed/$s.R2.trimmed.fq -out_good $outtrimmed/$s.good.trimmed $outtrimmed/$s.good.trimmed -derep 1
	gzip $outtrimmed/$s.good.trimmed_1.fastq
	gzip $outtrimmed/$s.good.trimmed_2.fastq
	rm $outtrimmed/$s.*.trimmed.fq
fi


####################################################################################
################     maping stage      #################
	#Map to reference
if [ "$mask" != "None" ];then
	bwa index $maskedReference
fi

if [ "$mask" != "None" ];
then
	bwa mem -K 100000000 -v 1 -t 4 $maskedReference \
	<(zcat $outtrimmed/$s.good.trimmed_1.fastq.gz) \
	<(zcat $outtrimmed/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $outp/$s.bam
else
	bwa mem -K 100000000 -v 1 -t 4 $referenceEBV \
	<(zcat $outtrimmed/$s.good.trimmed_1.fastq.gz) \
	<(zcat $outtrimmed/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $outp/$s.bam
fi

#Remove duplicates and unmapped reads
samtools view -b -F $FilterBinaryCode $outp/$s.bam > $outp/$s.mapped.bam  
	
#Sort bam
samtools sort $outp/$s.mapped.bam > $outp/$s.mapped.sorted.bam
	
#Drop out repetitive regions (according to coords file /home/cata/). 
outbam=$outp/$s.mapped.sorted.withoutrep.bam
if [ "$mask" != "None" ];
then
	samtools view -bL $interval $outp/$s.mapped.sorted.bam > $outbam
else
	outbam=$outp/$s.mapped.sorted.bam
fi

#Create bam index
samtools index $outbam

#Compute statistics for whole bam
samtools stats $outp/$s.bam > $outp/$s.stats
#Compute statistics for final bam (without repetitive regions nither unmapped reads)
samtools stats $outbam > $outbam.stats

#### end mapping stage (including bam statistics)
#################################################################################################



#################################################################################################
##################################### VARIANT CALLING STEP 
#Get vcf
if [ "$mask" != "None" ];
then
	bcftools mpileup -f $maskedReference $outbam |bcftools call -mv --ploidy 1 -o $outp/$s.calls.vcf
else
	bcftools mpileup -f $referenceEBV $outbam |bcftools call -mv --ploidy 1 -o $outp/$s.calls.vcf
fi

outvcf=$outp/$s.calls.vcf
outvcfbgz=$outp/$s.calls.vcf.gz
bgzip -c $outp/$s.calls.vcf > $outvcfbgz
tabix $outvcfbgz


#Cat VCF
if [ "$mask" != "None" ];
	then
	ovcfbgz=$outp/$s.NonRep.calls.vcf.gz
	outvcf=$outp/$s.NonRep.calls.vcf
	
	zgrep '^#' $outvcfbgz > header
	intersectBed -wa -a $outvcfbgz -b $interval > body.vcf  ### 
	#intersectBed -wa -v -a $outvcfbgz -b /data/EBV/analisis_NGS/Coordenadas/$type/repetitive.$type > body.vcf 
	cat header body.vcf > $outvcf
	bgzip -c $outp/$s.NonRep.calls.vcf > $ovcfbgz
	tabix $ovcfbgz
	rm header
	rm body.vcf
else
	ovcfbgz=$outvcfbgz
fi

#################### end of variant calling step ###############################
########################################################################################################


	#Coverange_0
	bedtools genomecov -ibam $outbam -bga | awk '$4==0' > $outp/$s.Cob0.bed    ## This file is Zero-based

	#Deletion
	deletionfile=$outp/$s.deletion.bed
	deletionfileZeroBased=$outp/$s.deletionZB.bed
	vcf2bed -n < $outvcf > $deletionfile  # Cata this file is One-based, so you can't mix it directly with bedfiles as Cob0.bed
	less $deletionfile | awk '{print $1, $2, $3}' > $outp/$s.deletion.c.bed
	less $outp/$s.deletion.c.bed |tr ' ' '\t' > $deletionfile
	rm $outp/$s.deletion.c.bed

	cat $deletionfile | while read chr start end; 
	do
		echo -e $chr'\t'$(($start-1))'\t'$(($end-1)) 
	done  > $deletionfileZeroBased

	
	
	############################################
	## get consensus sequence
	bedtools subtract -a $outp/$s.Cob0.bed -b $deletionfileZeroBased > $outp/$s.COB-DEL.bed
	NonZeroCoverageReference=$refdir/$s'_NonZeroCoverageReference.fa'
	#/home/cata/EBVtools/mask_reference.py -r $maskedReference -c $outp/$s.COB-DEL.bed -o $modified_reference

if [ "$mask" != "None" ];
then
	# mask again but incorpore zero coverage regions
	bedtools maskfasta -fi $maskedReference -bed $outp/$s.COB-DEL.bed -fo $NonZeroCoverageReference # Ok, the bedfile is zero-based
	outConsensus=$outp/$s.nonrep.nonzero.consensus.fa
else
	bedtools maskfasta -fi $referenceEBV -bed $outp/$s.COB-DEL.bed -fo $NonZeroCoverageReference # Ok, the bedfile is zero-based
	outConsensus=$outp/$s.nonzero.consensus.fa
fi
bcftools consensus -f $NonZeroCoverageReference $ovcfbgz > $outConsensus
