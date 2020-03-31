#!/bin/bash
#please use in this way:
# ./Map2Refpipeline.sh --sample='samplename' --inputpath=/path/to/samplefolder \
#	--outpath=/path/to/output_folder --reference=/path/to/reference_file.fa \
#	--maskfile=/path/to/masking_file.bed --FilterBinaryCode=1548

# please, notice that
# 1) --maskfile and --FilterBinaryCode are optional paramemeters
# 2) --maskfile is expected as a 3 column, tabulated bed file BUT ONE-BASED . 
set -e

for i in "$@"
do
# default values
maskfile='None'
FilterBinaryCode="1548"

case $i in
    -S=*|--sample=*)
    sample="${i#*=}"

    ;;
    -i=*|--inputpath=*)
    inputpath="${i#*=}"

    ;;
    -o=*|--outpath=*)
    outpath="${i#*=}"

    ;;
    -r=*|--reference=*)
    reference="${i#*=}"

    ;;
    -m=*|--maskfile=*)  # this expect one-based bedfile for masking reference
    maskfile="${i#*=}"

    ;;
    -F=*|--FilterBinaryCode=*)  
    FilterBinaryCode="${i#*=}"
    ;;
    *)

            # unknown option
    ;;
esac
	#default paramenters
done

####################################################
# create output path if do not exist
# create variables for writing outputs
outp=$outpath/$sample
mkdir -p $outp
outtrimmed=$outp/'trimmed'/$sample

s=$sample #alias
referenceEBV=$reference
mask=$maskfile

inputs=$outpath/$sample/inputs.log
echo 'runing parameters'>$inputs
echo 'sample: '$sample >>$inputs
echo 'input directory: '$inputpath >>$inputs
echo 'output direrctory: '$outpath >>$inputs
echo 'reference file: '$reference >>$inputs
echo 'barcode for filtering: '$FilterBinaryCode >>$inputs
echo 'Masking File: '$mask >>$inputs


refdir=$outp/refGenomes
intervaldir=$outp/intervals
maskedReference=$refdir/maskedReference.fa
statsdir=$outp/stats
bamfolder=$outp/bams
vcffolder=$outp/vcfs

mkdir -p $refdir
mkdir -p $intervaldir
mkdir -p $statsdir
mkdir -p $bamfolder
mkdir -p $vcffolder
#####################################################


###############################################################################
## Get zero-based bedfile for masking in the right way. 

if [ "$mask" != "None" ];then
zeroMaskFile=$intervaldir/maskZB.bed

	cat $mask | while read chr start end 
	do
		echo -e $chr'\t'$(($start-1))'\t'$(($end-1)) 
	done  > $zeroMaskFile

	# mask original reference
	bedtools maskfasta -fi $referenceEBV -bed $zeroMaskFile -fo $maskedReference
fi


if [ "$mask" != "None" ];then
	interval=$intervaldir/keeping_intervalZB.bed
	#Get complement of masing file and name it $interval
	len=$(awk '/^>/{if (l!="") ;; l=0; next}{l+=length($0)}END{print l}' $referenceEBV)
	refid=$(head -n1 $mask|cut -f1)
	echo -e $refid'\t'$len > $intervaldir/referenceLength.tsv

	## amplio uno a izquierda y uno a derecha porque el complement de bedtools no permite sacar los extremos
	maskforcomplement=$intervaldir/maskforcomplement.tmp
		cat $zeroMaskFile | while read chr start end 
		do
			echo -e $chr'\t'$(($start-1))'\t'$(($end+1))   
		done  > $maskforcomplement

	bedtools complement -i $maskforcomplement -g $intervaldir/referenceLength.tsv > $interval
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
	prinseq-lite -fastq $outtrimmed/$s.R1.trimmed.fq -fastq2 $outtrimmed/$s.R2.trimmed.fq -out_good $outtrimmed/$s.good.trimmed -out_bad $outtrimmed/$s.bad.trimmed $outtrimmed/$s.good.trimmed -derep 1
	gzip $outtrimmed/$s.good.trimmed_1.fastq
	gzip $outtrimmed/$s.good.trimmed_2.fastq
	# Liberamos espacio

	rm $outtrimmed/$s.R1.fq.gz
	rm $outtrimmed/$s.R1.trimmed.fq
	rm $outtrimmed/$s.R1.trimmed.fq.gz
	rm $outtrimmed/$s.bad.trimmed_1.fastq

	rm $outtrimmed/$s.R2.trimmed.fq
	rm $outtrimmed/$s.R2.fq.gz
	rm $outtrimmed/$s.R2.trimmed.fq.gz
	rm $outtrimmed/$s.bad.trimmed_2.fastq

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
	<(zcat $outtrimmed/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $bamfolder/$s.bam
else
	bwa mem -K 100000000 -v 1 -t 4 $referenceEBV \
	<(zcat $outtrimmed/$s.good.trimmed_1.fastq.gz) \
	<(zcat $outtrimmed/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $bamfolder/$s.bam
fi

#Remove duplicates and unmapped reads
samtools view -b -F $FilterBinaryCode $bamfolder/$s.bam > $bamfolder/$s.mapped.bam  
	
#Sort bam
outbam=$bamfolder/$s.mapped.sorted.bam
samtools sort $bamfolder/$s.mapped.bam > $outbam
	
#Drop out repetitive regions (according to coords file /home/cata/). 
#if [ "$mask" != "None" ];
#then
#	outbam=$bamfolder/$s.mapped.sorted.withoutrep.bam
#	samtools view -bL $interval $bamfolder/$s.mapped.sorted.bam > $outbam
#else
#	
#fi

#Create bam index
samtools index $outbam

#Compute statistics for whole bam
samtools stats $bamfolder/$s.bam > $statsdir/$s.stats
#Compute statistics for final bam (without repetitive regions nither unmapped reads)
bn=$(basename $outbam)

samtools stats $outbam > $statsdir/$bn.stats

# liberamos espacio
rm $bamfolder/$s.bam
rm $bamfolder/$s.mapped.bam


#### end mapping stage (including bam statistics)
#################################################################################################



#################################################################################################
##################################### VARIANT CALLING STEP 
#Get vcf
if [ "$mask" != "None" ];
then
	bcftools mpileup -f $maskedReference $outbam |bcftools call -mv --ploidy 1 -o $vcffolder/$s.calls.vcf
else
	bcftools mpileup -f $referenceEBV $outbam |bcftools call -mv --ploidy 1 -o $vcffolder/$s.calls.vcf
fi

outvcf=$vcffolder/$s.calls.vcf
outvcfbgz=$vcffolder/$s.calls.vcf.gz
bgzip -c $vcffolder/$s.calls.vcf > $outvcfbgz
tabix $outvcfbgz


#Cat VCF
if [ "$mask" != "None" ];
	then
	ovcfbgz=$vcffolder/$s.NonRep.calls.vcf.gz
	outvcf=$vcffolder/$s.NonRep.calls.vcf
	
	zgrep '^#' $outvcfbgz > $vcffolder/header
	intersectBed -wa -a $outvcfbgz -b $interval > $vcffolder/body.vcf  ### 
	#intersectBed -wa -v -a $outvcfbgz -b /data/EBV/analisis_NGS/Coordenadas/$type/repetitive.$type > body.vcf 
	cat $vcffolder/header $vcffolder/body.vcf > $outvcf
	bgzip -c $vcffolder/$s.NonRep.calls.vcf > $ovcfbgz
	tabix $ovcfbgz
	rm $vcffolder/header
	rm $vcffolder/body.vcf
else
	ovcfbgz=$outvcfbgz
fi

#################### end of variant calling step ###############################
########################################################################################################


	#Coverange_0
	bedtools genomecov -ibam $outbam -bga | awk '$4==0' > $intervaldir/$s.Cob0.bed    ## This file is Zero-based

	#Deletion
	deletionfile=$intervaldir/$s.deletion.bed
	vcf2bed -n < $outvcf > $deletionfile  # Cata this file is One-based, so you can't mix it directly with bedfiles as Cob0.bed
	less $deletionfile | awk '{print $1, $2, $3}' > $intervaldir/$s.deletion.c.bed
	less $intervaldir/$s.deletion.c.bed |tr ' ' '\t' > $deletionfile
	rm $intervaldir/$s.deletion.c.bed
	# this file is  zero based
	
	
	############################################
	## get consensus sequence
	bedtools subtract -a $intervaldir/$s.Cob0.bed -b $deletionfile > $intervaldir/$s.COB-DEL.bed
	NonZeroCoverageReference=$refdir/$s'_NonZeroCoverageReference.fa'
	

if [ "$mask" != "None" ];
then
	# mask again but incorpore zero coverage regions
	bedtools maskfasta -fi $maskedReference -bed $intervaldir/$s.COB-DEL.bed -fo $NonZeroCoverageReference # Ok, the bedfile is zero-based
	outConsensus=$outp/$s.nonrep.nonzero.consensus.fa
else
	bedtools maskfasta -fi $referenceEBV -bed $intervaldir/$s.COB-DEL.bed -fo $NonZeroCoverageReference # Ok, the bedfile is zero-based
	outConsensus=$outp/$s.nonzero.consensus.fa
fi



bcftools consensus -f $NonZeroCoverageReference $ovcfbgz > $outConsensus

# Change name in consensus sequence, change the sammplename into the reference.fa
refidentifier=$(head -n1 $outConsensus|cut -f2 -d'>')
sed -i "s/$refidentifier/$s/g" $outConsensus

