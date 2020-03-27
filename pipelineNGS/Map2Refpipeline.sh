set -e
#type='EBV1'



## ARMADO DEL DIRECTORIO ##
#Generamos una carpeta que contendrá las herramientas de análisis
# mkdir /home/cata/analisis_NGS

#Dentro de la carpeta "analisis_NGS" generamos las carpetas correspondientes:
# mkdir /home/cata/analisis_NGS/Coordenadas : contendrá las coordenadas de las regiones repetitivas que queremos eliminar para luego enmascarar con Ns
# mkdir /home/cata/analisis_NGS/Genoma_Referencia : contendrá el genoma de referencia en formato fasta, el genoma indexado para el mapeo de las reads con BWA y la referencia con Ns en la region repetitiva
# mkdir /home/cata/analisis_NGS/Trimmed : contendrá el analisis obtenido por fastp y PrintSeq
# mkdir /home/cata/analisis_NGS/Processed: contendrá el análisis posterior al preprocesamiento 

#Dentro de la carpeta "Coordenadas" generamos 2 carpetas correspondientes al tipo viral
#mkdir /home/cata/analisis_NGS/Coordenadas/'$type' : contendrá las coordenadas de la región repetitiva de EBV

#Dentro de la carpeta "Genoma_Referencia" generámos 2 carpetas corespondientes al tipo viral:
#mkdir /home/cata/analisis_NGS/Genoma_Referencia/'$type': contendrá el genoma de referencia del "type" y la referencia con Ns en la región repetitiva

#Indexar el genoma de referencia
#bwa index /home/cata/analisis_NGS/Genoma_Referencia´/'$type'/'$type'.fa 

#Dentro de la carpeta "Processed" generámos dos carpetas correspondientes al tipo viral
#mkdir /home/cata/analisis_NGS/Genoma_Referencia/'$type' :contendrá el bam y la consenso generada a partir de las reads pre-procesadas

sample=$1
inputdir=$2
outpath=$3
referenceEBV=$4
interval=$5
FilterBinaryCode='1548'

#referenceEBV='/data/EBV/analisis_NGS/Genoma_Referencia/'$type'/*.fa'  #localización de la secuencia de referencia indexada
#reference='/data/EBV/analisis_NGS/Genoma_Referencia/'$type'/*NNN.fasta' #localizacion de la secuencia de referencia enmascarada con Ns en las zonas repetitivas
#interval='/data/EBV/analisis_NGS/Coordenadas/'$type'/Coordenadas.'$type #coordenadas para eliminar del bam

## Mask reference
mkdir -p $outpath/$sample
maskedReference=$outpath/maskedReference.fa
bedtools maskfasta -fi $referenceEBV -bed $interval -fo $maskedReference

#samples=($(ls /data/EBV/ncbi/bySRAid/)) #localización de las muestras

#outpath='/data/EBV/analisis_NGS/Processed/'$type #localización del archivo final
#outtrimmed='/data/EBV/analisis_NGS/Trimmed' #localización de las reads pre y post trimming y eliminación de duplicados
outtrimmed=$outpath/'trimmed'


#cat $samplefile |while read s;
#	do
s=$sample
	echo 'processing sample id: '$s
	
	##Generamos una nueva carpeta dentro del directorio Trimmed que contendrá las lecturas resultantes del preprocesamiento de las muestras
	mkdir -p $outtrimmed  


################## PREPROCESSING #################
########### ONLY PERFORM THIS STEP IF PREPROCESSED FILE DO NOT EXIST ###########

if [ -f "$outtrimmed/$s.good.trimmed_1.fastq.gz" ]; then
	echo 'skiping preprocessing step'
else
	##Agrupamos los archivos descargados en dos grupos:
	zcat $inputdir/*_1.fastq.gz |gzip > $outtrimmed/$s.R1.fq.gz
	zcat $inputdir/*_2.fastq.gz |gzip > $outtrimmed/$s.R2.fq.gz

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
################     maping stage      #################

	#Map to reference
	bwa mem -K 100000000 -v 1 -t 4 $referenceEBV <(zcat $outtrimmed/$s.good.trimmed_1.fastq.gz) <(zcat $outtrimmed/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $outpath/$s/$s.bam
	
	#Remove duplicates and unmapped reads
	samtools view -b -F $FilterBinaryCode $outpath/$s/$s.bam > $outpath/$s/$s.mapped.bam  
	
	#Sort bam
	samtools sort $outpath/$s/$s.mapped.bam > $outpath/$s/$s.mapped.sorted.bam
	
	#Drop out repetitive regions (according to coords file /home/cata/). 
	samtools view -bL $interval $outpath/$s/$s.mapped.sorted.bam > $outpath/$s/$s.mapped.sorted.withoutrep.bam

	#Create bam index
	samtools index $outpath/$s/$s.mapped.sorted.withoutrep.bam
	
	#Compute statistics for whole bam
	samtools stats $outpath/$s/$s.bam > $outpath/$s/$s.stats

	#Compute statistics for final bam (without repetitive regions nither unmapped reads)
	samtools stats $outpath/$s/$s.mapped.sorted.withoutrep.bam > $outpath/$s/$s.mapped.sorted.withoutrep.stats
	
	#Get vcf
	bcftools mpileup -f $referenceEBV $outpath/$s/$s.mapped.sorted.withoutrep.bam |bcftools call -mv --ploidy 1 -o $outpath/$s/$s.calls.vcf
	
	bgzip -c $outpath/$s/$s.calls.vcf > $outpath/$s/$s.calls.vcf.gz
	tabix $outpath/$s/$s.calls.vcf.gz

	#Cat VCF
	zgrep '^#' $outpath/$s/$s.calls.vcf.gz > header
	intersectBed -wa -a $outpath/$s/$s.calls.vcf.gz -b $interval > body.vcf  ### 
	#intersectBed -wa -v -a $outpath/$s/$s.calls.vcf.gz -b /data/EBV/analisis_NGS/Coordenadas/$type/repetitive.$type > body.vcf 
	cat header body.vcf > $outpath/$s/$s.NonRep.calls.vcf
	bgzip -c $outpath/$s/$s.NonRep.calls.vcf > $outpath/$s/$s.NonRep.calls.vcf.gz
	tabix $outpath/$s/$s.NonRep.calls.vcf.gz
	rm header
	rm body.vcf

	#Coverange_0
	bedtools genomecov -ibam $outpath/$s/$s.mapped.sorted.bam -bga | awk '$4==0'| less > $outpath/$s/$s.Cob0.bed

	#Deletion
	vcf2bed -n < $outpath/$s/$s.NonRep.calls.vcf > $outpath/$s/$s.deletion.bed
	less $outpath/$s/$s.deletion.bed | awk '{print $1, $2, $3}' > $outpath/$s/$s.deletion.c.bed
	less $outpath/$s/$s.deletion.c.bed |tr ' ' '\t' > $outpath/$s/$s.deletion.bed
	rm $outpath/$s/$s.deletion.c.bed

	
	## get consensus sequence
	bedtools subtract -a $outpath/$s/$s.Cob0.bed -b $outpath/$s/$s.deletion.bed > $outpath/$s/$s.COB-DEL.bed
	NonZeroCoverageReference=$outpath/$s/$s'_NonZeroCoverageReference.fa'
	#/home/cata/EBVtools/mask_reference.py -r $maskedReference -c $outpath/$s/$s.COB-DEL.bed -o $modified_reference

	# mask again but incorpore zero coverage regions
	bedtools maskfasta -fi $maskedReference -bed $outpath/$s/$s.COB-DEL.bed -fo $NonZeroCoverageReference
	bcftools consensus -f $NonZeroCoverageReference $outpath/$s/$s.NonRep.calls.vcf.gz > $outpath/$s/$s.nonrep.consensus.fa
	
	
#done
