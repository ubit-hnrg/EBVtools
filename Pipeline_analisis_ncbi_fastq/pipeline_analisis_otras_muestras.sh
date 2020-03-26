set -e
type='EBV1'


## INSTALACIÓN DE PROGRAMAS ##

#Fastp:
#wget http://opengene.org/fastp/fastp
#chmod +x ./fastp

#PrintSeq (https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/prinseq.htm)
#wget -N http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz
#tar -zxvf prinseq-lite-0.20.4.tar.gz
#sudo cp -puv prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite && chmod +x /usr/local/bin/prinseq-lite
#sudo cp -puv prinseq-lite-0.20.4/prinseq-graphs.pl /usr/local/bin/prinseq-graphs && chmod +x /usr/local/bin/prinseq-graphs

#BWA (https://zoomadmin.com/HowToInstall/UbuntuPackage/bwa)
#sudo apt-get update -y
#sudo apt-get install -y bwa

#Bedtools (https://www.howtoinstall.co/es/ubuntu/xenial/bedtools)
#sudo apt-get update
#sudo apt-get install bedtools

#Bcftools mpileup (https://zoomadmin.com/HowToInstall/UbuntuPackage/bcftools)
#sudo apt-get update (lo hice en la descarga anterior)
#sudo apt-get install -y bcftools

#Samtools (https://zoomadmin.com/HowToInstall/UbuntuPackage/samtools)
#sudo apt-get update (lo hice en la descarga anterior)
#sudo apt-get install -y samtools

#vcf2bed
#   git clone https://github.com/bedops/bedops.git
# cd bedops
#make
#make install
# sudo cp bin/* /usr/local/bin


#Pandas
#sudo apt-get install python-pip
#sudo pip install numpy
#sudo pip install pandas

#Bio
#pip install biopython

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

referenceEBV='/data/EBV/analisis_NGS/Genoma_Referencia/'$type'/*.fa'  #localización de la secuencia de referencia indexada
reference='/data/EBV/analisis_NGS/Genoma_Referencia/'$type'/*NNN.fasta' #localizacion de la secuencia de referencia enmascarada con Ns en las zonas repetitivas
regionfile='/data/EBV/analisis_NGS/Coordenadas/'$type'/Coordenadas.'$type #coordenadas para eliminar del bam

samples=($(ls /data/EBV/ncbi/bySRAid/)) #localización de las muestras
declare -a samples=("ERS1791251")
#samples= : puedo declarar una muestra en particular para analizar

outpath='/data/EBV/analisis_NGS/Processed/'$type #localización del archivo final
outtrimmed='/data/EBV/analisis_NGS/Trimmed' #localización de las reads pre y post trimming y eliminación de duplicados

#sudo mkdir $outpath -p

for s in "${samples[@]}"
	do
	#echo $s
	
	##Generamos una nueva carpeta dentro del directorio Trimmed que contendrá las lecturas resultantes del preprocesamiento de las muestras
	mkdir $outtrimmed/$s  
	
	##Agrupamos los archivos descargados en dos grupos:
	zcat /data/EBV/ncbi/bySRAid/$s/*_1.fastq.gz |gzip > $outtrimmed/$s/$s.R1.fq.gz
	zcat /data/EBV/ncbi/bySRAid/$s/*_2.fastq.gz |gzip > $outtrimmed/$s/$s.R2.fq.gz

	##Trimmed y eliminacion de lecturas con baja calidad
	fastp --disable_adapter_trimming -i $outtrimmed/$s/$s.R1.fq.gz -I $outtrimmed/$s/$s.R2.fq.gz -o $outtrimmed/$s/$s.R1.trimmed.fq.gz \
        -O $outtrimmed/$s/$s.R2.trimmed.fq.gz -h \
        $outtrimmed/$s/$s.report_fastp.html -j $outtrimmed/$s/$s.report_fastp.json \
        --trim_front1=15 --trim_tail1=5 -e 28 --length_required 50;


	##Descomprimir las lecturas trimeadas : requisito necesario para entrar a prinseq
	gzip -c -d $outtrimmed/$s/$s.R1.trimmed.fq.gz > $outtrimmed/$s/$s.R1.trimmed.fq
	gzip -c -d $outtrimmed/$s/$s.R2.trimmed.fq.gz > $outtrimmed/$s/$s.R2.trimmed.fq
	
       	##Remove duplicates
	prinseq-lite -fastq $outtrimmed/$s/$s.R1.trimmed.fq -fastq2 $outtrimmed/$s/$s.R2.trimmed.fq -out_good $outtrimmed/$s/$s.good.trimmed $outtrimmed/$s/$s.good.trimmed -derep 1
	gzip $outtrimmed/$s/$s.good.trimmed_1.fastq
	gzip $outtrimmed/$s/$s.good.trimmed_2.fastq
	rm $outtrimmed/$s/$s.*.trimmed.fq
	
	##Generamos una carpeta dentro del outpath con el nombre de la muestra
	mkdir $outpath/$s
 
	## Se comentarán las lineas superiores una vez analizadas las lecturas con  uno de los genomas de referencia para evitar procesarlas nuevamente

	#Map to reference
	bwa mem -K 100000000 -v 1 -t 4 $referenceEBV <(zcat $outtrimmed/$s/$s.good.trimmed_1.fastq.gz) <(zcat $outtrimmed/$s/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $outpath/$s/$s.bam
	

	#Remove duplicates and unmapped reads
	samtools view -b -F 1548 $outpath/$s/$s.bam > $outpath/$s/$s.mapped.bam  
	
	#Sort bam
	samtools sort $outpath/$s/$s.mapped.bam > $outpath/$s/$s.mapped.sorted.bam
	
	#Drop out repetitive regions (according to coords file /home/cata/). 
	samtools view -bL $regionfile $outpath/$s/$s.mapped.sorted.bam > $outpath/$s/$s.mapped.sorted.withoutrep.bam

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
	intersectBed -wa -v -a $outpath/$s/$s.calls.vcf.gz -b /data/EBV/analisis_NGS/Coordenadas/$type/repetitive.$type > body.vcf 
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
	modified_reference=$outpath/$s/$s'_modified_reference.fa'
	/home/cata/EBVtools/mask_reference.py -r $reference -c $outpath/$s/$s.COB-DEL.bed -o $modified_reference
	bcftools consensus -f $modified_reference $outpath/$s/$s.NonRep.calls.vcf.gz > $outpath/$s/$s.nonrep.consensus.fa
	
	
done
