set -e
type='EBV2'
referenceEBV='/home/cata/EBV/Genoma_referencia/'$type'/*.fa'
reference='/home/cata/EBV/Genoma_referencia/'$type'/*NNN.fasta'
regionfile='/home/cata/EBV/coordenadas/'$type'/Coordenadas.'$type

samples=($(ls /mnt/datos/EBV/ncbi/fastq/))
declare -a samples=("ERX588930" "ERX588931" "ERX588933" "ERX588934" "ERX588935" "ERX588936" "ERX588937" "ERX588928" "ERX588926" "ERX218634" "ERX588892" "ERX588923" "ERX588924" "ERX588925" "ERX588929" "ERX218638" "ERX588893" "ERX218626" "ERX218636")
#samples= 'SRR9599785'
outpath='/mnt/datos/EBV/ncbi/Processed/'$type

outtrimmed='/mnt/datos/EBV/ncbi/Processed/'$type
#sudo mkdir $outpath -p



#s1='/home/ariel/corona/data/PRJNA601736/wuhan1_1.fq'
#s1='/home/ariel/corona/analysis/PRJNA601630/HKU-SZ-005b/HKU-SZ-005b_all.fastq'
#s2='/home/ariel/corona/data/PRJNA601736/wuhan1_2.fq'
for s in "${samples[@]}"
	do
	#echo $s
	mkdir $outpath/$s -p
	#zcat /home/cata/EBV/ANALISIS-OTRAS-SEQ/FASTQ/$s/*_1.fastq.gz |gzip > $outtrimmed/$s/$s.R1.fq.gz
	#zcat /home/cata/EBV/ANALISIS-OTRAS-SEQ/FASTQ/$s/*_2.fastq.gz |gzip > $outtrimmed/$s/$s.R2.fq.gz
	
	fastp --disable_adapter_trimming -i /mnt/datos/EBV/ncbi/fastq/$s/$s'_1.fastq.gz' -I /mnt/datos/EBV/ncbi/fastq/$s/$s'_2.fastq.gz' -o $outtrimmed/$s/$s.R1.trimmed.fq -O $outtrimmed/$s/$s.R2.trimmed.fq --trim_front1=15 --trim_tail1=5 -e 28 --length_required 50;

	##Remove duplicates

	prinseq-lite -fastq $outtrimmed/$s/$s.R1.trimmed.fq -fastq2 $outtrimmed/$s/$s.R2.trimmed.fq -out_good $outtrimmed/$s/$s.good.trimmed $outtrimmed/$s/$s.good.trimmed -derep 1
	sudo gzip $outtrimmed/$s/$s.good.trimmed_1.fastq
	sudo gzip $outtrimmed/$s/$s.good.trimmed_2.fastq
	
	## map to reference
	bwa mem -K 100000000 -v 1 -t 4 $referenceEBV <(zcat $outtrimmed/$s/$s.good.trimmed_1.fastq.gz) <(zcat $outtrimmed/$s/$s.good.trimmed_2.fastq.gz) | samtools view -b - > $outpath/$s/$s.bam
	

	## remove duplicates and unmapped reads
	samtools view -b -F 1548 $outpath/$s/$s.bam > $outpath/$s/$s.mapped.bam  
	
	## sort bam
	samtools sort $outpath/$s/$s.mapped.bam > $outpath/$s/$s.mapped.sorted.bam
	
	# drop out repetitive regions (according to coords file /home/cata/). 
	samtools view -bL $regionfile $outpath/$s/$s.mapped.sorted.bam > $outpath/$s/$s.mapped.sorted.withoutrep.bam

	## create bam index
	samtools index $outpath/$s/$s.mapped.sorted.withoutrep.bam
	
	# compute statistics for whole bam
	samtools stats $outpath/$s/$s.bam > $outpath/$s/$s.stats

	# compute statistics for final bam (without repetitive regions nither unmapped reads)
	samtools stats $outpath/$s/$s.mapped.sorted.withoutrep.bam > $outpath/$s/$s.mapped.sorted.withoutrep.stats
	
	## get vcf
	bcftools mpileup -f $referenceEBV $outpath/$s/$s.mapped.sorted.withoutrep.bam |bcftools call -mv --ploidy 1 -o $outpath/$s/$s.calls.vcf
	
	bgzip -c $outpath/$s/$s.calls.vcf > $outpath/$s/$s.calls.vcf.gz
	tabix $outpath/$s/$s.calls.vcf.gz

	#Cat VCF
	zgrep '^#' $outpath/$s/$s.calls.vcf.gz > header
	intersectBed -wa -v -a $outpath/$s/$s.calls.vcf.gz -b /home/cata/EBV/coordenadas/$type/repetitive.$type > body.vcf 
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
	/home/cata/EBV/src/mask_reference.py -r $reference -c $outpath/$s/$s.COB-DEL.bed -o $modified_reference
	bcftools consensus -f $modified_reference $outpath/$s/$s.NonRep.calls.vcf.gz > $outpath/$s/$s.nonrep.consensus.fa
	
	
done
