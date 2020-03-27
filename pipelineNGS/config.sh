## INSTALACIÓN DE PROGRAMAS ##

#Fastp:
wget http://opengene.org/fastp/fastp
chmod +x ./fastp

#PrintSeq (https://vcru.wisc.edu/simonlab/bioinformatics/programs/install/prinseq.htm)
wget -N http://downloads.sourceforge.net/project/prinseq/standalone/prinseq-lite-0.20.4.tar.gz
tar -zxvf prinseq-lite-0.20.4.tar.gz
sudo cp -puv prinseq-lite-0.20.4/prinseq-lite.pl /usr/local/bin/prinseq-lite && chmod +x /usr/local/bin/prinseq-lite
sudo cp -puv prinseq-lite-0.20.4/prinseq-graphs.pl /usr/local/bin/prinseq-graphs && chmod +x /usr/local/bin/prinseq-graphs

#BWA (https://zoomadmin.com/HowToInstall/UbuntuPackage/bwa)
sudo apt-get update -y
sudo apt-get install -y bwa

#Bedtools (https://www.howtoinstall.co/es/ubuntu/xenial/bedtools)
sudo apt-get update
sudo apt-get install bedtools

#Bcftools mpileup (https://zoomadmin.com/HowToInstall/UbuntuPackage/bcftools)
sudo apt-get update (lo hice en la descarga anterior)
sudo apt-get install -y bcftools

#Samtools (https://zoomadmin.com/HowToInstall/UbuntuPackage/samtools)
sudo apt-get update (lo hice en la descarga anterior)
sudo apt-get install -y samtools

#vcf2bed
git clone https://github.com/bedops/bedops.git
cd bedops
make
make install
sudo cp bin/* /usr/local/bin


#Pandas
sudo apt-get install python-pip
sudo pip install numpy
sudo pip install pandas

#Bio
pip install biopython
