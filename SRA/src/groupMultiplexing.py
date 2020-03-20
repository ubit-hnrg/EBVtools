import pandas as pd
import os
import argparse
import glob

def get_args():
    parser = argparse.ArgumentParser(description='given a list of SRA ids, get its metadata info and download the associated RUN ids')
    parser.add_argument('--inputfile','-f',dest = 'inputfile',help='SRA list with metadata in tsv format')
    parser.add_argument('--outpath','-o',dest = 'outpath',help='root for the directory organization of sampes')
    parser.add_argument('--runidpath','-r',dest = 'runidpath',help='root for the directory of RUN accession')
    

    args = parser.parse_args()

    ifile = args.inputfile
    opath = args.outpath
    runidpath = args.runidpath

    return(ifile,opath,runidpath)

def makesyml(ifile,opath,runidpath):
    df = pd.read_csv(ifile,sep = '\t')
    samplenames = df.sample_accession
    mapa = {}
    for s in samplenames:
        mapa.update({s:df[df.sample_accession == s].run_accession.values})
        spath = opath+s
        checklen = 2* len(mapa[s]) # consider foward and reveerse     
        fastqs = []
        for runid in mapa[s]:
            fastqs = glob.glob(runidpath+runid+'/*.gz')

        if(len(fastqs) == checklen):
            if(not os.path.exists(spath)):
                print(1)
                os.makedirs(spath)

            for fqsource in fastqs:
                fq = os.path.basename(fqsource)    
                sfile = spath + '/' + fq
                src = os.path.abspath(fqsource)
                os.symlink(src,sfile)
        else:
            print 'sample %s do not contain all required fastq-files yet'%s

def main():
    ifile,opath,runidpath = get_args()
    makesyml(ifile,opath,runidpath)

if __name__ == '__main__':
    main()
