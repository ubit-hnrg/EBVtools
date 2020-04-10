import pandas as pd
import os
import argparse
import glob

def get_args():
    parser = argparse.ArgumentParser(description='given a list of SRA ids, get its metadata info and download the associated RUN ids')
    parser.add_argument('--inputfile','-f',dest = 'inputfile',help='SRA list with metadata in tsv format')
    parser.add_argument('--outpath','-o',dest = 'outpath',help='root for the directory organization of sampes')
    parser.add_argument('--runidpath','-r',dest = 'runidpath',help='root for the directory of RUN accession')

    parser.add_argument('--paired', dest='paired', action='store_true')
    parser.add_argument('--no-paired', dest='paired', action='store_false')
    parser.set_defaults(paired=True)

    args = parser.parse_args()

    ifile = args.inputfile
    opath = args.outpath
    runidpath = args.runidpath
    paired = args.paired

    return(ifile,opath,runidpath,paired)

def makesyml(ifile,opath,runidpath,paired=True):
    df = pd.read_csv(ifile,sep = '\t')
    samplenames = df.sample_accession
    mapa = {}
    for s in samplenames:
        mapa.update({s:df[df.sample_accession == s].run_accession.values})
        spath = opath+s
        if(paired):
            checklen = 2* len(mapa[s]) # consider foward and reveerse
        else:
            checklen = len(mapa[s]) # consider only Foward
        fastqs = []
        for runid in mapa[s]:
            fastqs  = fastqs + glob.glob(runidpath+runid+'/*.gz')

        if(len(fastqs) == checklen):
            if(not os.path.exists(spath)):
                os.makedirs(spath)

            for fqsource in fastqs:
                fq = os.path.basename(fqsource)    
                sfile = spath + '/' + fq
                src = os.path.abspath(fqsource)
                if (not os.path.exists(sfile)):
                    print('symlink for sample %s runid %s'%(s,fq))
                    os.symlink(src,sfile)
        else:
            print 'sample %s do not contain all required fastq-files yet'%s

def main():
    ifile,opath,runidpath,paired = get_args()
    makesyml(ifile,opath,runidpath,paired)

if __name__ == '__main__':
    main()
