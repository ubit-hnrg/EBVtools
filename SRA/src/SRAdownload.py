import pandas as pd
import os
import argparse
import timeit
import glob

def get_args():
    parser = argparse.ArgumentParser(description='given a list of SRA ids, get its metadata info and download the associated RUN ids')
    parser.add_argument('--samplefile','-f',dest = 'samplefile',help='SRA list')
    parser.add_argument('--outpath','-o',dest = 'outpath',help='SRA list')
    
    parser.add_argument('--download', dest='download', action='store_true')
    parser.add_argument('--no-download', dest='download', action='store_false')

    parser.add_argument('--prefetch', dest='prefetch', action='store_true')
    parser.add_argument('--no-prefetch', dest='prefetch', action='store_false')

    parser.set_defaults(download=True)
    parser.set_defaults(prefetch=False)  # prefetch mode is more sequre, its download the sra file and then convert it locally to fastq u other format

    args = parser.parse_args()

    ifile = args.samplefile
    opath = args.outpath
    download = args.download
    prefetch = args.prefetch

    return(ifile,opath,download,prefetch)

def params(ifile,opath, Filter = True):
    if not os.path.exists(opath):
        os.makedirs(opath   )
    base = os.path.basename(ifile)
    base = base.split('.')[0]

    idlist = pd.read_csv(ifile,header = None)
    sras = idlist[0].unique()
    samplequery = ' '.join(sras)
    query=opath+'/'+base
    os.system('pysradb srr-to-srp %s --expand --saveto %s.tsv'%(samplequery,query))
    #os.system('pysradb srr-to-srp %s --detailed --expand --saveto %s.tsv'%(samplequery,query))
    df = pd.read_table('%s.tsv'%query)
    df.to_excel('%s.xlsx'%query)
  
    if(Filter):
        run_acc = df[df['organism_taxid '] == 10376].run_accession 
    else:
        run_acc = df.run_accession 
        #run_acc = df[df['experiment_accession'] == 'ERX588931'].run_accession 
    return(run_acc)

def download(run_acc,opath,gzip = True,prefetch=False):
    for racc in run_acc:
        outp=opath + '/' + racc
        logfile = outp + '/' + racc + '.log'
        ofiles = outp +'/' + racc + '*fastq'
        if not os.path.exists(outp):
            os.makedirs(outp)

        if len(glob.glob(ofiles+'.gz'))>0:
            print('skiping sample %s'%racc)
            continue
        print('downloading %s'%racc)    
                    
        if(prefetch):
            start = timeit.default_timer()
            os.system('prefetch -O %s %s 2>%s'%(outp,racc,logfile))
            srafile=outp + '/'+racc+'.sra'
            end = timeit.default_timer()
            os.system('fastq-dump --split-files --gzip -O %s %s 2>>%s'%(outp,srafile,logfile))
            
            endfastq = timeit.default_timer()
            dt1= (end -start)/60
            dt2 = (endfastq - end)/60

            with open(logfile, 'a') as file:
                file.write('downloading SRA prefetch time (minutes): %.2f'%dt1)
                file.write('\n')
                file.write('dumping locally to fastq (minutes): %.2f'%dt2)
        else:
            start = timeit.default_timer()
            os.system('fasterq-dump -S -O %s %s 2>%s'%(outp,racc,logfile))
            end = timeit.default_timer()
            if(gzip):
                os.system('gzip %s'%ofiles)
            endgzip = timeit.default_timer()
            dt1= (end -start)/60
            dt2 = (endgzip - end)/60

            with open(logfile, 'a') as file:
                file.write('downloading time (minutes): %.2f'%dt1)
                file.write('\n')
                file.write('gzip time (minutes): %.2f'%dt2)

    return()

def main():
    ifile,opath,download,prefetch = get_args()
    run_acc = params(ifile,opath)
    if(download):
        download(run_acc,opath,prefetch=prefetch)

if __name__ == '__main__':
    main()
