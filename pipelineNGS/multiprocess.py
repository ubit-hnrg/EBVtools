#Process level parallelism for shell commands
import glob
import subprocess as sp
import multiprocessing as mp
import argparse 

def myParser():
    parser = argparse.ArgumentParser(description='multiprocess version of Map2Refpipeline.sh')
    parser.add_argument('samplefile', dest='samplefile',actioon='store',type=str,
                    help='input path to samples names')
    parser.add_argument('--inputpath', dest='inputpath', action='store',type=str,
                    help='input path to sample folders')

    parser.add_argument('--outpath', dest='outpath', action='store',type=str,
                    help='output path to results')

    parser.add_argument('--reference', dest='reference', action='store',type=str,
                    help='reference fasta file')

    parser.add_argument('--maskfile', dest='maskfile', action='store',type=str,required=False,
                    help='mask file',default='None')

    parser.add_argument('--FilterBinaryCode', dest='FilterBinaryCode', action='store',type=str,
                    help='BinaryCode for filtering bam file',default='1548')

    args = parser.parse_args()
    return(args.samplefile,args.inputfile,args.outfile,args.reference,args.maskfile,args.FilterBinaryCode)


def work(sample,inputpath,outpath,reference,maskfile,FilterBinaryCode):
#    sp.call(['Map2Refpipeline.sh', '{}'.format(in_file), 'other', 'arguments'])
    sp.call(['/home/hnrg/repos/EBVtools/pipelineNGS/Map2Refpipeline.sh',
    '--sample=%s'%sample,
    '--inputpath=%s'%inputpath,
    '--outpath=%s'%outpath',
    '--reference=%s'%reference',
    '--maskfile=%s'%maskfile'
    '--FilterBinaryCode=%s'%FilterBinaryCode'
    ])
    return 0
 
if __name__ == '__main__':
    samplefile,inputpath,outpath,reference,maskfile,FilterBinaryCode = myParser()

    tasks = glob.glob(file_path)
    
    #Set up the parallel task pool to use all available processors
    count = mp.cpu_count()
    pool = mp.Pool(processes=count)
 
    #Run the jobs
    pool.map(work, tasks)