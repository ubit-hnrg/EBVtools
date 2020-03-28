#!/usr/bin/python
import pandas as pd
from Bio import SeqIO
from Bio import Seq
from Bio.Seq import MutableSeq
import argparse


parser = argparse.ArgumentParser(prog= 'cobertura0V2.py', description='gets coverage')
parser.add_argument('-r', '--reference', help='sequence to mask in fasta format')
parser.add_argument('-c', '--coordenadas', help='coordenadas de cobertura 0 en cada muestra')
parser.add_argument('-o', '--outfile', help='referencia con N en regiones repetitivas y n en zonas de cobertura igual a 0')

args = parser.parse_args()
refN = args.reference
coordenadas = args.coordenadas
outfile = args.outfile

######  main  ######
ref = list(SeqIO.parse(refN,"fasta"))[0]

# load cob-del coordinates
coord= pd.read_csv(coordenadas,sep ='\t',header =None)
coord.columns =["seq","start","end", "other"]

newrefseq = ref.seq.tomutable()
for i in range(coord.shape[0]):
    l = coord.end[i] - coord.start[i]
    newrefseq[(coord.start[i]):(coord.end[i])] = ''.join(["N"]*l)
ref.seq = newrefseq

SeqIO.write(ref,outfile,"fasta")
