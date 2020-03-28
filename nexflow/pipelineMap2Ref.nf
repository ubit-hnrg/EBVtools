#!/usr/bin/env nextflow

params.sf = ""
params.ipath = ""
params.opath = ""
params.ref = ""
params.mask = ""

samplefile = file(params.sf)
allSamples  = samplefile.readLines()
ipath = file(params.ipath)
opath = file(params.opath)
ref = file(params.ref)
mask = file(params.mask)

process pipeline {

    input:
    val samp from allSamples
    val ip from ipath
    val op from opath
    val reference from ref
    val maskf from mask


    output:
    stdout result 

    script:
    """
    /home/hnrg/repos/EBVtools/pipelineNGS/Map2Refpipeline.sh --sample=$samp --inputpath=$ip --outpath=$op --reference=$reference --mask=$maskf
    """

}


result.subscribe { println it }