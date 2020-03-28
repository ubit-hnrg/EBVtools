#!/usr/bin/env nextflow

params.command = ""
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
sh = params.command


process pipeline {

    input:
    val samp from allSamples
    val ip from ipath
    val op from opath
    val reference from ref
    val maskf from mask
    file sh


    output:
    stdout result 

    script:
    """
    $sh --sample=$samp --inputpath=$ip/$samp --outpath=$op --reference=$reference --maskfile=$maskf
    """

}


result.subscribe { println it }