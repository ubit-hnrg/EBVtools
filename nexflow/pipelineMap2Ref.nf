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
sh = file(params.cmd)


process pipeline {

    input:
    val samp from allSamples
    val ip from ipath
    val op from opath
    val reference from ref
    val maskf from mask
    val cmd from sh


    output:
    stdout result 

    script:
    """
    $cmd --sample=$samp --inputpath=$ip/$samp --outpath=$op --reference=$reference --maskfile=$maskf
    """

}


result.subscribe { println it }