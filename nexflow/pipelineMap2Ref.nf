#!/usr/bin/env nextflow

params.command = ""
params.sf = ""
params.ipath = ""
params.opath = ""
params.ref = ""
params.mask = ""
params.filtercode = ""


samplefile = file(params.sf)
allSamples  = samplefile.readLines()
ipath = file(params.ipath)
opath = file(params.opath)
ref = file(params.ref)
mask = file(params.mask)
sh = file(params.cmd)
fc = params.filtercode


process pipeline {

    errorStrategy 'ignore'

    input:
    val samp from allSamples
    val ip from ipath
    val op from opath
    val reference from ref
    val maskf from mask
    val cmd from sh
    val fcode from fc

    output:
    stdout result 

    script:
    """
    $cmd --sample=$samp --inputpath=$ip/$samp --outpath=$op --reference=$reference --maskfile=$maskf --FilterBinaryCode=$fcode
    """

}


result.subscribe { println it }