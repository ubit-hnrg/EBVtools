#!/usr/bin/env nextflow

params.in = "$HOME/input.txt"

samplefile = file(params.in)
allSamples  = samplefile.readLines()

process echoS {

    input:
    val str from allSamples

    output:
    stdout result 

    script:
    """
    echo $STR
    """

}


result.subscribe { println it }