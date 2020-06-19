#!/usr/bin/python

import pandas as pd
import sys, getopt
import os

from argparse import ArgumentParser
from pipeConfig import *



def main(argv):
    
    parser = ArgumentParser()
    parser.add_argument( dest="comparisons", help="write report to FILE", metavar="comparisonsFile(.csv)")

    # data locations
    #parser.add_argument( dest="dataDirectory", help="project base data location", metavar="Directory containing src and data directorird")
    #parser.add_argument( dest="config", help="location for intermediate files, must exist and be writeable", metavar="directory for intermediate files")

    parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")

    args = parser.parse_args()
    

    try:
        df = pd.read_csv (args.comparisons)
    except IOError:
        print("file " + filename + "doesn't exist or is malformed")
	
    try:
        dataDirectory
    except NameError:
	    print("dataDirectory not defined correctly in pipeConfig.py")		


    try:
        dataDirectory
    except NameError:
        print("dataDirectory not defined correctly in pipeConfig.py")

    try:
        intermediateDirectory
        if not os.path.exists(intermediateDirectory):
            os.makedirs(intermediateDirectory)
    except NameError:
        print("intermediateDirectory not defined correctly in pipeConfig.py")

    try:
        qcDirectory
        if not os.path.exists(qcDirectory):
            os.makedirs(qcDirectory)
    except NameError:
        print("qcDirectory not defined correctly in pipeConfig.py")

    try:
        tmpDirectory
    except NameError:
        print("tmp Directory not defined correctly in pipeConfig.py")




    currentPatient = ""
    
    ## original
    ##workflow = "run { [  control * [trim + align.using(type:'control') + markDuplicates.using(type:'control')  ], samples *  [trim + align.using(type:'test') + markDuplicates.using(type:'test')  ]] + samples * [pileUp] + samples *[annotation] + samples * [count] + samples *[reportGeneration] } \n\n "

    #workflow = "run { [  control * [trim + align.using(type:'control') + markDuplicates  ], samples *  [trim + align.using(type:'test') + markDuplicates  ]] + samples * [pileUp] + samples *[annotation] + samples * [count] + samples *[reportGeneration] } \n\n "

    # ADTex run
    workflow = "run { [  control * [trim + align.using(type:'control') + markDuplicates.using(type:'control')  ], samples *  [trim + align.using(type:'test') + markDuplicates.using(type:'test')  ]] + samples * [pileUp] + samples *[annotation] + samples *[adtex]  } \n\n "
    #workflow = "run { [  control * [trim + align.using(type:'control') + markDuplicates  ], samples *  [trim + align.using(type:'test') + markDuplicates  ]] + samples * [pileUp] + samples *[annotation] + samples *[adtex]  } \n\n "

    #//QC run
   # workflow = "run{ [control * [qc] + control * [trim +  align.using(type:'control') + alignmentMetrics.using(type:'control')]  , samples * [qc] + samples *  [trim +  align.using(type:'test') +  alignmentMetrics.using(type:'test') ]] }\n \n\n run{ [control * [qc] + control * [trim + qc], samples * [qc] + samples *  [trim + qc]]}"

    for index, row in df.iterrows():
        patient = row["patient"]
        sampleNumber = str(row["sampleId"])
        print(patient)
        print(sampleNumber)
        if patient != currentPatient:
            ## close the previous file
            if currentPatient!="":
                pipelineFile.write("\n]\n")
                pipelineFile.write(workflow)
		
            # open the new one
            pipelineFile  = open("minimalPipe-"+ patient, "w")
            pipelineFile.write("load \"pipeline\" \n\n")

            pipelineFile.write("baseDirectory=\""+ baseDirectory +"\"\n")
            pipelineFile.write("singularityBuilds=\"" + singularityBuilds + "\"\n")
            pipelineFile.write("bwaIndex =\"" + bwaIndex + "\"\n")
            pipelineFile.write("referenceDirectory =\"" + referenceDirectory + "\"\n")
            pipelineFile.write("hg19RefDirectory =\"" + hg19RefDirectory + "\"\n")
            pipelineFile.write("refData =\"" + refData + "\"\n")
            pipelineFile.write("resultsDirectory=\"" + resultsDirectory + "\"\n")
			
            pipelineFile.write("dataDirectory =\"" + dataDirectory + "\"\n")
            pipelineFile.write("qcDirectory =\"" + qcDirectory + "\"\n")
            pipelineFile.write("intermediateDirectory=\""+ intermediateDirectory +  "\"\n")
            pipelineFile.write("expandedIntermediate =\"" + expandedIntermediate + "\"\n")
            pipelineFile.write("tmpDirectory =\"" + tmpDirectory + "\"\n\n")




            pipelineFile.write("control = [ \n\t\""  + patient + sampleNumber + "\" :[ \""+ dataDirectory + patient  + sampleNumber + "_1.fastq.gz\",\"" + dataDirectory + patient  + sampleNumber + "_2.fastq.gz\"]\n]\n")
            pipelineFile.write("samples = [\n\t")

            currentPatient = patient

        else:
            pipelineFile.write(",\n")


        baseNumber = str(row["compareTo"])
        pipelineFile.write("\"" + patient + baseNumber +  "\" :[\"" + dataDirectory + patient  + baseNumber + "_1.fastq.gz\", \"" + dataDirectory + patient +  baseNumber + "_2.fastq.gz\"]")  

        
    pipelineFile.write("\n]\n")
    pipelineFile.write(workflow)

    
    #for stage in postfix:
    #    pipelineFile.write( "\"%_*\" * [" + stage + "] + ")
    #pipelineFile.write("\tbranches * [trim]  + branches * [align] \n") 
    #pipelineFile.write(" branches * [pileUp]\n")
    #pipelineFile.write("}\n" )


if __name__ == "__main__":
   main(sys.argv[1:])


