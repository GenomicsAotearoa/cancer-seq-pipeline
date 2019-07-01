
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
    
	
    # This needs to come from a config file rather than be hardcoded
    #dataDir="/nesi/project/uoa02461/data/raw/gastric/"
    

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
    except NameError:
        print("intermediateDirectory not defined correctly in pipeConfig.py")



    currentPatient = ""
    
    
    workflow = "run { [  control * [trim + align.using(type:'control') + removeDuplicates + removeSuplementary ], samples *  [trim + align.using(type:'test') + removeDuplicates + removeSuplementary  ]] + samples * [pileUp] + samples *[annotation] + count + reportGeneration } \n\n "
    #//QC run
    #//run{ [control * [qc] + control * [trim + align.using(type:'control')] + control * [alignmentMetrics.using(type:'control'), collectMetrics.using(type:'control')]  , samples * [qc] + samples *  [trim + align.using(type:'test')] + samples * [ alignmentMetrics.using(type:'control') ,collectMetrics.using(type:'test')]] }"

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
            pipelineFile.write("baseDir=\""+ baseDirectory +"\"\n")
 	    pipelineFile.write("singularityBuilds=\"" + baseDirectory + "/bin/\"\n")
            pipelineFile.write("bwaIndex =\"" + bwaIndex + "\"\n")
            pipelineFile.write("bwaIndexDir =\"" + bwaIndexDir + "\"\n")
            pipelineFile.write("dataDirectory =\"" + dataDirectory + "\"\n")


            pipelineFile.write("intermediateDirectory=\""+ intermediateDirectory + "\"\n\n")

            pipelineFile.write("control = [ \n\t"  + patient + sampleNumber + " :[ \""+ dataDirectory + patient  + sampleNumber + "_1.fastq.gz\",\"" + dataDirectory + patient  + sampleNumber + "_2.fastq.gz\"]\n]\n")
            pipelineFile.write("samples = [\n\t")

            currentPatient = patient

        else:
            pipelineFile.write(",\n")


        baseNumber = str(row["compareTo"])
        pipelineFile.write(patient + baseNumber + " : [\"" + dataDirectory + patient  + baseNumber + "_1.fastq.gz\", \"" + dataDirectory + patient +  baseNumber + "_2.fastq.gz\"]")  

        
    pipelineFile.write("\n]\n")
    pipelineFile.write(workflow)

    
    #for stage in postfix:
    #    pipelineFile.write( "\"%_*\" * [" + stage + "] + ")
    #pipelineFile.write("\tbranches * [trim]  + branches * [align] \n") 
    #pipelineFile.write(" branches * [pileUp]\n")
    #pipelineFile.write("}\n" )


if __name__ == "__main__":
   main(sys.argv[1:])


