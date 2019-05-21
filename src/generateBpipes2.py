import pandas as pd
import sys, getopt

from argparse import ArgumentParser




def main(argv):
    
    parser = ArgumentParser()
    parser.add_argument( dest="comparisons", help="write report to FILE", metavar="comparisonsFile(.csv)")
    parser.add_argument( dest="workflow", help="write report to FILE", metavar="workflow(.csv)")

    parser.add_argument("-q", "--quiet", action="store_false", dest="verbose", default=True, help="don't print status messages to stdout")

    args = parser.parse_args()

    try:
        df = pd.read_csv (args.comparisons)
    except IOError:
        print("file " + filename + "doesn't exist or is malformed")
        # doesn't exist
    try:
        workflow = pd.read_csv (args.workflow)
    except IOError:
        print("file " + filename + "doesn't exist or is malformed")
    

    currentPatient = ""
    

    workflow = "run { [  control * [trim + align.using(type:'control') + removeDuplicates + removeSuplementary ],   samples *  [trim + align.using(type:'test') + removeDuplicates + removeSuplementary ]] + samples * [  pileUp  ] }"

    for index, row in df.iterrows():
        patient = row["patient"]
        sampleNumber = str(row["sampleId"])

        if patient != currentPatient:
            ## close the previous file
            if currentPatient!="":
                pipelineFile.write("\n]\n")
                pipelineFile.write(workflow)

            # open the new one
            pipelineFile  = open("minimalPipe-"+ patient, "w")
            pipelineFile.write("load \"pipeline\" \n\n")
            pipelineFile.write("control = [ " +patient + sampleNumber + " : [\"" + patient  + sampleNumber + "_1.fastq.gz\",\"" + patient  + sampleNumber + "_2.fastq.gz\"]\n]\n")
            pipelineFile.write("samples = [\n")

            currentPatient = patient

        else:
            pipelineFile.write(",\n")


        baseNumber = str(row["compareTo"])
        pipelineFile.write(patient + baseNumber + " : [\"" + patient  + baseNumber + "_1.fastq.gz\", \"" + patient +  baseNumber + "_2.fastq.gz\"]")  

        
    pipelineFile.write("\n]\n")
    pipelineFile.write(workflow)

    
    #for stage in postfix:
    #    pipelineFile.write( "\"%_*\" * [" + stage + "] + ")
    #pipelineFile.write("\tbranches * [trim]  + branches * [align] \n") 
    #pipelineFile.write(" branches * [pileUp]\n")
    #pipelineFile.write("}\n" )


if __name__ == "__main__":
   main(sys.argv[1:])


