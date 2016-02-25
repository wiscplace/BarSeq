#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
program: SplitBarSeq.py

purpose: Split Fastq (DNA Bar-Seq data) by sampleID.

input  : Bar-Seq Fastq file, Sequence_ID file

    Bar-Seq Fastq file - standard fastq file 
    
    Sequence_ID file   - plain text w/ 4 columns, provided by sequencing facility,
    only columns 1 & 4 matter. This identifies the experiment which is associated with 
    a particular user
                         
            Sample ID       Seq Index ID    Sequence Tag    User
            1248 T0 A       A1      AATAGGCGCT    kevin
            1248 T0 B       B1      AGCGTATGTC    kevin
            1248 T0 C       C1      AGTATGCACC    kevin

output  : Separate Fastq files for each sample, based on Sequence Tag
            
author  : mplace
date    : Feb 25th 2016

"""
#from Bio import SeqIO
from collections import defaultdict
import argparse
import itertools
import os
import sys
        
def matchSeqId( read, seqID ):
    """
    Match a sequence read to an experiment with a seqID tag.
    seqID tag must match start of string and be an EXACT match.
    seqID tag is assumed to be the 1st 10 bp of read.
    """
    seqTag = read[:10]       # assumes 1st 10 bp are the sampleID (experiment ID)
    if seqTag in seqID:
        return seqID[seqTag]
        
    return None
                        
def processFastq( fastq, seqID ):
    """
    Open and read fastq file read by read, calling function matchSeqId
    If match found write read (all 4 lines) to sample fastq file.
    """
    names = set()
    for keys, val in seqID.items():            # unique user names list
        names.add(val)

    for i in names:                            # create empty user fastq files for output
        newFastq = i + ".fastq"
        os.mknod(newFastq)        
        
    with open( fastq, "rU" ) as data:
        for header, line2, sign, qual in itertools.zip_longest(*[data]*4):
            read    = line2.rstrip()
            sample = matchSeqId(read, seqID)  #exact matching
            if sample:
                with open(sample + '.fastq','a') as out:
                    out.writelines([header, read, "\n", sign, qual])            
                            
def getSeqIdTag( seqIDFile ):
    """
    Get SeqID Up_Tag Decode information, put in dictionary
    Key = SampleID, value = sampleID_tag_biotechID
    example:
        key = CGTATCATGG value = 1248T6B_CGTATCATGG_D3
    """
    seqUpTag = dict()
        
    with open( seqIDFile , 'r') as f:
        for line in f:
            if line.startswith('Sample'):
                continue
            line = line.rstrip()
            items = line.split()
            sampleID = "".join(items[0:3])
            sampleID = sampleID + "_" + items[4] + "_" + items[3]
            if items[4] in seqUpTag:
                print("Duplicate seqUpTag %s" %(items[4]))
            else:
                seqUpTag[items[4]] = sampleID
               
    return seqUpTag
                    
def main():
    cmdparser = argparse.ArgumentParser(description="Split DNA Bar-Seq fastq file by sampleID into separate fastq files.",
                                        usage='%(prog)s -f fastq -s Seq_ID Up_Tag file', prog='SplitBarSeq.py'  )                                  
    cmdparser.add_argument('-f', '--Fastq' , action='store'     , dest='FASTQ' , help='BarSeq fastq file'         , metavar='')  
    cmdparser.add_argument('-s', '--Seqid' , action='store'     , dest='SEQID' , help='Seq_ID file (experiment ID)', metavar='')
    cmdparser.add_argument('-i', '--info'  , action='store_true', dest='INFO'  , help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['INFO']:
        print("\n  SplitBarSeq.py -f fastq file -s Seq_ID Up_Tag file")
        print("\n  Purpose: Create a new fastq file for each sampleID.")
        print("\n  Input  : BarSeq fastq file, Seq_ID file")
        print()
        print("        Seq_ID Up_Tag file:")
        print("        Sample ID       Seq Index ID    Sequence Tag    User")
        print("        1248 T0 A       A1      AATAGGCGCT    kevin ")
        print("        1248 T0 B       B1      AGCGTATGTC    maria")
        print()
        print("\n  Usage  : SplitBarSeq.py -f data.fastq -s SeqID.txt")        
        print("  ")       
        print("\tTo see Python Docs for this program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport SplitBarSeq")
        print("\thelp(SplitBarSeq)")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    # evaluate command line arguments & check if they exist
    if cmdResults['FASTQ']:
        fastq = cmdResults['FASTQ']
        if not os.path.exists(fastq):
            print("\n\t*** The %s fastq file could not be found! ***\n\n" %(fastq) )
            cmdparser.print_help()
            sys.exit(1)
        
    if cmdResults['SEQID']:
        seqIDFile = cmdResults['SEQID']
        if not os.path.exists(seqIDFile):
            print("\n\t*** The %s seqID Decode file could not be found! ***\n\n" %(seqIDFile) )
            cmdparser.print_help()
            sys.exit(1)            
    
    # Start processing the data
    seqID = getSeqIdTag(seqIDFile)   # looks like 1248T1B ['TATCGGAAGC', 'maria']  
    processFastq( fastq, seqID )     # open fastq and separate reads by user


if __name__ == "__main__":
    main()           
