#!/home/mplace/anaconda3/bin/python
"""
program: SplitBarSeqByUser.py

purpose: Split Fastq (DNA Bar-Seq data) by user.

input  : Bar-Seq Fastq file, Sequence_ID file

    Bar-Seq Fastq file - standard fastq file 
    
    Sequence_ID file   - plain text w/ 4 columns, provided by sequencing facility,
    only columns 1 & 4 matter. This identifies the experiment which is associated with 
    a particular user
                         
            Sample ID       Seq Index ID    Sequence Tag    User
            1248 T0 A       A1      AATAGGCGCT    kevin
            1248 T0 B       B1      AGCGTATGTC    kevin
            1248 T0 C       C1      AGTATGCACC    kevin

output  : Separate Fastq files for each user, based on Sequence Tag
            
author  : mplace
date    : Feb 15th 2016
Notes   : use of itertools speeds the intial script up by 4 times.
"""
#from Bio import SeqIO
from collections import defaultdict
import argparse
import itertools
import os
import sys
        
def matchSeqId( read, seqID ):
    """
    Match a sequence read to an experiment with a seqID tag
    seqID tag must match start of string and be an EXACT match.
    """
    for key, value in seqID.items():
        if read.startswith(value[0]):
            return value[1]
    return None
                        
def processFastq( fastq, seqID ):
    """
    Open and read fastq file read by read, calling function matchSeqId
    If match found count and store result.
    """
    names = set()
    for keys, val in seqID.items():            # unique user names list
        names.add(val[1])    

    for i in names:                            # create empty user fastq files for output
        newFastq = i + ".fastq"
        os.mknod(newFastq)        
        
    with open( fastq, "rU" ) as data:
        for line1, line2, line3, line4 in itertools.zip_longest(*[data]*4):
            header = line1.rstrip()
            read    = line2.rstrip()
            sign   = line3.rstrip()
            qual   = line4.rstrip()
            user = matchSeqId(read, seqID)  #exact matching
            if user:
                with open(user + '.fastq','a') as out:
                    out.writelines([header,"\n", read, "\n", sign, "\n", qual, "\n"])            
                            
def getSeqIdTag( seqIDFile ):
    """
    Get SeqID Up_Tag Decode information, put in dictionary
    Key = SampleID, value = list[ sequence, user_name]
    """
    seqUpTag = defaultdict(list)
        
    with open( seqIDFile , 'r') as f:
        for line in f:
            if line.startswith('Sample'):
                continue
            line = line.rstrip()
            items = line.split()
            sampleID = ["".join(items[0:3])]    
            if sampleID[0] in seqUpTag:
                print("Duplicate sampleID %s" %(sampleID[0]))
            else:
                seqUpTag[sampleID[0]].append(items[4])
                seqUpTag[sampleID[0]].append(items[5])
                
    return seqUpTag
                    
def main():
    cmdparser = argparse.ArgumentParser(description="Split Fastq (DNA Bar-Seq data) by user.",
                                        usage='%(prog)s -f fastq -s Seq_ID Up_Tag file', prog='SplitBarSeqByUser.py'  )                                  
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
        print("\n  SplitBarSeqByUser.py -f fastq file -s Seq_ID Up_Tag file")
        print("\n  Purpose: Create a new fastq file for each users data.")
        print("\n  Input  : BarSeq fastq file, Seq_ID file")
        print()
        print("        Seq_ID Up_Tag file:")
        print("        Sample ID       Seq Index ID    Sequence Tag    User")
        print("        1248 T0 A       A1      AATAGGCGCT    kevin ")
        print("        1248 T0 B       B1      AGCGTATGTC    maria")
        print()
        print("\n  Usage  : BarSeq.py -f data.fastq -s SeqID.txt")        
        print("  ")       
        print("\tTo see Python Docs for this program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport SplitBarSeqByUser")
        print("\thelp(SplitBarSeqByUser)")
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
        
    cwd = os.getcwd()                                     # current working dir
    
    # Start processing the data
    seqID = getSeqIdTag(seqIDFile)   # looks like 1248T1B ['TATCGGAAGC', 'maria']
    processFastq( fastq, seqID )     # open fastq and separate reads by user


if __name__ == "__main__":
    main()           