#!/home/mplace/anaconda3/bin/python
"""
program: BarSeq.py

purpose: Produce 2 tables (exact and 1 mismatch) of gene counts from DNA Bar-Seq data.
         The Bar-Seq fastq file is expected to be for a single user, 
         if you need to split reads, use: SplitBarSeqByUser.py.

input  : Bar-Seq Fastq file, GENE UP_tag Decode file, Sequence_ID UP_tag Decode file

    Bar-Seq Fastq file - standard fastq file , assumed to be one users set of experiments
    GENE UP_tag Decode file - plain text filw w/ 2 columns identifying each yeast 
    gene with a DNA sequence
        
            Gene_ID UP_tag_Seq
            YAL001C ACTATATGTGAAGGCATGGC
            YAL002W ATACTGACAGCACGCATGGC
            YAL003W GACATATCAGCATACACGGC
    
    Sequence_ID UP_tag file   - plain text w/ 4 columns, provided by sequencing facility,
    only columns 1 & 4 matter. File should have only one user's name
                         
            Sample ID       Seq Index ID    Sequence Tag    User
            1248 T0 A       A1      AATAGGCGCT    kevin
            1248 T0 B       B1      AGCGTATGTC    kevin
            1248 T0 C       C1      AGTATGCACC    kevin

output  : 2 Text files w/ Table, with counts per gene for each experiment, 1 for exact matches, 1 allowing 1 mismatch

            GeneName    Experiment1    Experiment2
            YAL001C      read_count     read_count
            
author  : mplace
date    : Feb 5th 2016
"""
from Bio import SeqIO
from collections import defaultdict
import argparse
import os
import sys

class BarSeq( object ):
        
    def __init__( self, fastq, decode, seqID, cwd  ):
        """
        parameters:
        fastq  = BarSeq fastq file, expected to be for a single user
        decode = Gene UP_tag file 
        cwd    = current working directory
        seqID  = Experiment ID_tag file
        """
        self.fastq     = fastq  
        self.dir       = cwd
        self.seqid     = seqID
        self.decode    = decode
        self.geneList  = []                    # Gene names 

        self.geneUpTag = self.getGeneUpTag()   # key = gene names, values = list of sequence tags
        self.seqIdTag  = self.getSeqIdTag()    # key = sampleID, values = list [ sequnce, user ]
        self.exact_Results = defaultdict(dict) #  dict(key=experiment) value = Dict(key = Gene, value = Counts )
        self.fuzzy_Results = defaultdict(dict) # for results with 1 mismatch  
        
    def matchSeqId( self, read ):
        """
        Match a sequence read to an experiment with a seqID tag
        seqID tag must match start of string and be an EXACT match.
        """
        for key, value in self.seqIdTag.items():
            if read.startswith(value[0]):
                return key
        return
                    
    def matchGene( self, read ):
        """
        Match GENE UPTAG to sequence read,
        Tag must be an exact match
        Assumes the Gene part of the sequence begins at the 29th base of the read.
        If not you will need to change the base number in geneSeq = read[28:]
        to an appropriate value.
        """
        for key, value in self.geneUpTag.items():
            geneSeq = read[28:]
            if geneSeq.startswith(value[0]):             
                return key
        return
    
    def matchGeneFuzzy( self, read ):
        """
        Match GENE UPTAG allowing 1 mismatch using hamming distance (only considers substitution).
        Discard reads with 2 or more N's.
        sum(c1!=c2 for c1,c2 in zip(value[0], geneSeq ))
        """
        hammingDist = 0
        for key, value in self.geneUpTag.items():
            end = 28 + len(value[0])
            geneSeq = read[28:end]
            if geneSeq.count('N') < 2:           # discard reads with too many N's
                hammingDist = sum(c1!=c2 for c1,c2 in zip(value[0], geneSeq ))
                if hammingDist == 1:           
                    return key
            
        return None
            
    def processFastq( self ):
        """
        Open and read fastq file read by read, calling functions matchSEQID
        and matchGene functions.  If match found count and store result.
        """
        file = open(self.fastq, 'rU')
        for rec in SeqIO.parse(file, 'fastq'):
            read       = str(rec.seq)
            # EXACT MATCHING STARTS HERE
            SeqID =  self.matchSeqId(read)  # returns sample/experiment identifer
            if SeqID:
                geneName = self.matchGene(read)
                if geneName:
                    if SeqID in self.exact_Results:
                        if geneName in self.exact_Results[SeqID]:
                            self.exact_Results[SeqID][geneName] += 1
                        else:
                            self.exact_Results[SeqID][geneName] = 1
                    else:
                        self.exact_Results[SeqID][geneName] = 1              
            # CHECK FOR 1 mismatch if exact matching fails
                else:
                    fuzzy = self.matchGeneFuzzy(read)
                    if fuzzy:
                        geneName = fuzzy
                        if SeqID in self.fuzzy_Results:
                            if geneName in self.fuzzy_Results[SeqID]:
                                self.fuzzy_Results[SeqID][geneName] += 1
                            else:
                                self.fuzzy_Results[SeqID][geneName] = 1
                        else:
                            self.fuzzy_Results[SeqID][geneName] = 1                   
        file.close()
    
    def mergeCounts( self ):
        """
        Sum exact match gene counts and fuzzy counts producing total counts.
        If fuzzy has no count, just use the exact count.
        """
        for seqid in self.exact_Results.keys():
            for gene, count in self.exact_Results[seqid].items():
                if seqid in self.fuzzy_Results:
                    if gene in self.fuzzy_Results[seqid] and gene in self.exact_Results[seqid]:
                        self.fuzzy_Results[seqid][gene] =  int(count) + int(self.fuzzy_Results[seqid][gene])
                    else:
                        if gene in self.exact_Results[seqid] and gene not in self.fuzzy_Results[seqid]:
                            self.fuzzy_Results[seqid][gene] = int(count)
        """for key in self.fuzzy_Results.keys():
            print("key: %s" %(key))
            for k in self.fuzzy_Results[key].keys():
                print("second key %s value %s" %(k, self.fuzzy_Results[key][k]))
        """        

    def writeTable( self, data, outName ):
        """
        Write gene counts to file as a tab delimited table.
        """ 
        sampleList = list(data.keys())                  # get list to maintain order when writing table      
        sampleList.sort()

        with open( outName, 'w') as out:
            out.write('gene')
            for name in sampleList:
                newName = name + "-" + self.seqIdTag[name][0]
                out.write("\t%s" %(newName))
            out.write("\n")
            
            for idx in range(0, len(self.geneList)):
                out.write("%s\t" %( self.geneList[idx]))
                for sample in sampleList:
                    if self.geneList[idx] not in self.exact_Results[sample]:
                        out.write('0\t')
                    else:
                        out.write('%s\t' %(self.exact_Results[sample][self.geneList[idx]]))
                out.write("\n")

    def getSeqIdTag( self ):
        """
        Get SeqID UP_tag Decode information, put in dictionary
        Key = SampleID, value = list[ sequence, user_name]
        """
        seqUpTag = defaultdict(list)
        
        with open( self.seqid , 'r') as f:
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
    
    def getGeneUpTag( self ):
        """
        Get GENE UP_tag Decode information, put in dictionary, 
        key = GENE,  value = list of sequence tags
        """
        geneUpTag = defaultdict(list)
        
        with open( self.decode, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('Y'):
                    items = line.split()
                    geneUpTag[items[0]].append(items[1])
                    if items[0] not in self.geneList:
                        self.geneList.append(items[0])
        
        return geneUpTag            
                    
def main():
    cmdparser = argparse.ArgumentParser(description="Produce a table of gene counts from DNA Bar-Seq data.",
                                        usage='%(prog)s -f fastq -d UP_tag_Decode file -s Seq_ID file', prog='BarSeq.py'  )                                  
    cmdparser.add_argument('-f', '--Fastq' , action='store'     , dest='FASTQ' , help='BarSeq fastq file'         , metavar='')
    cmdparser.add_argument('-d', '--Decode', action='store'     , dest='DECODE', help='GENE UP_tag_Decode file'        , metavar='')    
    cmdparser.add_argument('-s', '--Seqid' , action='store'     , dest='SEQID' , help='Seq_ID file (experiment ID)', metavar='')
    cmdparser.add_argument('-i', '--info'  , action='store_true', dest='INFO'  , help='Print a more detailed description of program.')
    cmdResults = vars(cmdparser.parse_args())
    
    # if no args print help
    if len(sys.argv) == 1:
        print("")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['INFO']:
        print("\n  BarSeq.py -f fastq file -d GENE UP_tag Decode file -s Seq_ID UP_tag Decode file")
        print("\n  Purpose: Produce a table of gene counts from DNA Bar-Seq data")
        print("\n  Input  : BarSeq fastq file, GENE UP_tag Decode file, Seq_ID UP_tag Decode file")
        print("\n    GENE UP_tag Decode file:")
        print("        Gene_ID UP_tag_Seq" )
        print("        YAL001C ACTATATGTGAAGGCATGGC")
        print()
        print("        Seq_ID UP_tag Decode file:")
        print("        Only one user name allowed.")
        print("        Sample ID       Seq Index ID    Sequence Tag    User")
        print("        1248 T0 A       A1      AATAGGCGCT    kevin ")
        print("        1248 T0 B       B1      AGCGTATGTC    kevin")
        print()
        print("\n  Usage  : BarSeq.py -f data.fastq -d GENE_UP_tag_Decode.txt -s SeqID.txt")        
        print("  ")       
        print("\tTo see Python Docs for this program:")
        print("\n\tOpen python console and enter")
        print("\timport sys")
        print("\tsys.path.append('/full/path/to/script')")
        print("\timport BarSeq")
        print("\thelp(BarSeq)")
        print("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    # evaluate command line arguments & check if they exist
    if cmdResults['FASTQ']:
        fastq = cmdResults['FASTQ']
        if not os.path.exists(fastq):
            print("\n\t*** The %s fastq file could not be found! ***\n\n" %(fastq) )
            cmdparser.print_help()
            sys.exit(1)
        
    if cmdResults['DECODE']:
        geneDecode = cmdResults['DECODE']
        if not os.path.exists(geneDecode):
            print("\n\t*** The %s gene Decode file could not be found! ***\n\n" %(geneDecode) )
            cmdparser.print_help()
            sys.exit(1)
        
    if cmdResults['SEQID']:
        seqID = cmdResults['SEQID']
        if not os.path.exists(seqID):
            print("\n\t*** The %s seqID Decode file could not be found! ***\n\n" %(seqID) )
            cmdparser.print_help()
            sys.exit(1)            
        
    cwd = os.getcwd()                                     # current working dir
    
    # Start processing the data
    data = BarSeq(fastq, geneDecode, seqID, cwd)
    
    #print(data.geneUpTag)
    #print(data.seqIdTag)    
    data.processFastq()
    data.writeTable(data.exact_Results,"Exact-Match.table")
    data.mergeCounts()
    data.writeTable(data.fuzzy_Results, "1-MisMatch.table")

if __name__ == "__main__":
    main()           