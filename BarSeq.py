#!/home/GLBRCORG/mplace/anaconda3/bin/python
"""
program: BarSeq.py

purpose: Produce 2 tables (exact and 1 mismatch) of gene counts from DNA Bar-Seq data.
         The Bar-Seq fastq file is expected to be for a single user, 
         if you need to split the fastq file use: SplitBarSeqByUser.py.

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
            
method  : Step 1) Contruct a BarSeq object, containing the fastq file, Gene UP_tag and Sequence_ID UP_tag information, and 
                  containers for the count results.
          Step 2) Fastq file is precessed one read at a time.
                  First the sampleID (Experiment info) is identified using matchSeqID function : this matches the read start
                  with the Sequence_ID UP_tag.
                  If a match is found, then matchGene is called, this trims the first 28 bases off the read and sets the length
                  to the length of a Gene UP_tag and searches the geneUpTag dict using the sequence as a key.
                  If a match is found, the gene count is incremented.
                  If no match then a fuzzy search is used to try and match a gene with 1 base difference.
                  
author  : mplace
date    : Feb 5th 2016
version : 1.0
"""
from collections import defaultdict
import argparse
import itertools
import os
import sys
import GENE_UpTag_Dictionary as gn # a dict, key = exact sequence tags and all possible 1-mismatch sequence combinations, value = gene name

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
        self.exact_Results = dict()            #  dict(key=experiment) value = Dict(key = Gene, value = Counts )
        self.fuzzy_Results = dict()            # for results with 1 mismatch
        self.geneUpTag = self.getGeneUpTag()   # key = sequence tag values = gene name
        self.seqIdTag  = self.getSeqIdTag()    # key = sampleID, values = list [ sequence, user ]
        self.missingTag    = defaultdict(set) # key = sequence tag, values = gene name, tags found in exact match
         
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
        Match GENE UPTAG to sequence read, Tag must be an exact match
        Assumes the Gene part of the sequence begins at the 29th base of the read.
        This is true for ~99% of the Moby Up Tag sequences.
        Since the "gene tag" in decode file is of variable length, 
        we loop through the known ranges, 10-22.
        """
        for idx in range(10,23):
            end = 28 + idx
            geneSeq = read[28:end]
            if geneSeq in self.geneUpTag: 
                return self.geneUpTag[geneSeq]
        return 
    
    def matchGeneFuzzy( self, read ):
        """
        Match GENE UPTAG allowing 1 mismatch.
        """
        # walk down the sequence from the first start position
        for outidx in range(25,31):
            start = outidx 
            for idx in range(10,23):      # extract a sequence window
                end = start + idx
                geneSeq = read[start:end]
                if geneSeq in gn.tags: 
                    return gn.tags[geneSeq]
        return 
       
    def processFastq( self ):
        """
        Read in fastq file, calling functions matchGene and matchGeneFuzzy
        If match found count and store result.
        If no match found write sequence (all 4 lines) to Unmatched_reads.fastq
        """        
        with open( self.fastq, "rU" ) as data:
            for line1, line2, line3, line4 in itertools.zip_longest(*[data]*4):
                read    = line2.rstrip()                  # This is the sequence
                SeqID =  self.matchSeqId(read)            # returns sample/experiment identifer
                
                if SeqID in self.exact_Results:           
                    geneName = self.matchGene(read)       # Look for an Exact match  
                    if geneName:
                        if geneName in self.exact_Results[SeqID]:
                            self.exact_Results[SeqID][geneName] += 1
                        else:
                            self.exact_Results[SeqID][geneName] = 1
                    else:
                        fuzzyName = self.matchGeneFuzzy(read)
                        if fuzzyName:
                            if fuzzyName in self.fuzzy_Results[SeqID]:
                                self.fuzzy_Results[SeqID][fuzzyName] += 1
                            else:
                                self.fuzzy_Results[SeqID][fuzzyName] = 1
                        else:
                            with open('Unmatched_reads.fastq', 'a') as out:
                                out.write("SEQID: %s" %(SeqID) )
                                out.write(line1)
                                out.write(line2)
                                out.write(line3)
                                out.write(line4)
             
    def mergeCounts( self ):
        """
        Sum exact match gene counts and fuzzy counts producing total counts.
        If fuzzy has no count, just use the exact count.
        """
        
        for gene in self.geneList:
            for seqid in self.fuzzy_Results.keys():
                if gene in self.fuzzy_Results[seqid] and gene in self.exact_Results[seqid]:       # if gene in both
                    self.fuzzy_Results[seqid][gene] += self.exact_Results[seqid][gene]
                elif gene in self.fuzzy_Results[seqid] and gene not in self.exact_Results[seqid]: # if gene not in exact
                    self.exact_Results[seqid][gene] = 0
                elif gene not in self.fuzzy_Results[seqid] and gene in self.exact_Results[seqid]: # if gene not in fuzzy
                    self.fuzzy_Results[seqid][gene] = self.exact_Results[seqid][gene]
                else:                                                                             # if gene not in both
                    self.fuzzy_Results[seqid][gene] = 0
                    self.exact_Results[seqid][gene] = 0
                    self.missingTag[seqid].add(gene)

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
                    out.write('%s\t' %(data[sample][self.geneList[idx]]))
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
                self.exact_Results[sampleID[0]] = dict()
                self.fuzzy_Results[sampleID[0]] = dict()
                if sampleID[0] in seqUpTag:
                    print("Duplicate sampleID %s" %(sampleID[0]) )
                else:
                    seqUpTag[sampleID[0]].append(items[4])
                    seqUpTag[sampleID[0]].append(items[5])
                
        return seqUpTag
    
    def getGeneUpTag( self ):
        """
        Get GENE UP_tag Decode information, put in dictionary, 
        key = sequence uptag  value = gene
        """
        #geneUpTag = defaultdict(list)
        geneUpTag = dict()

        with open( self.decode, 'r') as f:
            for line in f:
                line = line.rstrip()
                if line.startswith('Y'):
                    items = line.split()
                    #geneUpTag[items[0]].append(items[1])
                    geneUpTag[items[1]] = items[0]
                    # Store a list of GENE Names
                    if items[0] not in self.geneList:
                        self.geneList.append(items[0])
        
        return geneUpTag            

    def writeMissingGenes( self ):
        """
        Write a table of sampleID's and the gene names which had no match.
        """       
        sampleID = list(self.seqIdTag.keys())
        sampleID.sort()
        
        with open('Genes_not_reported.txt', 'w') as out:
            for sample in sampleID:
                out.write(sample + ":")
                out.write("\t".join(self.missingTag[sample]))
                out.write('\n')
              
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
        print ("\n")
        cmdparser.print_help()
        sys.exit(1)
    
    if cmdResults['INFO']:
        print ("\n  BarSeq.py -f fastq file -d GENE UP_tag Decode file -s Seq_ID UP_tag Decode file")
        print ("\n  Purpose: Produce a table of gene counts from DNA Bar-Seq data")
        print ("\n  Input  : BarSeq fastq file, GENE UP_tag Decode file, Seq_ID UP_tag Decode file")
        print ("\n    GENE UP_tag Decode file:")
        print ("        Gene_ID UP_tag_Seq" )
        print ("        YAL001C ACTATATGTGAAGGCATGGC")
        print ("\n")
        print ("        Seq_ID UP_tag Decode file:")
        print ("        Only one user name allowed.")
        print ("        Sample ID       Seq Index ID    Sequence Tag    User")
        print ("        1248 T0 A       A1      AATAGGCGCT    kevin ")
        print ("        1248 T0 B       B1      AGCGTATGTC    kevin")
        print ("\n")
        print ("\n  Usage  : BarSeq.py -f data.fastq -d GENE_UP_tag_Decode.txt -s SeqID.txt")        
        print ("  ")     
        print ("\tTo see Python Docs for this program:")
        print ("\n\tOpen python console and enter")
        print ("\timport sys")
        print ("\tsys.path.append('/full/path/to/script')")
        print ("\timport BarSeq-fuzzy")
        print ("\thelp(BarSeq-fuzzy)")
        print ("\n\tSee Mike Place for any problems or suggestions.")
        sys.exit(1)
    
    # evaluate command line arguments & check if they exist
    if cmdResults['FASTQ']:
        fastq = cmdResults['FASTQ']
        if not os.path.exists(fastq):
            print ("\n\t*** The %s fastq file could not be found! ***\n\n" %(fastq) )
            cmdparser.print_help()
            sys.exit(1)
        
    if cmdResults['DECODE']:
        geneDecode = cmdResults['DECODE']
        if not os.path.exists(geneDecode):
            print ("\n\t*** The %s gene Decode file could not be found! ***\n\n" %(geneDecode)) 
            cmdparser.print_help()
            sys.exit(1)
        
    if cmdResults['SEQID']:
        seqID = cmdResults['SEQID']
        if not os.path.exists(seqID):
            print ("\n\t*** The %s seqID Decode file could not be found! ***\n\n" %(seqID)) 
            cmdparser.print_help()
            sys.exit(1)            
    
    cwd = os.getcwd()
    # Start processing the data
    data = BarSeq(fastq, geneDecode, seqID, cwd)
    data.processFastq()   
    data.mergeCounts()
    data.writeTable(data.exact_Results,"Exact-Match.table")
    data.writeTable(data.fuzzy_Results, "1-MisMatch.table")
    data.writeMissingGenes()
    
if __name__ == "__main__":
    main()           
