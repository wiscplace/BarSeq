BarSeq.py -f fastq file -d GENE UP_tag Decode file -s Seq_ID UP_tag Decode file

    Purpose: Produce a table of gene counts from DNA Bar-Seq data

    Input  : BarSeq fastq file, GENE UP_tag Decode file, Seq_ID UP_tag Decode file

    GENE UP_tag Decode file:
        Gene_ID UP_tag_Seq
        YAL001C ACTATATGTGAAGGCATGGC

        Seq_ID UP_tag Decode file:
        Only one user name allowed.
        Sample ID       Seq Index ID    Sequence Tag    User
        1248 T0 A       A1      AATAGGCGCT    kevin 
        1248 T0 B       B1      AGCGTATGTC    kevin


    Usage  : BarSeq-fast.py -f data.fastq -d GENE_UP_tag_Decode.txt -s SeqID.txt
  
	To see Python Docs for this program:

	Open python console and enter
	import sys
	sys.path.append('/full/path/to/script')
	import BarSeq
	help(BarSeq)

	See Mike Place for any problems or suggestions.


SplitBarSeqByUser.py -f fastq file -s Seq_ID Up_Tag file

    Purpose: Create a new fastq file for each users data.

    Input  : BarSeq fastq file, Seq_ID file
             File format is important, user almost never provides
             the decode file in the same format, so you will
             most likely have to manipulate it w/ vim or perl.
             
    example1:
        Seq_ID Up_Tag file:
        Sample ID       Seq Index ID    Sequence Tag    User
        1248 T0 A       A1      AATAGGCGCT    kevin 
        1248 T0 B       B1      AGCGTATGTC    maria

    example2:
        Sample Name XX  Index Tag   Maria
        WA Lib 1    XX  AATAGGCGCT  Maria
        WA Lib 2    XX  TACAGTTGCG  Maria
        PAP Lib 1   XX  ATCCTAGCAG  Maria
        PAP Lib 2   XX  GATTAGCCTC  Maria
    
    Usage  : SplitBarSeqByUser.py -f data.fastq -s SeqID.txt
  
	To see Python Docs for this program:

	Open python console and enter
	import sys
	sys.path.append('/full/path/to/script')
	import SplitBarSeqByUser
	help(SplitBarSeqByUser)

	See Mike Place for any problems or suggestions.

GENE_UpTag_Dictionary.py

    Dictionary for MoBY Sequence UP_tags

SplitBarSeq.py -f fastq file -s Seq_ID Up_Tag file

    Purpose: Create a new fastq file for each sampleID.

    Input  : BarSeq fastq file, Seq_ID file

    Seq_ID Up_Tag file:
    Sample ID       Seq Index ID    Sequence Tag    User
    1248 T0 A       A1      AATAGGCGCT    kevin 
    1248 T0 B       B1      AGCGTATGTC    maria

    Usage  : SplitBarSeq.py -f data.fastq -s SeqID.txt

tag_dictionary.py 

    Purpose: Generate all possible single mismatches for a sequence

    Input  : Gene Up Tag file

GENE_UPTAG_Decode.txt
    
    Purpose:  Provides a list of S.cerevisiae gene names and the associated 
    gene uptags for the MOBY collection.

    format:
    
        Gene_ID UP_Tag_Seq_ID
        YDR225W TGCTCTCCAT
        YIL077C AGTATCACTA
        YML058W TCAGGCTATT
        YMR115W GCATCTGGCT
        YNR059W GGATATATCG
        YDL161W CGGATTCATCT



