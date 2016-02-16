BarSeq-fast.py -f fastq file -d GENE UP_tag Decode file -s Seq_ID UP_tag Decode file

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


SplitBarSeqByUser-fast.py -f fastq file -s Seq_ID Up_Tag file

    Purpose: Create a new fastq file for each users data.

    Input  : BarSeq fastq file, Seq_ID file

        Seq_ID Up_Tag file:
        Sample ID       Seq Index ID    Sequence Tag    User
        1248 T0 A       A1      AATAGGCGCT    kevin 
        1248 T0 B       B1      AGCGTATGTC    maria


    Usage  : SplitBarSeqByUser.py -f data.fastq -s SeqID.txt
  
	To see Python Docs for this program:

	Open python console and enter
	import sys
	sys.path.append('/full/path/to/script')
	import SplitBarSeqByUser
	help(SplitBarSeqByUser)

	See Mike Place for any problems or suggestions.


NOTE:
BarSeq.py and SplitBarSeqByUser.py -- original implementation, very slow, do not use.

