#!/home/mplace/anaconda3/bin/python

"""
Created on Sat Feb 20 19:47:27 2016

@author: mplace
"""
import sys

bases = ['G', 'C', 'A', 'T' , 'N']
geneTags = dict()

print("tags = {")

with open(sys.argv[1], 'r') as tags:
    for item in tags:
        item = item.rstrip()
        columns = item.split()
        seq = list(columns[1])

        for idx in range(len(seq)):
            newSeq = []
            newSeq = seq[:]   
            for base in bases:
                newSeq[idx] = base
                outSeq = "".join(newSeq)              
                if outSeq not in geneTags:
                    geneTags[outSeq] = columns[0]
                                
newline = 0
for key, value in geneTags.items():
    print("\'%s\':\'%s\',"  %(key, value), end=" ")
    newline += 1
    if newline == 10:
        print()
        newline = 0
            

print("}")

print(len(geneTags))