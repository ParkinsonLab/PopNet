'''
Created on Dec 10, 2013

@author: javi
'''
from Bio import SeqIO
import os
import sys
from os import listdir
from os.path import isfile, join, isdir

def prefix(name):
    while( not name.endswith(".")):
        name = name[:-1]
    name = name[:-1]
    return name

if __name__ == '__main__':
    directory = "/data/javi/Toxo/consensus"
    os.chdir(directory)
    try:
        onlyfiles = [ f for f in listdir(directory) if (isfile(join(directory,f)) and (f.endswith(".fq"))) ]
    except Exception:
        print("\n%s is not a valid file." % f)
        sys.exit()
    
    "files identified. starting..."
    for f in onlyfiles:
        with open(f, "r") as input:
            outputname = prefix(f) + ".fasta"
            count = SeqIO.convert(f, "fastq", outputname, "fasta")
            print("Converted %i reads in %s" % (count, f))
    print("done.")