'''
Created on Sep 23, 2013

@author: javi
'''

'''FastA to FastQ!'''
import os
from os import listdir
from os.path import isfile, join
import re
from .FqBlock import FqBlock
from Bio import SeqIO

def format(directory, f):
    with open(f, 'r') as source:
        data = source.read()
    re.findall(">.+?\n(.+)(?:[>]|$)")

if __name__ == '__main__':
     #modify these as needed
    directory = ''
    
    #Do not modify
    os.chdir(directory)  
     
    if not (os.path.isdir(directory + "/Formatted")):
        os.mkdir("Formatted")
    
    chopdirectory = directory + "/Formatted"
        
    log = open(chopdirectory + "/log.txt", "w")
    log.write("Run inputs are: \n%s" % (directory))
    
    onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
    for f in onlyfiles:
        print("\nProcessing %s ..." % f)
        log.write("\nProcessing %s ... " % f)
        
        try: 
            if f.split(".")[1] == "fasta":
                format(directory, f)
            else:
                log.write("\n%s is not a fq file" % f)
        except:
            log.write("\n%s is not a fq file" % f)
 
              
    print("\nend of script")