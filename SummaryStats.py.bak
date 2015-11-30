'''
Created on Dec 13, 2013

@author: javi
'''
import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
import traceback
import copy
from subprocess import call
import itertools

#Recursive Pairwise Comparison. Returns two values as a tuple: (differences, total) after making pairwise comparisons
#of all possible pair given in the symbolList. Intended for use to calculate pi.
def RPC(symbolList):
    if len(symbolList) == 1:
        return (0,0)
    diff = 0
    total = 0
    next = RPC(symbolList[1:])
    
    for item in symbolList[1:]:
        total += 1
        if not symbolList[0] == item:
            diff += 1
    return (next[0] + diff, next[1] + total)

if __name__ == '__main__':
    #directory = "/data/javi/Toxo/BaseData/tempz"
    directory = "/home/javi/testzone/stattest"
    
    filename = "hardresults.txt"
    os.chdir(directory)
    with open(filename, "r") as source:
        #processing header and generating namelist
        totalcount = 0
        header = re.search("Chromosome\tRefPosition\t(.*?)\n", source.readline()).group(1)
        headersplit = re.split("\t", header)
        samplenames = []
        samplecount = 0
        linelength = 33
        for name in headersplit:
            samplenames.append(name)
            samplecount += 1
        
        #This section includes ME49 in the name list. When this script expands, it can turn
        #into whichever genome is used as alignment reference. This is because there is 
        #relationship to ME49, but it's not reflected per se in the data section. 
        #practically, it adds a "-" to every slot for ME49. 
        refname = "ME49"
        samplenames.append(refname)
        #sample count does not increase by 1 because the dash is added manually.
        
        print "reading.."
        #read data
        tree = []
        chrName = ""
        same = False
        first = True
        data = re.findall("(?m)^(.+?)\t([0-9]+?)\t(.*?)$", source.read())
        for line in data:
            if re.search("(?i)chr", line[0]):
                
                #mechanism for changing chromosomes              
                same = (line[0] == chrName)
                if not same:
                    chrName = line[0]  
                    chrTree = {}
                    tree.append((chrName, chrTree)) 
                    
                #processes 1 characters per sample
                position = int(line[1])
                symbols = re.split("\t",line[2])
                #for ME49
                symbols.append("-")
                
                #checks if all symbols are valid.
                for x in range(samplecount):
                    if not re.match("[ATCNG.-]",symbols[x]): 
                        print "error wrong symbol" + symbols[x]

                #performs calculation for pi
                prepi = RPC(symbols)
                print prepi
                pi = float(prepi[0]) / float(prepi[1])
                chrTree[position] = pi
            
        print "writing.."    
        for chr in tree:
            with open("%s_nexus.nex"%chr[0], "w") as output:            
                #write the header portion
                output.write("#Pi and Theta at each SNP site\n")
                output.write("%s\tpi\ttheta\n"%chr[0])
                for e in sorted(chr[1].items()):
                    output.write("%i\t%f\n" % (e[0], e[1]))
        print "done\n"
        
        
        
        
        
    
