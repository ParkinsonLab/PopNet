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

def nexusOutput(tree, samplenames):
    samplecount = len(samplenames)    
    print("writing..")    
    for chr in tree:
        with open("%s_nexus.nex"%chr[0], "w") as output:            
            #write the header portion
            output.write("#NEXUS\nBEGIN taxa;\n\tDIMENSIONS ntax=%i;\nTAXLABELS\n"%(samplecount+1))
            for name in samplenames:
                output.write("[%i]    %s\n"%(samplenames.index(name), name))
            output.write(";\nEND;\nBEGIN characters;\n\tDIMENSIONS nchar=%i;\n\tFORMAT\n\t\tdatatype=standard\n\t\tmissing=.\n\t\tsymbols=\"A T G C N - \"\n\t\tgap=?\n\t\tlabels\n\t\tinterleave\n\t;\n\tMATRIX\n"%chr[1])
            for tb in chr[2]:
                for x in range(samplecount+1):
                    output.write("\t" + samplenames[x] + "\t" + tb[samplenames[x]] + "\n")
                output.write("\n")
            output.write("\n")
            output.write("\t;\nEND;")
    print("done\n%i"%totalcount)

def nexusMatrixOutput(matrix, sampleList):
    samplecount = len(sampleList)    
    print("writing..")
    for chr in matrix:
        data = matrix[chr]
        with open("%s_matrix_nexus.nex"%chr, "w") as output:
            output.write("#NEXUS\nBEGIN taxa;\n\tDIMENSIONS ntax=%i;\nTAXLABELS\n"%(samplecount))
            #for the taxa names thinga at the top.
            for index, name in enumerate(sampleList):
                output.write("[%i]    %s\n"%(index, name))
            output.write(";\nEND [taxa];\n\nBEGIN distances;\nDIMENSIONS ntax=%i;\
            \nFORMAT labels diagonal triangle=both;\nMATRIX\n"%samplecount)
            #the actual data
            for index, x in enumerate(sampleList):
                output.write("[%i] %s\t" % (index, x))
                for y in sampleList:
                    output.write("%i "%data[x][y])
                output.write("\n")
            output.write(";\nEND [distances];")
    print("done")

def clusterMatrixOutput(matrix, sampleList):
    samplecount = len(sampleList)    
    print("writing..")
    for chr in matrix:
        data = matrix[chr]
        with open("%s_matrix.txt"%chr, "w") as output:
            output.write("Strains\t")
            #for the taxa names thinga at the top.
            for index, name in enumerate(sampleList):
                output.write("%s\t"%(name))
            output.write("\n")
            #the actual data
            for index, x in enumerate(sampleList):
                output.write("%s" % (x))
                for y in sampleList:
                    output.write("\t%i"%data[x][y])
                output.write("\n")
    print("done")
    
if __name__ == '__main__':
    directory = "/data/new/javi/yeast/results"
    filename = "results.txt"
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
        #commented out for yeast dataset
#         refname = "ME49"
#         samplenames.append(refname)
        #sample count does not increase by 1 because the dash is added manually.
        
        print("reading..")
        #read data
        tree = []
        block = {}
        for name in samplenames:
            block[name] = ""
        charcount = 0
        thisblock = copy.deepcopy(block)
        chrName = ""
        same = True
        first = True
        full = False
        chrTree = []
        data = re.findall("(?m)^(.+?)\t[0-9]+?\t(.*?)$", source.read())
        for line in data:
            if re.search("(?i)chr", line[0]):
                #checks for chr name, and if not the same start a new block. The full check
                #(33 char per line) is also performed here.
                #exception made for the first one because headache..
                
                if first:
                    chrName = line[0]
                    first = False
                
                same = (line[0] == chrName)
                if full or not same: 
                    chrTree.append(thisblock)
                    if not same:
                        tree.append((chrName, totalcount, chrTree))
                        totalcount = 0
                        chrTree = []  
                        chrName = line[0] 
                        charcount = 0
                    full = False
                    thisblock = copy.deepcopy(block)  
                    charcount = 0
                    
                #processes 1 characters per sample
                symbols = re.split("\\t",line[1])
                for x in range(samplecount):
                    if not re.match("[ATCNG.-]",symbols[x]): 
                        print("error wrong symbol" + symbols[x])
                    thisblock[samplenames[x]] += symbols[x]
                
                
#                 thisblock[refname] += "-" #adds the "-" for the reference sequence
                
                charcount += 1
                totalcount += 1
                if charcount >= 33:
                    full = True
                    charcount = 0
        chrTree.append(thisblock)
        tree.append((chrName, totalcount, chrTree))
    nexusOutput(tree, samplenames)

        

        
        
        
        
        
    
