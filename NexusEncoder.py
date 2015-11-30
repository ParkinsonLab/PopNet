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

def toMCLTab(matrix, outpath):
    print("Writing MCL tab file")
    with open(outpath, "w") as tab:
        for index, name in enumerate(sorted(matrix)): #writes the doms string, as well as the tab file
            tab.write("%d %s\n"%(index, name))                     
    
def toMCL(matrix, outpath):
    print("Encoding matrix to MCL")
    dimension = len(matrix) #dimensions match the number of samples
    with open(outpath, "w") as output:                    
        #The header portion of each matrix
        output.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"%(dimension, dimension))
        domString = "(mcldoms\n"
        keycount = 0
        for index, name in enumerate(sorted(matrix)): #writes the doms string
            domString += "%d "%index                  
            keycount+=1
        domString += "$\n)\n"
        output.write(domString)
            
        #The data portion
        output.write("(mclmatrix\nbegin\n")  
        for xindex, xinfoTuple in enumerate(sorted(matrix.items())):
            currString = "%d\t"%xindex
            for index, infoTuple in enumerate(sorted(xinfoTuple[1].items())): #0 is the xname
                currString += "%i:%f "%(index, infoTuple[1]) #0 is the yname
            currString += "$\n"
            output.write(currString)
        output.write(")\n")
    

def nexusOutput(tree):
    print("Encoding to Nexus format")
    #Tree is {"genome": info}. items returns a list of items, of which index 0 is the genome tuple. [1] is the info table.
    data = list(tree.items())[0][1]
    exampleDict = data[0]
    samplecount = len(exampleDict)
    
    #calculating total length
    lineLength = len(list(exampleDict.items())[0][1])
    lastLineLength = len(list(data[-1].items())[0][1])
    totalLength = lineLength * (len(data) - 1) + lastLineLength    
    print("writing..")
    with open("Genome_nexus.nex", "w") as output:
        #write the header portion
        output.write("#NEXUS\nBEGIN taxa;\n\tDIMENSIONS ntax=%i;\nTAXLABELS\n"%(samplecount))         
        for index, name in enumerate(sorted(exampleDict)):
            output.write("[%i]    %s\n"%(index+1, name))
        output.write(";\nEND;\nBEGIN characters;\n\tDIMENSIONS nchar=%i;\n\tFORMAT\n\t\tdatatype=standard\n\t\tmissing=.\n\t\tsymbols=\"A T G C N \"\n\t\tgap=?\n\t\tlabels\n\t\tinterleave\n\t;\n\tMATRIX\n"%totalLength)
        for table in data:
            for key, info in sorted(table.items()):
                output.write("\t" + key + "\t" + info + "\n")
            output.write("\n")
        output.write("\t;\nEND;")
    print("done")

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
        with open("%s_matrix.cluster"%chr, "w") as output:
            output.write("Strains\t")
            #for the taxa names thinga at the top.
            for index, name in enumerate(sampleList):
                output.write("%s\t"%(name))
            output.write("\n")
            #the actual data
            for index, x in enumerate(sampleList):
                output.write("%s\t" % (x))
                for y in sampleList:
                    output.write("%i\t"%data[x][y])
                output.write("\n")
    print("done")

def toNexusTree(tree, sampleList):
    baseSampleDict = {}
    for samp in sampleList:
        baseSampleDict[samp] = ""
    sampleDict = copy.deepcopy(baseSampleDict)
    results = []
    for chr in tree:
        for position in chr:
            for key, info in list(chr[position].items()):
                sampleDict[key] += info
            if len(list(sampleDict.items())[0]) >= 33:
                results.append(sampleDict)
                sampleDict = copy.deepcopy(baseSampleDict)
            
    
if __name__ == '__main__':
    directory = "/data/javi/Toxo/BaseData/tempz"
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
            
            
#             thisblock[refname] += "-" #adds the "-" for the reference sequence
            
            charcount += 1
            totalcount += 1
            if charcount >= 33:
                full = True
                charcount = 0
        chrTree.append(thisblock)
        tree.append((chrName, totalcount, chrTree))
    nexusOutput(tree, samplenames)

        

        
        
        
        
        
    
