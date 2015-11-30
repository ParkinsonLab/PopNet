import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
import traceback
import copy
from subprocess import call
import itertools


''' 
SimilarityParser aims to record the degree of similarity between one chosen target and three other chosen types. 
In this case, the numbers are fixed (one and three) due to the coloring method it's designed to work with. There is no
intrinsic limit, I suppose, but if the hitList contains more than three things you won't be able to visualize properly later on.
This script workes off the persistentMatrix.txt generated by the SNPSorter script as well as the tab file. You actually 
need to retain a copy of the tab file and then rename is tab.tab. Please refer to that for further documentation on the input file.

In general, the output file will look like the persistent results generated by SNPSorter (i.e. blocks with identical results are not condensed, and
stats are not recorded). This is again due to the way I've decided to visualize it later. In fact, there is no longer any need to condense those in the first place,
but it's unnecessary to change it either. 

To run this script, the variables at the beginning will be set. It will only process that one file designated by fname. The
results are recorded in a new file. name depends on the input file name, preceeded by the target name and followed by ".sim".
This script must be examined manually before running to make sure you understand what's going on, because it's not all that
robust and there are caveats.
'''
#specify the matrix folder
path = ""
fname = ""
tabName = ""
    
#The names here Must Fully match the names present in the tab file. You
#Might as well copy them over. 
typeI = ""
typeIII = ""
typeII = ""
    
hitList = [typeI, typeIII, typeII]
target = ""
    
def parse():
    #Convert converts the specified names in target and hitList into the column numbers on the matrix.
    #This way, we know what we are looking for inside the matrix. 
    if typeI == "" or path == "":
        print("Error: hit list or path not specified.")
        sys.exit()
    
    convert = {}
    os.chdir(path)
    with open(fname, "r") as infile, open("%s.sim"%(target), "w") as outfile, open(tabName, "r") as tabs:
        #Separates the file by chromosome. The processing logic goes by the tree used in SNP Sorter. 
        #Things are sorted by chromosome, then matrix, then line. 
        chrs = re.findall("(?s)(@.+?)\n(.+?)\n(?=$|@)", infile.read())
        tabLines = re.split("\n", tabs.read())
        tabLines = tabLines[:-1] 
        #The last line contains a newline that I didn't see how to get rid of, and split will generate an empty last element
        #Hence I just got rid of the last element. 
        
        #This section writes the beginning bit listing the targets
        #The order listed here will be the same as the order of the entries
        #for each matrix. It's mainly used for coloring at the next step.
        outfile.write("$\n")
        for hit in hitList:
            outfile.write("%s\n"%hit)
        outfile.write("$\n")
        
        #Convertion dictionary is populated.
        for line in tabLines:
            lineSplit = re.split(" ", line)
            convert[lineSplit[0]] = lineSplit[1]
        
        #The data processing part. 
        for chr in chrs:
            count = 0
            #Writes in the chr name
            outfile.write("%s\n"%chr[0])
            
            #separate the matrices with in chromosome chunk
            matrices = re.findall("(?s)begin\n(.+?)\n[)]", chr[1])
            for matrix in matrices:
                #split each matrix into lines
                outfile.write("#%d\n"%count)
                matrixTable = {}
                matrixSplit = re.split("\n", matrix)
                for line in matrixSplit:
                    #identify the header, and then the rest. 
                    #this format is used because it's not certain how many lines there will be
                    #and everything excepet the header is in a regular format. 
                    #I could probably go a little fancier in the regex but I've decided not to
                    #because I don't like them all that much. 
                    lineSplit = re.match("(?m)^(\d+?)\t(.+?)[$]", line)
                    header = convert[lineSplit.group(1)]
                    matrixTable[header] = {}
                    info = re.findall("(\d+?)[:]([0-9.]+?) ", lineSplit.group(2))
                    for cell in info:
                        #builds this nested dictionary containing all the information within the matrix
                        matrixTable[header][convert[cell[0]]] = cell[1]
                for type in hitList:
                    tally = 0
                    for strain in type:
                        tally += float(matrixTable[target][strain])
                    tally = tally / len(type)
                    #looks for specific points with in specific lines
                    #depending on the target and hitList. 
                    outfile.write("%s - %s: %s\n"%(target, type, tally))
                #this count here supplies the block number. It ends up being useless but oh well. 
                count+=1
    print("done!")


def setTypes(hitList_param, target_param): #(type_1_name, type_2_name, type_3_name), target_name
    global hitList
    hitList = hitList_param
    global target
    target = target_param

def setLocation(folder, file, tab): #"path_to_folder", "data file name", "tab file name"
    global path 
    path = folder
    global fname
    fname = file
    global tabName 
    tabName = tab

def loadToColors(path):
    results = {}
    filePattern = '(?s)^[$]\n(.+?)[$](.+)$'
    linePattern = '^(.+?)\s-\s(.+?):\s(.+)$'

    for filename in os.listdir(path):
        if str(filename).endswith(".sim"):
            with open("{}/{}".format(path, filename), 'r') as input:
                strain = {}
                #the file name of each sim file must be strainName.sim. Fix this manually for now.
                strName = filename[:-4]
                data = re.match(filePattern, input.read()).groups()[1]
                #by looking at the file I see that it's Type I, III, and then II
                chrs = re.split('@', data)[1:]
                for chr in chrs:
                    blocks = re.split('#', chr)
                    colors = []
                    chrName = "@" + blocks[0].rstrip("\n")
                    
                    for block in blocks[1:]:
                        lines = re.split("\n", block)
                        blockNum = lines[0]
                        
                        red = int(float(re.match(linePattern, lines[1]).groups()[2]) * 255)
                        green = int(float(re.match(linePattern, lines[2]).groups()[2]) * 255)
                        blue = int(float(re.match(linePattern, lines[3]).groups()[2]) * 255)
                        
                        hexPattern = '{:0>2X}'*3
                        hex = '#' + hexPattern.format(red, green, blue)
                        colors.append(hex)
                    
                    strain[chrName] = condense(colors)
                results[strName] = aggregate(strain)
    return results

def condense(list):
    chrBlocks = []
    color = ''
    length = 0
    
    for block in list:
        if block is not color:
            if length > 0:
                chrBlocks.append((0, length, color))
            color = block
            length = 1
        else:
            length += 1
    
    chrBlocks.append((0, length, color))
    return chrBlocks
                        
def aggregate(matrix):
    results = [(0, 5, '#000000')]
    import ChrNameSorter as cns
    for name, chr in sorted(matrix.items(), key=lambda x: cns.getValue(x[0])):
        results.append((0, 1, '#000000'))
        results += chr
    return results
                                

    
if __name__ == '__main__':
    path = "/data/new/javi/yeast/results/matrix"
    fname = "persistentMatrix.txt"
    tabname = "persistentMatrix.tab"
    setLocation(path, fname, tabname)
    
        
    global typeI 
    typeI = ['PW5', 'KYOKAI7', 'YJM269', 'YPS163', 'T7', 'EC9-8', 'Y10', 'UC5', 'ZTW1']
    global typeII 
    typeII = ['BY4742', 'BY4741', 'CLIB324', 'W303', 'FL100', 'CEN.PK113-7D', 'SIGMA1278B', 'S288C']
    global typeIII 
    typeIII = ['VIN13', 'T73', 'CBS7960', 'RM11-1A', 'EC1118', 'JAY291', 'AWRI1631', 'M22', 'CLIB215', 'VL3', 'LALVINQA23', 'AWRI796', 'FOSTERSO', 'FOSTERSB', 'YJM789', 'CLIB382']
    
    types = ['BY4742', 'BY4741', 'CLIB324', 'W303', 'FL100', 'CEN.PK113-7D', 'SIGMA1278B', 'S288C', 'PW5', 'KYOKAI7', 'YJM269', 'YPS163', 'T7', 'EC9-8', 'Y10', 'UC5', 'ZTW1'\
, 'VIN13', 'T73', 'CBS7960', 'RM11-1A', 'EC1118', 'JAY291', 'AWRI1631', 'M22', 'CLIB215', 'VL3', 'LALVINQA23', 'AWRI796', 'FOSTERSO', 'FOSTERSB', 'YJM789', 'CLIB382']
    for type in types:
        print('parsing for type {}'.format(type))
        setTypes([typeI, typeIII, typeII], type)
        parse()
    print("Done!")
        
        
                    
                        
                    
                        