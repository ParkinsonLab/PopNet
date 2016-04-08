'''
Created on Oct 11, 2013

@author: javi
'''
import os
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
import traceback
import copy
import subprocess as subp
import itertools
import ChrTranslator as ct
import AutoGrouper as ag
import logging
import StringIO
import numpy as np
import math
from multiprocessing import Pool

#dataTree will be organized as  a nested dictionary:
#dataTree{'Chromosome name':{'Nucleotide position':{'Sample name':'The SNP'}}}
#Method inputs: File f, Dictionary dataTree. File f should end in either .snps or .vcf. 
def addData(f, sampleName, dataTree, minCoverage, reference, organism):
#Add a failsafe for file type
    try:
        if not (f.name.endswith(".snps") or f.name.endswith(".vcf")):
            raise Exception("Invalid file type in %s" % f.name)
    except Exception:
        sys.exit()  
        
#parse that coverage file if it exists
    minCoverageData = None
    if minCoverage:
        minCoverageData = re.findall("(?m)^(.*?)\t(.*?)\t.*?\n", minCoverage.read())
    
#come up with a sample name and type
    if f.name.endswith('.snps'):
        type = 'snps'
    else:
        type = 'vcf'    
#Parsing Data
    
    global disregardCount
    totalCount = 0
    disregardCount = 0
    

    #Remove the header
    try:
        data = f.read()
        noHeader = re.search(fetchRegexPattern(f.name), data).group(1)
    except Exception as e:
        print("header parsing error in %s: %s" % (f.name, str(e)))
        sys.exit()
    
    #Divide the remaining data into lines
    rawLines = re.split("\n", noHeader)
    parsed = []
    #Parsing lines
    for line in rawLines[:-1]:
        totalCount += 1
        temp = organizeLine(line, sampleName, type, minCoverageData, organism)
        if temp is not None and re.search('CHR', temp[0]):
            parsed.append(temp)
        elif temp is None:
            disregardCount += 1
      
#Add the parsed data to the dataTree      
    for dataPoint in parsed:
         
        #Ensures that the branch this data line refers exists.
        #If not, create it 
        currentLevel = dataTree
        for branchName in dataPoint[:2]:
            if branchName not in currentLevel:
                currentLevel[branchName] = {}
            currentLevel = currentLevel[branchName]
        #Adds the SNP value to the tree
        currentLevel[dataPoint[-3]] = dataPoint[-2]
        if reference not in currentLevel:
            currentLevel[reference] = dataPoint[-1]
    
    print("{0!s} disregarded out of {1!s} total".format(disregardCount, totalCount))
    
    return dataTree



def fetchRegexPattern(name):
    print("fetching RegEx for %s" % name)
    if name.endswith(".snps"):
        return "(?s)^.*?=\n(.*)"
    else:
        return "(?sm)^([^#].*)"

#Parsed Data Definition: returns as list of lists with the following items:
#[0]Reference Chromosome name
#[1]Reference Nucleotide position
#[2]Sample name
#[3]Nucleotide value
def organizeLine(rawLine, name, type, minCoverage, organism):
    
    lineSplit = re.split("\s+",rawLine)
    try:
        if type is 'snps':
            chr = ct.translate(lineSplit[14].upper(), organism=organism)#for plasmo aligned from 3D7
            ref = lineSplit[2].upper()
            snp = lineSplit[3].upper()
            indel = len(ref) > 1 or len(snp) > 1 or re.search('[^AGCT]', ref) or re.search('[^AGCT]', snp)
            pos = int(lineSplit[1])
            if not indel:
                return [chr, pos, name, snp, ref]
            else:
                return None
        else:
            chr = ct.translate(lineSplit[0].upper(), organism=organism).upper()
            indel = len(lineSplit[3]) > 1 or len(lineSplit[4]) > 1 or re.search('[^AGCT]', lineSplit[3]) or re.search('[^AGCT]', lineSplit[4])
            hetero = re.search("1/1", lineSplit[9])
            quality = float(lineSplit[5])
            pos = int(lineSplit[1])
            snp = lineSplit[4].upper()
            ref = lineSplit[3].upper()
            #choice of 5 as a quality min is a bit arbitrary
            if hetero and not indel and quality > 0:
                return [chr, pos, name, snp, ref] 
            else:
                return None
    except Exception as e:
        print("Illegal Line in file %s: %s" % (name, e))
        print(lineSplit)
        traceback.print_exc()
        sys.exit()
    
#This outputs the data structure into something we can read.
def record(dataTree, output, fileList):
#Writes a header into the output file    
    headerString = "Chromosome\tRefPosition"
    for headerSample in fileList:
        headerString += "\t%s"%headerSample
    headerString += "\n"
    output.write(headerString)

#Records the rest of the data. Format is Chromsome name    Reference Position    Sample 1 allele.. etc. 
#One line per Reference position that has SNP in at least one file
#Post processing will be done after. In a separate file probably. 

    for chrName in sorted(dataTree):
        chrString = chrName
        chrBranch = dataTree[chrName]
        chrBranchSortedList = sorted(dataTree[chrName])
        for position in chrBranchSortedList:
            posString = chrString + "\t%d"%position
            posBranch = chrBranch[position]
            for sample in fileList:
                if sample in posBranch:
                    posString = posString + "\t%s"%posBranch[sample]
                else:
                    posString = posString + "\t-"
            sampString = posString + "\n"
            output.write(sampString)
                
#Assumes a sorted fileList, with type 1 reference being first and type 3 being second. 
def resort(dataTree, fileList):
    print("resorting..")
    for chrName in dataTree:
        chrBranch = dataTree[chrName]
        
        for position in chrBranch:           
            posBranch = chrBranch[position]
            temp = {}
            for n in fileList:                                      #This section populates temp such that each element in fileList
                if n in posBranch:                                  #is compared with the ones before it. "Genotypes" are assigned based on
                    for p in fileList[:fileList.index(n)]:          #the sample's index position in fileList (GT1 and VEG are always 1 and 2)
                        if p in posBranch:
                            if (posBranch[p]==posBranch[n] and n not in temp):
                                temp[n] = fileList.index(p)
                    if n not in temp:
                        temp[n] = fileList.index(n)
            chrBranch[position] = temp
    return dataTree        

def calculateComposition(resortedData, sampleList, output):
    
    #Writes a header into the output file    
    headerString = "Chromosome\tRefPosition"
    for headerSample in sampleList:
        headerString += "\t%s"%headerSample
    headerString += "\n"
    output.write(headerString)
    
    #The rest of the logic. Calculates the percent of SNPs belonging to each (previous) type.
    #Typing priority will be given to the samples that appear earlier on the list. 
    #Since this requires a resorted list, the references are necessary, and will always be first and second. 
    for chrName in sorted(resortedData):
        currString = chrName
        aggregateData = {}
        chrBranch = resortedData[chrName]
        for name in sampleList:
            aggregateData[name] = []
        
        for position in chrBranch:
            posBranch = chrBranch[position]
            for sample in posBranch:
                aggregateData[sample].append(posBranch[sample])
        
        for name in sampleList:
            stats = []
            if len(aggregateData[name])>0: 
                total = float(len(aggregateData[name]))  
                for x in range(len(sampleList)):
                    try:
                        percent = aggregateData[name].count(x) / total * 100
                    except ZeroDivisionError:
                        percent = 0
                    finally:
                        stats.append(percent)
                sampString = ""
                for n in stats:
                    if not n==0:
                        sampString += "Type %d:%d, "%(stats.index(n), n)
                currString += "\t%s |"%sampString
            else:
                currString += "\tNONE |"
        output.write("%s\n"%currString)
                
def snpDensity(dataTree, outpath, sampleList, section_length):
    with open(outpath, "w") as output:
        #Writes a header into the output file    
        headerString = "Chromosome\tRefPosition"
        for headerSample in sampleList:
            headerString += "\t%s"%headerSample
        headerString += "\n"
        output.write(headerString)

        #create a new tree, define branches
        posDict = {}
        for name in sampleList:
            posDict[name] = 0
                    
        for chrName in dataTree:
            if re.search("(?i)chr", chrName):
                chrString = chrName
                chrBranch = dataTree[chrName]
                length = sorted(chrBranch.keys())[-1]
                newList = [copy.deepcopy(posDict) for x in range(length/section_length + 1)] #snps per 1000 base pairs
                for position in sorted(chrBranch):
                    group = position/section_length
                    posBranch = chrBranch[position]
                    newBranch = newList[group]
                    for sample in sampleList:
                        if sample in posBranch:
                            newBranch[sample] += 1
                
                for index, branch in enumerate(newList):
                    temppos = index  * section_length
                    posString = chrString + "\t%d"% (temppos)
                    for sample in sampleList:
                        posString += "\t%d"%branch[sample]
                    posString += "\n"
                    output.write(posString)
                
#calculates how many snps are shared betweenthe strains.
#uses the datatree and compared to all other results at that position
#sampleList is modified to contain ME49 at position 0. 
#Anything that                  
def calculateMatrix(dataTree, sampleList, section_length):
    
    matrixDict = {}
    modSampleList = sampleList
    #Use only if having ME49 as part of the list!
#    modSampleList.insert(0,"ME49")
    for chrName in dataTree:
        if re.search("(?i)chr", chrName):
            chrBranch = dataTree[chrName]
            
            simpleMatrix = {}
            for x in modSampleList:
                simpleMatrix[x] = {}
                for y in modSampleList:
                    simpleMatrix[x][y] = 0
            chrMatrix = {}
            index = int(sorted(chrBranch.keys())[-1])/section_length #This values specifies 1 matrix per 10kb.  
            for x in range(index+1):
                chrMatrix[x] = copy.deepcopy(simpleMatrix)
            
            
            for position in chrBranch:
                chrMatrixBranch = chrMatrix[position/section_length] #This one matches the above number. and the value from snpdensity
                posBranch = chrBranch[position]
                for x in modSampleList:
                    for y in modSampleList:
                        if x in posBranch and y in posBranch:
                            if (posBranch[x] == posBranch[y]):
                                chrMatrixBranch[x][y] += 1
                        elif x not in posBranch and y not in posBranch:
                                chrMatrixBranch[x][y] += 1
                                
#Transform Matrix
#Put the lowest value as 0 and the highest as 1, rescale all values. 
            for index, branch in chrMatrix.items():
                sampleKey = branch.keys()[0]
                highest = max(branch[sampleKey][sampleKey], 1)
                templist = []
                for x in branch.values():
                    templist += x.values()
                lowest = min(templist)
                if highest == lowest:
                    print("bad matrix at {0} {1}".format(chrName, str(index)))
                for x in branch.keys():
                    for y in branch[x].keys():
                        if branch[x][y] > highest:
                            print("Over limit at {0} {1} {2} {3}".format(chrName, str(index), x, y))
                        branch[x][y] = float(branch[x][y] - lowest) / float(max(highest - lowest, 1))
    
            matrixDict[chrName] = chrMatrix
    return matrixDict


def recordTab(sampleList, tabpath):
    #writes the single tab file
    with open(tabpath, 'w') as persistentTab:
        for index, key in enumerate(sampleList):
            persistentTab.write("{0} {1}\n".format(str(index), key))
                    
###helper for recordMatrix
def mcl_unpack(param_tuple):
    currString = param_tuple[0]
    tabpath = param_tuple[1]
    i = param_tuple[2]
    pi = param_tuple[3]
    
    return currString, ag.mcl(currString, tabpath, i, pi, True)

#This method converts the raw data into the mcl matrix specifications.
#The input is a datatree containing all the information for building the matricies. 
#The matrix is always going to be square
#The given matrix really shouldn't contain a value against itself. (center diagonal should be empty)   
#So now we have a matrix file for every 10kb of the chromosome. In order to preserve the sanity of the folder,
#the scheme is to create a tempfile containing each matrix, send that one off to be analyzed, and the deleted.
#a separate persistent file will be created containing all the information in a single file, from the entire thing. 
def recordMatrix(matrixDict, sampleList, tabpath, persistentMatrixName, persistentResultName, i, pi):
    
    ##helpers
    def generateParams(currString_list, tabpath, i, pi):
        
        return [(currString, tabpath, i, pi) for currString in currString_list]   
    ###
                  
    with open(persistentMatrixName, "w") as persistentFile, open(persistentResultName, "w") as persistentResult:
        print('SNPSorter: Applied I value is {}'.format(i))       
        for chrName in matrixDict:
            if re.search("(?i)chr", chrName): #Only real chromosomes allowed. 
                persistentFile.write("@%s\n"%chrName)
                persistentResult.write("@%s\n"%chrName)
                chrBranch = matrixDict[chrName]
                currString_list = [buildMatrix(chrBranch[pindex], sampleList) for pindex in sorted(chrBranch)]
                param_list = generateParams(currString_list, tabpath, i, pi)
                
                pool = Pool()
                results_list = pool.map(mcl_unpack, param_list)
                    
            for index, (currString, result) in enumerate(results_list):            
                persistentFile.write("#%d\n%s\n"%(index, currString))
                persistentResult.write("#%d\n%s\n"%(index, result))

    analyzeMatrix(persistentResultName)


'''(dict, list) -> string
builds the matrix (.mci) string to be used in mcl
'''
def buildMatrix(tree, sampleList):
    
    currString = ""
    dimension = len(tree) #dimensions match the number of samples
    #now we make the files in here, one separate file per chromosome per 10kb.
     
    #The header portion of each matrix
    currString += "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"%(dimension, dimension)
    currString += "(mcldoms\n"
    for index, key in enumerate(sampleList): #writes the doms string, as well as the tab file
        currString += "%d "%index                      
    currString += "$\n)\n"
        
    #The data portion
    currString += "(mclmatrix\nbegin\n"
    xcount = 0
    
    #This parts allows the matrix to divide by the first value,
    #thus normalizing all input to between zero and 1
    for x in sampleList:
        currString += "%d\t"%xcount
        branch = tree[x]
        ycount = 0
        for y in sampleList:
            #value is pre-normalized during matrix construction
            currString += "%s:%f "%(ycount, branch[y])
            ycount+=1
        currString += "$\n"
        xcount+=1
    currString += ')'
    
    return currString

#Reads a persistentResults file, parses it, and then interprets it. 
#Input is the file name
#Output is a summary of the types of clutering present in the results
def analyzeMatrix(results):
    with open(results, "r") as data, open(results + ".analyzed", "w") as output:
        matrixStats = {}
  
        for chr in re.findall("(?s)(@.*?)\n(.*?\n\n)(?=@|$)", data.read()):            
            output.write("%s\n"%chr[0])
            posData = re.findall("(?s)[#]([0-9]*?)\n(.*?)\n\n(?=#|$)", chr[1])#Pos data is a tuple in the form of (position, matrix) Matrixes are just evaluated as strings because we just
            tree = {}                                                       #need to know whether they are the same or not. 
            for position in posData:
                entry = position[1]
                if entry not in tree:
                    tree[entry] = []
                tree[entry].append(int(position[0]))    
            for entry in arrangeMatrixResults(tree):
                if entry[2] in matrixStats:
                    matrixStats[entry[2]] += (entry[1] - entry[0] + 1)
                else:
                    matrixStats[entry[2]] = (entry[1] - entry[0] + 1) 
                output.write("%d - %d \n%s\n###\n"%(entry[0], entry[1], entry[2]))         
        
        sortedInvStats = sorted({(v,k) for k, v in matrixStats.items()}, reverse=True)
        
        #stats
        total_clusters = 0
        single_sections = 0
        total = 0
        for x in range(len(sortedInvStats)):
            pattern = sortedInvStats[x][1]
            count = sortedInvStats[x][0]
            
            output.write("@@%s--\n%d\n"%(pattern, count))
            
            #calculate some statistics!
            total += 1
            n = len(re.split("\n", pattern))
            total_clusters += n  
            if n == 1: single_sections += 1
        import decimal as dc
        dc.getcontext().prec = 2
        print("Analysis Results:\n Average {0} clusters over {1} sections.\n{2} unclustered region detected, representing {3}% of total".format(\
str(dc.Decimal(total_clusters)/dc.Decimal(total)), str(total), str(single_sections), str(dc.Decimal(single_sections)/dc.Decimal(total))))
            
                                
        
                  
def toRange(numList):
        #This is a bit of a convoluted function, but it basically sorts the list of numbers into a list of tuples representing the ranges
        #for ex. [1,2,4,5,9] -> [(1,2), (4,5), (9,9)]
        #this is then transformed to have one entry per line and returned as a string to be written in the parent function
        return [(t[0][1], t[-1][1]) for t in (tuple(g[1]) for g in itertools.groupby(enumerate(sorted(numList)), lambda i_x: i_x[0] - i_x[1]))]
        
            
def arrangeMatrixResults(tree): 
    result = []
    for entry in tree:   
        for value in sorted(toRange(tree[entry])):
            result.append(value + (entry,))
    return sorted(result)  

#I'm just going to have it require a persistent matrix.txt
def remcl():
    matrixPrefix = "persistentMatrix"
    outputPrefix = "persistentResult"
    with open("%s.txt"%matrixPrefix, "r") as data, open("%s.tab"%matrixPrefix, "r") as tab, open("%s.txt"%outputPrefix, "w") as output:
        for chr in re.findall("(?s)(@.*?)\n(.*?\n\n)(?=@|$)", data.read()):            
            output.write("%s\n"%chr[0])
            allMatrix = re.findall("(?s)(#[0-9]+?\n)([(].*?[)])\n\n(?=#|\Z)", chr[1])
            for matrix in allMatrix:
                with open("%s.mci"%matrixPrefix, "w") as currMatrixFile:
                    currMatrixFile.write(matrix[1])
                result = ag.mcl(matrixPrefix)
                output.write(matrix[0] + result[1] + "\n")

#records densities of all strains
def multiDensity(dataTree, outfile, section_length):
    results = {}
    text = ""
    
    for chrName, chr in dataTree.items():
        text += "@{}\n".format(chrName)
        data = {}
        for position, posBranch in chr.items():
            binNum = int(position / section_length)
            if binNum not in data: data[binNum] = 0
            data[binNum] += 1
            
        results[chrName] = [x[1] for x in sorted(data.items())]
           
#printing only         
        for number, bin in data.items():
            text += str(number) + "\t" + str(bin) + "\n"
               
    with open(outfile, "w") as output:
        output.write(text)
            
    return results

#fills the empty spots in the dataTree with the reference sequence, so it functions like data from Griggloader.
#the reference is guaranteed to be in every position.
def fillDataTree(dataTree, sampleList, reference):

    
    for name, chr in dataTree.items():
#         print("filling data for chromosome {}".format(name))
        for position in chr.values():
            for sample in sampleList:
                if sample not in position:
                    position[sample] = position[reference]
                    
    return dataTree

            
if __name__ == '__main__':
    
# #modify these as needed
#     #Works with all files under a directory. 
#     #directory = raw_input("Please specify working directory, in full \n")
#     #forceRefs = raw_input("Resort to reference types? 0 for No, 1 for Yes.")
#     #directory = "/home/javi/testzone/snps"
# #     directory = "/data/new/javi/yeast/results"
# #     directory = '/data/new/javi/yeast/data'
#     #directory = "/data/new/javi/plasmo/vcfs"
#     #directory = "data/javi/Toxo/NewSeqs/temp"
# #     directory = "/data/javi/Toxo/64Genomes/Filtered"
#     directory = "/data/javi/Toxo/Hybrid"
#     
#     #Use this if using on scinet
# #    directory = "/scratch/j/jparkins/xescape/SNPSort"
#     forceRefs = 0
#      
#      
#     #Move working directory to specified
#     os.chdir(directory)  
#      
#     #Create a log for this run. Always called log.txt so will be replace per run.    
#     log = open("log.txt", "w")
#     log.write("\nRun inputs are: %s" % (directory))
#   
# #Locates all vcf and snps files in the directory and adds them to the dataTree one by one.   
#     try:
#         onlyfiles = [ f for f in listdir(directory) if (isfile(join(directory,f)) and (f.endswith(".snps") or f.endswith(".vcf"))) ]
#     except Exception:
#         print("\n%s is not a valid file." % f)
#         sys.exit()
#     
# # #Attempt to sort the samples so the references are processed first
# # #Only when reference types are being mapped to.
# # #The Type 1 reference is GT1. Type 3 reference is VEG.
# #     temp = []
# #     sampleList = []
# #     if int(forceRefs) == 1:
# #         RefA = ""
# #         RefB = ""
# #         for g in onlyfiles:
# #             currName = re.split(".", g)[0].upper()
# #             if re.match("GT1", currName):
# #                 RefA = currName
# #             elif re.match("VEG", currName):
# #                 RefB = currName
# #             else:
# #                 sampleList.append(currName)
# #         sampleList.insert(0,RefB)
# #         sampleList.insert(0,RefA)
# #     else:
# #         for g in onlyfiles:
# #             currName =re.split(".", g)[0].upper()
# #             sampleList.append(currName)
# 
# #a simpler sample list, tailor to your own needs:
# 
#    
# # fileList should already be sorted, if necessary, for resorting. 
# # vcf files Must be accompanied by a min coverage file. 
#     dataTree = {}   #actually a dictionary
#     sampleList = []
#      
#     reference = "ME49"
#     if reference not in sampleList: sampleList.append(reference)
#       
#     for f in onlyfiles:
#         print("\nProcessing %s ..." % f)
#         log.write("\nProcessing %s ... " % f)
#         #sample name pattern here
#  #        pattern = '^(.+?)[_].*' #for the old plasmodium stuff
#         pattern = '^(.+?)[\.].*'
#          
#         sampleName = re.match(pattern, f).group(1).upper()
#         if sampleName not in sampleList: 
#             if f.endswith("vcf"):
#                 with open("%s_coverage.min"%(re.split("\.", f)[0]), "r") as minCoverage, open(f, "r") as data:
#                     dataTree = addData(data, sampleName, dataTree, minCoverage, reference)
#             else:
#                 with open(f, "r") as data:
#                     dataTree = addData(data, sampleName, dataTree, None, reference)
#             sampleList.append(sampleName)
#         else:
#             print("Duplicate for {0}".format(sampleName))
#     with open("results.txt", "w") as results:
#         record(dataTree, results, sampleList)
#      
#          
# #     if int(forceRefs) == 1:
# #         print "resorting"
# #         resorted = resort(dataTree, sampleList)
# #         with open("resortedResults.txt", "w") as resortedResults:
# #             record(dataTree, resortedResults, sampleList)
# #              
# #         print "calculating composition"
# #         with open("percents.txt", "w") as percentageOutput:
# #             calculateComposition(dataTree,sampleList, percentageOutput)
#              
#     print("calculating density")
#  
#     snpDensity(dataTree,"density.txt",sampleList)
# #    multiDensity(dataTree, "density.txt")
#            
#     import DriftDetection as dd
#     dataTree = dd.scan(dataTree)
#       
#     print('filling data tree')
#     dataTree = fillDataTree(dataTree, sampleList, reference)
#    
#     print("generating matrix")     
#     matrixDir = directory + "/matrix"
#     if not isdir(matrixDir):    
#         os.mkdir(matrixDir)
#     matrix = calculateMatrix(dataTree, sampleList)
#     recordMatrix(matrix)
#       
#     print('encoding to nexus')     
#     import GriggsLoader as gl
#     import NexusEncoder as ne
#     treetuple = (dataTree, sampleList)
#     ne.nexusOutput(gl.aggregateForNexus(treetuple))
    
#         #for reanalyzing Only!
#     matrixDir = directory + "/matrix"
#     os.chdir(matrixDir)
#     print("remcl")
#     remcl()
#     print("Reanalyzing Matrix")
#     analyzeMatrix("persistentResult.txt")
         
#       for griggsloader
    import TabularLoader as gl
    import NexusEncoder as ne
    os.chdir('/data/javi/Toxo/64Genomes/Filtered/matrix')
    data = gl.load("/data/javi/Toxo/64Genomes/OrderedSNPV6.txt", '/data/javi/Toxo/64Genomes/Filtered/matrix/exclude.txt')
    ne.nexusOutput(gl.aggregateForNexus(data))
#testing
#    data = GriggsLoader.load("/data/javi/Toxo/64Genomes/test.txt") 

#     dataTree = data[0]
#     sampleList = data[1]
#     with open("results.txt", "w") as results:
#         record(dataTree, results, sampleList)
#      
#     densityFile = "/data/javi/Toxo/64Genomes/Filtered/matrix/density.txt"
#     snpDensity(dataTree, densityFile, sampleList)
          
#     matrixDir = directory + "/matrix"
#     if not isdir(matrixDir):    
#         os.mkdir(matrixDir)
#     matrix = calculateMatrix(dataTree, sampleList)
#     recordMatrix(matrix)    
    print("\nend of script")



#Iterate through all the files, parse, and merge into a single nested dictionary. The order is, from the top level,
#Chromosome, Position, Source, and Nucleotide. (Delegated to addData())

#Sort the dictionary by Chromosome, then position, then source (Delegated to addData())

#Output the results, and do the necessary calculations (Delegated to record())\


###############################

#No longer needed. Left here for reference. Delete when it works. 
# def parseVCF(f):
#     
#     #First, remove the header information
#     try:
#         noHeader = re.search("(?sm)^([^#].*)", f.open()).group(1)
#     except Exception:
#         print "header parsing error in %s" % f.name
#     
#     #Divide the remaining data into lines
#     rawLines = re.split("\n", noHeader)
#     results = []
#     
#     #Parsing lines
#     for line in rawLines:
#     
# def parseSNPS(f):
#     
#     #First, remove the header information
#     try:
#         noHeader = re.match("(?s)^.*?=\n(.*)", f.read()).group(1)
#     except Exception:
#         print "header parsing error in %s" % f.name
#         sys.exit()
#     
#     #Divide the remaining data into lines
#     rawLines = re.split("/n", noHeader)
#     results = []
#     #Parsing lines
#     for line in rawLines:
#         lineSplit = re.split("\W+",rawLines)
#         results.append([lineSplit[11], lineSplit[1], f.name[:5], lineSplit[3]])
#     
#     return results    
