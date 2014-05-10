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
from subprocess import call
import itertools





#dataTree will be organized as  a nested dictionary:
#dataTree{'Chromosome name':{'Nucleotide position':{'Sample name':'The SNP'}}}
#Method inputs: File f, Dictionary dataTree. File f should end in either .snps or .vcf. 
def addData(f, dataTree, minCoverage):
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
    
    
#Parsing Data
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
        temp = organizeLine(line, f.name, minCoverageData)
        if temp is not None:
            parsed.append(temp)
      
#Add the parsed data to the dataTree      
    for dataPoint in parsed:
         
        #Ensures that the branch this data line refers exists.
        #If not, create it 
        currentLevel = dataTree
        for branchName in dataPoint[:-2]:
            if branchName not in currentLevel:
                currentLevel[branchName] = {}
            currentLevel = currentLevel[branchName]
        #Adds the SNP value to the tree
        currentLevel[dataPoint[-2]] = dataPoint[-1]
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
def organizeLine(rawLine, name, minCoverage):
    
    lineSplit = re.split("\s+",rawLine)
    try:
        if name.endswith(".snps"):
            return [lineSplit[14].upper(), int(lineSplit[1]), name[:5].upper(), lineSplit[3].upper()]
        else:
            indel = (len(lineSplit[3]) > 1 or len(lineSplit[4]) > 1)
            hetero = re.search("1/1", lineSplit[9])
            quality = (lineSplit[0], lineSplit[1]) in minCoverage
            if hetero and not indel and not quality:
                return [lineSplit[0].upper(), int(lineSplit[1]), name[:5].upper(), lineSplit[4].upper()] 
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
                
def snpDensity(dataTree, output, sampleList):

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
            newList = [copy.deepcopy(posDict) for x in range(length/100000 + 1)] #snps per 100,000 base pairs
            for position in chrBranch:
                group = position/100000
                posBranch = chrBranch[position]
                newBranch = newList[group]
                for sample in sampleList:
                    if sample in posBranch:
                        newBranch[sample] += 1
            
            for branch in newList:
                temppos = newList.index(branch)  * 100000
                posString = chrString + "\t%d"% (temppos)
                for sample in sampleList:
                    posString += "\t%d"%branch[sample]
                posString += "\n"
                output.write(posString)
                
#calculates how many snps are shared betweenthe strains.
#uses the datatree and compared to all other results at that position
#sampleList is modified to contain ME49 at position 0. 
#Anything that                  
def calculateMatrix(dataTree, sampleList):
    
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
            index = int(sorted(chrBranch.keys())[-1])/10000 #This values specifies 1 matrix per 10kb.  
            for x in range(index+1):
                chrMatrix[x] = copy.deepcopy(simpleMatrix)
            
            
            for position in chrBranch:
                chrMatrixBranch = chrMatrix[position/10000] #This one matches the above number
                posBranch = chrBranch[position]
                for x in modSampleList:
                    if x in posBranch:
                        for y in posBranch:
                            if (posBranch[x] == posBranch[y]):
                                chrMatrixBranch[x][y] += 1
# Use only if having ME49!
#                     else:
#                         for y in modSampleList:
#                             if y not in posBranch:
#                                 chrMatrixBranch[x][y] += 1
    
            matrixDict[chrName] = chrMatrix
    return matrixDict

#This method converts the raw data into the mcl matrix specifications.
#The input is a datatree containing all the information for building the matricies. 
#The matrix is always going to be square
#The given matrix really shouldn't contain a value against itself. (center diagonal should be empty)   
#So now we have a matrix file for every 10kb of the chromosome. In order to preserve the sanity of the folder,
#the scheme is to create a tempfile containing each matrix, send that one off to be analyzed, and the deleted.
#a separate persistent file will be created containing all the information in a single file, from the entire thing. 
def recordMatrix(matrixDict):
    persistentMatrixName = "matrix/persistentMatrix.txt"
    persistentResultName = "matrix/persistentResult.txt"   
    persistentDensityName = "matrix/persistentDensity.txt"
    persistentTabName = "matrix/persistentTab.txt"
    with open(persistentMatrixName, "w") as persistentFile, open(persistentResultName, "w") as persistentResult, open(persistentDensityName, "w") as persistentDensity, open(persistentTabName, "w") as persistentTab:
        for chrName in matrixDict:
            if re.search("(?i)chr", chrName): #Only real chromosomes allowed. 
                persistentFile.write("@%s\n"%chrName)
                persistentResult.write("@%s\n"%chrName)
                persistentDensity.write("@%s\n"%chrName)
                persistentTab.write("@%s\n"%chrName)
                chrBranch = matrixDict[chrName]
                
                fileList = []
                tabList = []
                total = 0
                for index in chrBranch:
                    indexBranch = chrBranch[index]
                    dimension = len(list(indexBranch.keys())) #dimensions match the number of samples
                    #now we make the files in here, one separate file per chromosome per 10kb.
                    matrixPrefix = "matrix/%s.%d"%(chrName, index)
                    with open(matrixPrefix + ".mci", "w") as output, open(matrixPrefix + ".tab", "w") as tab:
                        fileList.append(output.name)
                        tabList.append(tab.name)                        
                        #The header portion of each matrix
                        output.write("(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"%(dimension, dimension))
                        domString = "(mcldoms\n"
                        keycount = 0
                        for key in indexBranch: #writes the doms string, as well as the tab file
                            domString += "%d "%keycount
                            tab.write("%d %s\n"%(keycount, key))
                            persistentTab.write("%d %s\n"%(keycount, key))                        
                            keycount+=1
                        domString += "$\n)\n"
                        output.write(domString)
                        persistentTab.write("\n")
                            
                        #The data portion
                        output.write("(mclmatrix\nbegin\n")
                        xcount = 0
                        
                        #This parts allows the matrix to divide by the first value,
                        #thus normalizing all input to between zero and 1
                        tempPos = list(indexBranch.keys())[1]
                        total = max(float(indexBranch[tempPos][tempPos]), float(1))
                        
                        for x in indexBranch:
                            currString = "%d\t"%xcount
                            branch = indexBranch[x]
                            ycount = 0
                            for y in branch:
                                #this one is for not normalized. 
                                #currString += "%s:%f "%(ycount, branch[y])
                                #use for normalization
#                                currString += "%s:%f "%(ycount, branch[y]/total)
                                currString += "%s:%f "%(ycount, branch[y])
                                ycount+=1
                            currString += "$\n"
                            xcount+=1
                            output.write(currString)
                        output.write(")\n")
                    result = mcl(matrixPrefix)
                    persistentFile.write("#%d\n%s\n"%(index, result[0]))
                    persistentResult.write("#%d\n%s\n"%(index, result[1]))
                    persistentDensity.write("#%d\n%d\n\n"%(index, total))
    analyzeMatrix(persistentResultName)
    analyzeMatrix(persistentDensityName)

#This just calls mcl in shell and runs the basic command for it, using the previously generated matrix as input.  
#iValue is currently fixed! Small modification neede5d for it to produce iterant runs with different values.       
def mcl(matrixPrefix):
    iValue = 11 
    matrixName = matrixPrefix + ".mci"
    tabName = matrixPrefix + ".tab"
#    while (iValue <= 4):    
    modName = "%s.%f.result"%(matrixName, iValue)
    dumpName = "%s.%f.dump"%(matrixName, iValue)
    piValue = 8
    #ASSUMES IVALUE IS AN INT!!
    call(["mcl", matrixName, "-use-tab", tabName, "-I", "%d"%iValue, "-o", modName, "-pi", "%d"%piValue]) #The -I option is the inflation value. Play around with it. 
#    call(["mcxdump", "-icl", modName, "-tabr", tabName, "-o", dumpName])  #Don't need this for now because use-tab is doing the trick. 
#        iValue += 0.5    
    with open(matrixName, "r") as a, open(modName, "r") as b:
        data = a.read()
        result = b.read()
    os.remove("%s/%s"%(os.getcwd(), matrixName))
    os.remove("%s/%s"%(os.getcwd(), modName))
    os.remove("%s/%s"%(os.getcwd(), tabName))
    
    return (data,result)

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
                
                #Debug only!
                #print entry[2], matrixStats[entry[2]]
#                 if entry[2] == "3045.\t3142.\tME49\tTGSKN\tP89.S\tARI.V\nVEG.S\tGT1.S":
#                     print "control hit"
#                 if entry[2] == "3142.\tME49\tARI.V\nP89.S\tVEG.S\tGT1.S\n3045.\tTGSKN":
#                     print "hit zzz", matrixStats[entry[2]];  
                output.write("%d - %d \n%s\n###\n"%(entry[0], entry[1], entry[2]))         
        
        sortedInvStats = sorted({(v,k) for k, v in matrixStats.items()}, reverse=True)
        
        #debug!
        #temp = {(v,k) for k, v in matrixStats.iteritems()}
        #print len(temp)
        #print sortedInvStats, len(sortedInvStats), len(matrixStats)
        for x in range(len(sortedInvStats)):
            output.write("@@%s--\n%d\n"%(sortedInvStats[x][1], sortedInvStats[x][0]))
            
            #debug!
            if sortedInvStats[1] == "3142.\tME49\tARI.V\nP89.S\tVEG.S\tGT1.S\n3045.\tTGSKN":
                print("wrote this!")
            
                                
        
                  
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
                result = mcl(matrixPrefix)
                output.write(matrix[0] + result[1] + "\n")

                  
if __name__ == '__main__':
#modify these as needed
    #Works with all files under a directory. 
    #directory = raw_input("Please specify working directory, in full \n")
    #forceRefs = raw_input("Resort to reference types? 0 for No, 1 for Yes.")
    #directory = "/home/javi/testzone/snps"
    #directory = "/data/javi/Toxo/BaseData/tempz"
    #directory = "data/javi/Toxo/NewSeqs/temp"
    directory = "/data/javi/Toxo/64Genomes/Filtered"
    
    #Use this if using on scinet
#    directory = "/scratch/j/jparkins/xescape/SNPSort"
#    forceRefs = "1"
     
     
    #Move working directory to specified
    os.chdir(directory)  
     
    #Create a log for this run. Always called log.txt so will be replace per run.    
    log = open("log.txt", "w")
    log.write("\nRun inputs are: %s" % (directory))
  
# #Locates all vcf and snps files in the directory and adds them to the dataTree one by one.   
#     try:
#         onlyfiles = [ f for f in listdir(directory) if (isfile(join(directory,f)) and (f.endswith(".snps") or f.endswith(".vcf"))) ]
#     except Exception:
#         print "\n%s is not a valid file." % f
#         sys.exit()
  
# #Attempt to sort the samples so the references are processed first
# #Only when reference types are being mapped to.
# #The Type 1 reference is GT1. Type 3 reference is VEG.
#     temp = []
#     sampleList = []
#     if int(forceRefs) == 1:
#         RefA = ""
#         RefB = ""
#         for g in onlyfiles:
#             currName = g[:5].upper()
#             if re.match("GT1", currName):
#                 RefA = currName
#             elif re.match("VEG", currName):
#                 RefB = currName
#             else:
#                 sampleList.append(currName)
#         sampleList.insert(0,RefB)
#         sampleList.insert(0,RefA)
#     else:
#         for g in onlyfiles:
#             currName = g[:5].upper()
#             sampleList.append(currName)
  
# #fileList should already be sorted, if necessary, for resorting. 
# #vcf files Must be accompanied by a min coverage file. 
#     dataTree = {}   #actually a dictionary 
#     for f in onlyfiles:
#         print "\nProcessing %s ..." % f
#         log.write("\nProcessing %s ... " % f) 
#         if f.endswith("vcf"):
#             with open("%s_coverage.min"%(re.split("\.", f)[0]), "r") as minCoverage, open(f, "r") as data:
#                 dataTree = addData(data, dataTree, minCoverage)
#         else:
#             with open(f, "r") as data:
#                 dataTree = addData(data, dataTree, None)
#         
#     with open("results.txt", "w") as results:
#         record(dataTree, results, sampleList)
#     
#         
#     if int(forceRefs) == 1:
#         print "resorting"#                 this bit determines what to do with "drift". continue to skip line. use the
#                 correctDrift to correct to previous match. This also kind of assumes that
#                 drift doesn't really happen at positions with real SNPs. 
#         resorted = resort(dataTree, sampleList)
#         with open("resortedResults.txt", "w") as resortedResults:
#             record(dataTree, resortedResults, sampleList)
#             
#         print "calculating composition"
#         with open("percents.txt", "w") as percentageOutput:
#             calculateComposition(dataTree,sampleList, percentageOutput)
#             
#         print "calculating density"
#         with open("density.txt", "w") as density:
#             snpDensity(dataTree,density,sampleList)
#             
#         print "generating matrix"
#             
#         matrixDir = directory + "/matrix"
#         if not isdir(matrixDir):    
#             os.mkdir(matrixDir)
#         matrix = calculateMatrix(dataTree, sampleList)
#         recordMatrix(matrix)
         

#         #for reanalyzing Only!
#     matrixDir = directory + "/matrix"
#     os.chdir(matrixDir)
#     print "remcl"
#     remcl()
#     print "Reanalyzing Matrix"
#     analyzeMatrix("persistentResult.txt")
         
#       for griggsloader
    import GriggsLoader
    data = GriggsLoader.load("/data/javi/Toxo/64Genomes/OrderedSNPV6.txt")
    dataTree = data[0]
    sampleList = data[1]
    with open("results.txt", "w") as results:
        record(dataTree, results, sampleList)
          
    matrixDir = directory + "/matrix"
    if not isdir(matrixDir):    
        os.mkdir(matrixDir)
    matrix = calculateMatrix(dataTree, sampleList)
    recordMatrix(matrix)    
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
