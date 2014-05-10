'''
Created on Feb 16, 2014

@author: javi
'''
import re
import copy
import os
import NexusEncoder
from collections import Counter as counter

def load(path):
    print("loading..")
    with open(path) as source:
        nameList = [_f for _f in re.split("\t",source.readline().replace("\r\n", "").replace("\n", "")) if _f] #gets names, split, filter for preceeding tabs. Name will be in order.
        rawData = re.split("\n", source.read())
        tree = {}
        chrBranch = {}
        posBranch = {}
        prevLineSplit = None
            
        for line in rawData[:-1]:
            lineSplit = re.split("\t", line.replace("\n", "").replace("\r", ""))
            if lineSplit[0] not in tree:
                chrBranch = {}
                tree[lineSplit[0]] = chrBranch
            
            if int(lineSplit[1]) not in chrBranch:
                posBranch = {}
                chrBranch[int(lineSplit[1])] = posBranch
            
            if isDrift(lineSplit):
                continue
            
            if hasDrift(lineSplit) and prevLineSplit:
#                 this bit determines what to do with "drift". continue to skip line. use the
#                 correctDrift to correct to previous match. This also kind of assumes that
#                 drift doesn't really happen at positions with real SNPs. 
                 
                lineSplit = correctDrift(lineSplit, prevLineSplit)
                
            for snp, samp in zip(lineSplit[3:], nameList):
#                print "%s\t%s"%(samp, snp)
                posBranch[samp] = snp
                
            prevLineSplit = lineSplit
                
    print("loading done")
    return (tree, nameList)

'''(list of chars) -> boolean
see if this position is valid: has at least one valid SNPs shared by more than two strains'''
def isDrift(line):
    return sorted(counter(line).items(), key=lambda x: x[1], reverse=True)[1][1] < 3
    
'''(list of chars) -> boolean
see if there are any SNPs which aren't shared by at least 2 other strains.
Those SNPs are considered drifts and shouldn't be included'''
def hasDrift(line):
    snps = line[3:]
    for element in set(snps):
        if snps.count(element) < 3:
            return True
    return False

'''(list of chars) -> list of chars
picks out the positions where "drift" occurs, and corrects them. Their values
will be reassigned according to which strains the targets were similar in the
previous line. If those strains show different SNPs now, we'll use the majority'''
def correctDrift(line, prevLine):
    driftIndices = identifyDrift(line)
    for drift in driftIndices:
        candidates = findSimilar(prevLine, drift)
        majority = findMajority([line[x] for x in candidates])
#         print(line[drift] + " "  + majority)
        line[drift] = majority
    return line

'''(list of chars) -> list of ints
helper function for correctDrift. ids the elements that are drifts by index'''
def identifyDrift(line):
    MININUM_SHARED = 3
    snps = line[3:]
    drifts = []
    indices = []
    for element in set(snps):
        if snps.count(element) < MININUM_SHARED:
            drifts.append(element)
    for index, element in enumerate(line):
        if element in drifts:
            indices.append(index)
    return indices

'''(list of chars, index) -> list of ints
find the indices of the elements that are identical to the
element at index'''
def findSimilar(line, index):
    indices = []
    for otherIndex, element in enumerate(line):
        if element == line[index]:
            indices.append(otherIndex)
    return indices

'''(list of chars) -> string (of a char)
finds the most common element here'''
def findMajority(list):
    candidates = set(list)
    counts = [(list.count(x), x) for x in candidates]
    return sorted(counts, reverse=True)[0][1]
    
def aggregateForNexus(treeTuple):
    print("begin aggregate")
    results = []
    baseBlock = {}
    for name in treeTuple[1]:
        baseBlock[name] = ""
    block = copy.deepcopy(baseBlock)
    
    for chrs in list(treeTuple[0].items()):
        for pos in sorted(chrs[1].items()):
            for name, info in list(pos[1].items()):
                block[name] += info
                last = name
            if len(block[last]) >= 33:
                results.append(block)
                block = copy.deepcopy(baseBlock)
        results.append(block)
    print("aggregate done") 
    return {"Genome" : results}

if __name__ == '__main__':    
    os.chdir("/data/javi/Toxo/64Genomes")
    #OrderedSNPV6 
    treeTup =  load("/data/javi/Toxo/64Genomes/OrderedSNPV6.txt")  
    NexusEncoder.nexusOutput(aggregateForNexus(treeTup))
    
#     with open("density.txt", "w") as densityOutput:
#     	import SnpSorter
#     	SnpSorter.snpDensity(treeTup[0], densityOutput, treeTup[1])
    
    print("end of GriggsLoader script")
            