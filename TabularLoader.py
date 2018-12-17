'''
Created on Feb 16, 2014

@author: javi
'''
import re
import copy
import os
import NexusEncoder
from collections import Counter as counter
import ChrTranslator as ct
import pandas
import numpy as np


def loadToPandas(path, reference, organismm, excludePath=None):
    
    #load df
    df = pandas.DataFrame.from_csv(path, sep='\t')
    
    #filter
    msk = df.apply(filter, axis = 1)
    df_filtered = df.loc[msk]
    
    if not reference == None and df.columns[0] == 'REF':
        df.rename(index='str', columns={'REF': reference})
    
    return df_filtered, df.columns

def filter(x):
    '''
    x is the row. Return 1 if you want to keep.
    '''

    #filter for drift, i.e. only 1
    #Is the oTHER WAY: FALSE is DRIFT to generate a mask for loc
    def filterDrift(x):
        s = list(set(x))
        if len(s) > 2:
            return True
        elif len(s) < 2:
            return False
        else:
            if len(x[x == s[0]]) == 1 or len(x[x == s[1]]) == 1:
                return False
            else:
                return True
    
    return filterDrift(x)        

    
    
#     
# 
# def load(path, reference, organism, excludePath=None):
#     print("loading..")
#     with open(path) as source:
#         nameList = [_f.upper() for _f in re.split("\t",source.readline().replace("\r\n", "").replace("\n", ""))[2:] if _f] #gets names, split, filter for preceeding tabs. Name will be in order.
#         rawData = re.split("\n", source.read())
#         tree = {}
#         chrBranch = {}
#         
#         posBranch = {}
#         prevLineSplit = None
#         accessory = 0
#         
#         
#         
#         if excludePath is not None:
#             with open(excludePath, 'r') as tempfile:
#                 exclude = re.split("\n", tempfile.read().upper())
#         else:
#             exclude = None
#             
#         if not reference == 'None':
#             nameList.pop(0)
#             nameList.insert(0,reference)
# 
#         grey_dict = {n:[] for n in nameList}
#         for index, line in enumerate(rawData[:-1]):
#                 
#             lineSplit = re.split("\t", line.replace("\n", "").replace("\r", ""))
#             chr = ct.translate(lineSplit[0].upper(), mode=organism)
#             
# #             if isAccessory(line): 
# #                 accessory += 1
# #                 continue
#             
#             if chr not in tree:
#                 chrBranch = {}
#                 tree[chr] = chrBranch
#             
#             pos = int(lineSplit[1])
#             if pos not in chrBranch:
#                 posBranch = {}
#                 chrBranch[pos] = posBranch
#             
# #             if isDrift(lineSplit):
# #                 continue
#             
# #             if hasDrift(lineSplit) and prevLineSplit:
# # #                 this bit determines what to do with "drift". continue to skip line. use the
# # #                 correctDrift to correct to previous match. This also kind of assumes that
# # #                 drift doesn't really happen at positions with real SNPs. 
# #                  
# #                 lineSplit = correctDrift(lineSplit, prevLineSplit)
#             
#                 
#             for snp, samp in zip(lineSplit[2:], nameList):
# #                print "%s\t%s"%(samp, snp)
#                 if not (exclude and samp in exclude):
#                     posBranch[samp] = snp
#                 
#                 if snp == '-':
#                     grey_dict[samp].append((chr, pos))
#                 
#             prevLineSplit = lineSplit
#     
#     if exclude:
#         prunedNameList = [x for x in nameList if x not in exclude]
#         return (tree, prunedNameList)
#             
#     print('{} accessory SNPs detected'.format(accessory))
#     print("loading done")
#     return (tree, nameList, grey_dict)
# 
# def isAccessory(line):
#     for e in line:
#         if e == '-':
#             return True
#             
#     return False
# 
# '''(list of chars) -> boolean
# see if this position is valid: has at least one valid SNPs shared by more than two strains
# expecting a numpy array
# '''
# def isDrift(line):
#     min_strains = 2
#     
# #     return sorted(counter(line).items(), key=lambda x: x[1], reverse=True)[1][1] <= minStrains #this last value is the number of required strains + 1
#     #new version that uses numpy:
#     return np.partition(line, -2)[-2] >= min_strains
#     
# '''(list of chars) -> boolean
# see if there are any SNPs which aren't shared by at least 2 other strains.
# Those SNPs are considered drifts and shouldn't be included'''
# def hasDrift(line):
#     minStrains = 2
#     snps = line[3:]
#     for element in set(snps):
#         if snps.count(element) <= minStrains:
#             return True
#     return False
# 
# '''(list of chars) -> list of chars
# picks out the positions where "drift" occurs, and corrects them. Their values
# will be reassigned according to which strains the targets were similar in the
# previous line. If those strains show different SNPs now, we'll use the majority'''
# def correctDrift(line, prevLine):
#     driftIndices = identifyDrift(line)
#     for drift in driftIndices:
#         candidates = findSimilar(prevLine, drift)
#         majority = findMajority([line[x] for x in candidates])
# #         print(line[drift] + " "  + majority)
#         line[drift] = majority
#     return line
# 
# '''(list of chars) -> list of ints
# helper function for correctDrift. ids the elements that are drifts by index'''
# def identifyDrift(line):
#     MININUM_SHARED = 3
#     snps = line[3:]
#     drifts = []
#     indices = []
#     for element in set(snps):
#         if snps.count(element) < MININUM_SHARED:
#             drifts.append(element)
#     for index, element in enumerate(line):
#         if element in drifts:
#             indices.append(index)
#     return indices
# 
# '''(list of chars, index) -> list of ints
# find the indices of the elements that are identical to the
# element at index'''
# def findSimilar(line, index):
#     indices = []
#     for otherIndex, element in enumerate(line):
#         if element == line[index]:
#             indices.append(otherIndex)
#     return indices
# 
# '''(list of chars) -> string (of a char)
# finds the most common element here'''
# def findMajority(list):
#     candidates = set(list)
#     counts = [(list.count(x), x) for x in candidates]
#     return sorted(counts, reverse=True)[0][1]
#     
# def aggregateForNexus(treeTuple):
#     print("begin aggregate")
#     results = []
#     baseBlock = {}
#     for name in treeTuple[1]:
#         baseBlock[name] = ""
#     block = copy.deepcopy(baseBlock)
#     
#     for chrs in list(treeTuple[0].items()):
#         for pos in sorted(chrs[1].items()):
#             for name, info in list(pos[1].items()):
#                 if not re.match('[AGCTN-]$', info):
#                     print("Illegal Character {3} at {0}, {1}, {2}".format(name, str(pos[0]), chrs[0], info))
#                     block[name] += 'N'
#                 else:
#                     block[name] += info
#                 last = name
#             if len(block[last]) >= 33:
#                 results.append(block)
#                 block = copy.deepcopy(baseBlock)
#         results.append(block)
#     print("aggregate done") 
#     return {"Genome" : results}
# 
# def convertGreyDict(grey_dict, blocksize):
#     '''simply converts the position of the grey dict into block numbers according to block size'''
#     new_dict = {}
#     
#     for sample in grey_dict:
#         new_dict[sample] = []
#         for item in grey_dict[sample]:
#             chr = item[0]
#             pos = item[1]
#             result = (chr, pos // blocksize)
#             if result not in grey_dict[sample]:
#                 new_dict[sample].append(result)
#     
#     return new_dict

                
    
if __name__ == '__main__':    
#     os.chdir("/data/javi/Toxo/64Genomes")
#     #OrderedSNPV6 
#     treeTup =  load("/data/javi/Toxo/64Genomes/OrderedSNPV6.txt")  
#     NexusEncoder.nexusOutput(aggregateForNexus(treeTup))
#     
# #     with open("density.txt", "w") as densityOutput:
# #     	import SnpSorter
# #     	SnpSorter.snpDensity(treeTup[0], densityOutput, treeTup[1])

    
    print("end of GriggsLoader script")
            
