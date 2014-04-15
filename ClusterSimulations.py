'''
Created on Apr 7, 2014

@author: Javi
'''

import numpy as np
import MCLCounter
'''The cluster simulation script aims to perform simulation based functions on clustering patterns. 
This category is started with comparing the clustering patterns of random strain-pairs.
However, other similar functins may evolve from this, and should be included in this script.'''


'''(clusterTree, list of names(currently expects 2) -> Tree of numbers, arranged by chromosome
This function determines whether all the included strains cluster together at every position. 
It iterates through all clusters, compares, and gives out either 1 or 0 as a result. 
Output: data for each chromosome is transformed to a list of 1s and 0s. 
Format: {chr:[0, 1...]}'''
def strainPairComparison(clusters, strains):
    results = []
    for block in clusters:
        if identical(block, strains):
            results.append(1)
        else:
            results.append(0)
    return prune(results)
        
'''helper function for strainPairComparison. 
determines if the strains are identical in this block'''
def identical(clusters, strains):        
    for line in clusters:
        if set(strains).issubset(set(line)):
            return True
    return False
    
'''(List) -> list of tuples representing blocks
iterates through a list using a moving window, and deterine whether a particular
section is considered identical or not. Requires at least 5 consecutive identities
to begin an identity block. '''
def prune(numbers):
    numarray = np.array(numbers)
    results = []
    inBlock = False
    windowSize = 10
    enterValue = 10 #miin number of 1s required to enter a block. 
    exitValue = 3 #min number of 0s required to exit the block. Low because we're looking ahead
    tolerance = 2 #to exit a block, the first (tolerance) elements must not be 1s. 
    start = None
    
    
    for x in range(len(numbers) - windowSize):
        value = sum(numarray[x:x+windowSize])
        if inBlock:
            if value <= exitValue and sum(numarray[x:x+tolerance]) <= tolerance:
                inBlock = False
                results.append([start, x])
        else:
            if value >= enterValue:
                inBlock = True
                start = x
    if inBlock: results.append([start, x])
    
    return np.array(results)
    
            

'''(matrix) -> dictionary of key sites: relevant strains
gives a NxN matrix containing identity information between strains, find the important
blocks. Input should be from prune(). 

IDENTITY MATRIX MUST BE A {CHRNAME: 1D ARRAY of the datatype specified in toIdentity
'''
def findKeyBlocks(identityMatrix):
    tolerance = 5
    minMembers = 4
    keyBlocks = {}
    
    
    for name, chr in identityMatrix.items():
    #chr is of shape (#ofpairs, ['data'])
        dataList = chr['data'].tolist()
        data = np.array(dataList)
        assigned = np.zeros_like(data)
        chrKeyBlocks = {}
        keyBlocks[name] = chrKeyBlocks
        for block in data:       
            if not assigned[data==block].any():
                mask1 = np.abs(data[:,0] - block[0]) < tolerance
                mask2 = np.abs(data[:,1] - block[1]) < tolerance
                corrElements = mask1 * mask2
                if np.sum(corrElements) > minMembers:
                    corrList = np.sort(np.where(corrElements)[0])
                    chrKeyBlocks[(block[0], block[1])] = chr['name'][corrList]
                    assigned[corrElements] = 1
    return keyBlocks


'''(clusterTree) -> keyBlocks
Essentially the runner function to go from clusters to keyBlocks. I imagine there'll be
more to this script so this seems to be a good dividing step. Handles the construction 
of the identity matrix as well'''
def toIdentity(clusterTree, strains):            
    identityMatrix = {}
    datatype = np.dtype([('name', list), ('data', list)])
    for name, chr in clusterTree.items():
        print("processing {0}..".format(name))
        thisChr = []
        for x in range(len(strains)-1):
            for item in strains[x+1:]:
#                 print(" ".join([strains[x], item]))
                twoStrains = strainPairComparison(chr, [strains[x], item])
                for block in twoStrains:
                    thisChr.append(([strains[x], item], block))
        chrarr = np.array(thisChr, dtype=datatype)
        identityMatrix[name] = chrarr
    return findKeyBlocks(identityMatrix)


def recordKeyBlocks(keyBlocks, outfile):
    print("recording...")
    with open(outfile, "wb") as output:
        for chrName, chr in sorted(keyBlocks.items()):
            output.write(bytes("@{0}\n".format(chrName)).encode('utf-8'))
            for name, strains in sorted(chr.items()):
#                 strainList = []
#                 for pair in strains:
#                     for item in pair:
#                         if item not in strainList:
#                             strainList.append(item)
                output.write(bytes("@({0}, {1})\n{2}\n\n".format(name[0], name[1], "\n".join(sorted(build(strains))))).encode('utf-8'))

def build(pairs):
    lines = []
    for pair in pairs:
        found = False
        for line in lines:
            if found: break
            if pair[0] in line or pair[1] in line:
                found = True
                if pair[0] not in line:
                    line.add(pair[0])
                elif pair[1] not in line:
                    line.add(pair[1])
        if not found:
            lines.append(set(pair))
            
    return [" ".join(x) for x in lines]
            
    
def toSimilarity(clusterTree, keyBlocks):
    print("converting to similarity mode.. ")
    compactTree = {}
    for chrName, chr in keyBlocks.items():
        compactChr = {}
        compactTree[chrName] = compactChr
        for coords, strains in sorted(chr.items()):
            strainList = []
            for pair in strains:
                for item in pair:
                    if item not in strainList:
                        strainList.append(item)
            for x in range(coords[0],coords[1]):
                compactChr[x] = " ".join(strainList)

    
    print('constructing new tree...')
    resultsTree = {}
    for chrName, chr in clusterTree.items():
        resultsChr = []
        resultsTree[chrName] = resultsChr
        for index, block in enumerate(chr):
            if index in compactTree[chrName]:
                resultsChr.append(compactTree[chrName][index])
            else:
                resultsChr.append("A")
    return resultsTree
    
    
                
'''(clusterTree) -> Dict of key sites: relevant strains
This is the main function for the simulation. Interesting 
sites where many strains share identity will be picked out and
listed, along with the relevant strains.'''
def simulation(clusters):
    pass


if __name__ == '__main__':
#     filepath = '/home/javi/testzone/Griggs Stuff/persistentResults.txt'
#     tabpath = '/home/javi/testzone/Griggs Stuff/persistentMatrix.tab'
#     outpath = '/home/javi/testzone/Griggs Stuff/keyblocks'

    filepath = '/data/javi/Toxo/64Genomes/Counting/persistentResult.txt'
    tabpath = '/data/javi/Toxo/64Genomes/Counting/persistentMatrix.tab'
    outpath = '/data/javi/Toxo/64Genomes/Counting/keyblocks.txt'

#     filepath = r"F:\Documents\ProjectData\testzone\Griggs Stuff\persistentResults.txt"
#     tabpath = r"F:\Documents\ProjectData\testzone\Griggs Stuff\persistentMatrix.tab"
#     outpath = r"F:\Documents\ProjectData\testzone\Griggs Stuff\keyblocks.txt"
    
#     filepath = r"F:\Documents\ProjectData\64Genomes\Counting\persistentResult.txt"
#     tabpath = r"F:\Documents\ProjectData\64Genomes\Counting\persistentMatrix.tab"
#     outpath = r"F:\Documents\ProjectData\64Genomes\Counting\keyblocks-simple.txt"
    
    clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(filepath, tabpath)[0])
    sampleList = MCLCounter.loadSampleList(tabpath)
    keyblocks = toIdentity(clusterTree, sampleList)
    recordKeyBlocks(keyblocks, outpath)
    
#     import ClusterPattern as cp
#     colorTree = cp.calculateColor(toSimilarity(clusterTree, keyblocks))
#     print('writing..')
#     cp.write(colorTree, outpath)
    
    
    print("End of ClusterSimulations Script")

    