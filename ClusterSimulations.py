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
    return results
        
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
    results = np.array()
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
    return sorted(results)
    
            

'''(matrix) -> dictionary of key sites: relevant strains
gives a NxN matrix containing identity information between strains, find the important
blocks. Input should be from prune(). 

IDENTITY MATRIX MUST BE A {CHRNAME: 1D ARRAY of the datatype specified in toIdentity
'''
def findKeyBlocks(identityMatrix):
    tolerance = 15
    minMembers = 4
    keyBlocks = {}
    
    
    for name, chr in identityMatrix.items():
    #chr is of shape (#ofpairs, ['data'])
        data = np.array(chr['data'].tolist())
        assigned = []
        chrKeyBlocks = {}
        keyBlocks[name] = chrKeyBlocks
        for strain in chr:
            for element in strain['data']:          
                if element not in assigned:
                    mask1 = np.abs(data[:,:,0] - element[0]) < tolerance
                    mask2 = np.abs(data[:,:,1] - element[1]) < tolerance
                    corrElements = mask1 * mask2
                    if np.sum(corrElements) > minMembers:
                        
                        chrKeyBlocks[element] = chr
                        corrList = np.sort(np.where(corrElements)[0])
                        assigned += chr[corrElements]
    return keyBlocks


'''(clusterTree) -> keyBlocks
Essentially the runner function to go from clusters to keyBlocks. I imagine there'll be
more to this script so this seems to be a good dividing step. Handles the construction 
of the identity matrix as well'''
def toIdentity(clusterTree, strains):            
    identityMatrix = {}
    datatype = np.dtype([('name', tuple), ('data', list)])
    for name, chr in clusterTree.items():
        thisChr = []
        for x in range(len(strains)):
            for item in strains[x:]:
                thisChr.append(((x, item), strainPairComparison(clusterTree, strains)))
        chrarr = np.array(thisChr, dtype=datatype)
    identityMatrix[name] = thisChr
    return findKeyBlocks(identityMatrix)

def recordKeyBlocks(keyBlocks, outfile):
    with open(outfile, "wb") as output:
        for chrName, chr in keyBlocks:
            output.write(bytes("@{0}\n".format(chrName)))
            for name, strains in chr.items():
                output.write(bytes("@({0}, {1})\n{2}\n".format(name[0], name[1], "\n".join(strains))))
    
            
'''(clusterTree) -> Dict of key sites: relevant strains
This is the main function for the simulation. Interesting 
sites where many strains share identity will be picked out and
listed, along with the relevant strains.'''
def simulation(clusters):
    pass


if __name__ == '__main__':
    filepath = '/home/javi/testzone/Griggs Stuff/persistentResults.txt'
    tabpath = '/home/javi/testzone/Griggs Stuff/persistentMatrix.tab'
    outpath = '/home/javi/testzone/Griggs Stuff/keyblocks'
    clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(filepath, tabpath)[0])
    sampleList = MCLCounter.loadSampleList(tabpath)
    recordKeyBlocks(toIdentity(clusterTree, sampleList), outpath)
    

    