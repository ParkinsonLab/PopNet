'''
Created on Apr 7, 2014

@author: Javi
'''
'''The cluster simulation script aims to perform simulation based functions on clustering patterns. 
This category is started with comparing the clustering patterns of random strain-pairs.
However, other similar functins may evolve from this, and should be included in this script.'''


'''(clusterTree, list of names(currently expects 2) -> Tree of numbers, arranged by chromosome
This function determines whether all the included strains cluster together at every position. 
It iterates through all clusters, compares, and gives out either 1 or 0 as a result. 
Output: data for each chromosome is transformed to a list of 1s and 0s. 
Format: {chr:[0, 1...]}'''
def strainPairComparison(clusters, strains):
    for chr in clusters:

'''helper function for strainPairComparison. 
determines if the strains are identical in this block'''
def identical(clusters, strains):        
    for line in clusters:
        if set(strains).issubset(set(line)):
            return True
    
'''(List) -> list of tuples representing blocks
iterates through a list using a moving window, and deterine whether a particular
section is considered identical or not. Requires at least 5 consecutive identities
to begin an identity block. '''
def prune(numbers):
    pass

'''(matrix) -> dictionary of key sites: relevant strains
gives a NxN matrix containing identity information between strains, find the important
blocks. Input should be from prune(). '''
def findKeyBlocks(identityMatrix):
    pass

'''(clusterTree) -> Dict of key sites: relevant strains
This is the main function for the simulation. Interesting 
sites where many strains share identity will be picked out and
listed, along with the relevant strains.'''
def simulation(clusters):
    pass




if __name__ == '__main__':
    pass

    