'''
Created on May 15, 2014

@author: javi
'''

'''this is intended to build the maximum likelihood similarity contigs, 
with the purpose to eliminate, as much as possible, the effects of drift
on the generated image, and to correctly identify a species' closest relatives

Expected input: dataTree: {chrName:[clusterInfoByIndex]} from GriggsLoader, strainInQuestion
Output: {chrName:[closestNeighbour]}
'''

'''([[str]]) -> [str]
returns the list of strains that the given strain clusters with'''
def getCandidates(lines, strain):
    pass

'''([[str]]) -> (str, int)
finds the longest contig, starting from the first element.
returns the name and length of the contig'''
def build(lines):
    pass

'''([[str]]) -> [[str]]
replaces the names of each strain with the name of the group it belongs to.'''
def renameToGroup(lines, groups):
    pass

'''([lines]) -> [str]
main coordinating function for processing information. One chr at a time'''
def processChromosome(clusterInfoList):
    pass


if __name__ == '__main__':
    pass