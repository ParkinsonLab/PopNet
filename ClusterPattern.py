'''
Created on Apr 8, 2014

@author: Javi

This class depends on IntactClusters for now because it borrows some functions. 
I'm thinking about refactoring some common functions and reorganizing the whole thing,
but another day maybe?

For now, this script contains functions that all the occurring clustering types between a
small subset of strains in the context of the 64 genomes data. It writes to a file in a
similar format to Similarity and IntactClusters (hence the dependence), and connects to the
ClusterPatternParser on the java side. 
'''

import IntactClusters
import MCLCounter

'''(list of names, list of lines(2D)) -> string
turns the given clustering pattern into a string pattern representing all associations
that are present. Note that the clusters represent
all the lines of a clustering pattern. Tuples are sorted to maintain
consistency'''
def toPattern(strains, clusters):
    pairs = []
    for line in clusters:
        if len(line) > 1:
            pairs += interestingPairs(strains, line)
    
    if len(pairs) == 0: return ""
    pairset = set(pairs)
    stringPairs = ["".join([x[0], x[1]]) for x in pairset]
    pattern = "".join(sorted(stringPairs))
    return pattern        
    
def interestingPairs(strains, line):
    pairs = []
    if len(line)<2:
        return pairs
    
    if line[0] in strains:
        for element in line[1:]:
            if element in strains:
                pairs.append(tuple(sorted((line[0], element))))
    
    return pairs + interestingPairs(strains, line[1:])
        
    
    
def createColorTable(patternTree):
    allPatterns = set()
    maxVal = 16777216
    
    for chrName, chr in list(patternTree.items()):
        for pattern in chr:
            if pattern not in allPatterns:
                allPatterns.add(pattern)
    
    interval = maxVal / len(allPatterns)
    results = {}
    
    for index, pattern in enumerate(allPatterns):
        results[pattern] = int(interval * index)
    
    return results

    
    
'''(Nested list of patterns, by Chr) -> Nested list of color values, by Chr
takes the translated tree and convert each pattern to a value'''
def calculateColor(patternTree):
    colorTable = createColorTable(patternTree)
    newTree = {}
    
    for chrName, chr in list(patternTree.items()):
        newChr = []
        newTree[chrName] = newChr
        for pattern in chr:
            newChr.append(colorTable[pattern])
    
    return newTree

'''(Nested list of clusters, by Chr , list all intersting strains) -> Nested list of patterns, by Chr
translates the clusters into patterns'''
def translate(dataTree, strains):
    newTree = {}
    
    for chrName, chr in list(dataTree.items()):
        newChr = []
        newTree[chrName] = newChr
        for cluster in chr:
            newChr.append(toPattern(strains, cluster))
    return newTree

def write(resultsTree, outfile):
    #something to put on the first part of each line
    selfName = "GENOME"
    itemName = "ME49Cluster"
    
    with open(outfile, "wb") as output:
        #write the hit lit as the keys within the first block 
        output.write(bytes('$\n{0}\n$\n'.format(itemName), 'utf-8'))
        
        #write the rest
        for chrName, chr in resultsTree.items():
            output.write(bytes('%s\n' % chrName, "utf-8"))
            count = 0
            for block in chr:
                output.write(bytes('#%d\n%s - %s: %d\n'%(count, selfName, itemName, block), "utf-8"))
                count+=1
                
if __name__ == '__main__':
    directory = 'F:\Documents\ProjectData\\64Genomes\Counting' 
    filepath = 'F:\Documents\ProjectData\\64Genomes\Counting\persistentResult.txt'
    tabpath = 'F:\Documents\ProjectData\\64Genomes\Counting\persistentMatrix.tab'
    outpath = 'F:\Documents\ProjectData\\64Genomes\Counting\pattern.txt'
    strains = ['ME49', 'RAY', 'PRU', 'ARI', 'B73', 'B41']
    clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(filepath, tabpath)[0])
    results = calculateColor(translate(clusterTree, strains))
    write(results, outpath)
    print("end of ClusterPattern script")