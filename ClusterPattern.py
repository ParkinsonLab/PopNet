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
    
    if len(pairs) == 0: return frozenset("None")
    pairset = set(pairs)
    stringPairs = ["".join([x[0], x[1]]) for x in pairset]
    pattern = frozenset(sorted(stringPairs))
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
    import random
    
    
    colorIter = iter([16737894, 1337452, 10040217, 6063411, 6710937, 6037279, 9737364, 16750899, 50637, 16777215])
    #quick fix to get the same colors
    groups = ['GUYS', 'TgH', 'ME49', 'VEG', 'p89', 'TgCats', 'GAL-DOM1/2', 'GAL-DOM10', 'GT1', 'MISC']
    results = {}
    for group in groups:
        results[frozenset([group])] = colorIter.next()
    
#     allPatterns = set()
#     maxVal = 16777215
#     
#     for chrName, chr in patternTree.items():
#         for pattern in chr:
#             if pattern not in allPatterns:
#                 allPatterns.add(pattern)
#     
# #     if len(allPatterns)%2 == 1:
# #         allPatterns.update([frozenset(['FILLER1'])])
# 
#     
#     interval = 148638
#     results = {}
#     
#     for index, pattern in enumerate(allPatterns):
#         try:
#             results[pattern] = next(colorIter)
#         except StopIteration:
#             pass
            
    printColorTable(results)
    return results

'''table is a dict of frozenset:value
just a helper function because i don't like their print'''
def printColorTable(table):
    maxKeyLength = max([len(x) for x in table.keys()])
    maxValueLength = max([len(str(x)) for x in table.values()])
    print(maxKeyLength, maxValueLength)
    rowFormat = "{:<{klen}}\t{:<{vlen}}"
    result = "\n".join([rowFormat.format(x[0], hex(x[1]), klen=maxKeyLength, vlen=maxValueLength) for x in table.items()])
    print(result)
    
'''(dict) -> display
Displays the color table in a readable format, with group names
and their colors'''
def translateColorTable(table):
    import binascii as asc
    hexValues = {}
    for name, value in table.items():
        hexValues[name] = asc.hexlify(value)


   
'''(Nested list of patterns, by Chr) -> Nested list of color values, by Chr
takes the translated tree and convert each pattern to a value'''
def calculateColor(patternTree):
    colorTable = createColorTable(patternTree)
    newTree = {}
    
    for chrName, chr in list(patternTree.items()):
        newChr = []
        newTree[chrName] = newChr
        for pattern in chr:
            newChr.append((colorTable[pattern], next(iter(pattern))))
    
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
    itemName = "Cluster"
    
    with open(outfile, "wb") as output:
        #write the hit lit as the keys within the first block 
        output.write(bytes('$\n{0}\n$\n'.format(itemName)).encode('utf-8'))
        
        #write the rest
        for chrName, chr in resultsTree.items():
            output.write(bytes('%s\n' % chrName).encode("utf-8"))
            count = 0
            for block in chr:
                output.write(bytes('#%d\n%s - %s: %d\n'%(count, selfName, block[1], block[0])).encode("utf-8"))
                count+=1
                
if __name__ == '__main__':
#     directory = 'F:\Documents\ProjectData\\64Genomes\Counting' 
#     filepath = 'F:\Documents\ProjectData\\64Genomes\Counting\persistentResult.txt'
#     tabpath = 'F:\Documents\ProjectData\\64Genomes\Counting\persistentMatrix.tab'
#     outpath = 'F:\Documents\ProjectData\\64Genomes\Counting\pattern.txt'

    filepath = '/data/javi/Toxo/64Genomes/Counting/persistentResult.txt'
    tabpath = '/data/javi/Toxo/64Genomes/Counting/persistentMatrix.tab'
    outpath = '/data/javi/Toxo/64Genomes/Counting/pattern.txt'
#     strains = ['ME49', 'RAY', 'PRU', 'ARI', 'B73', 'B41']
#     strains = ['CAST', 'TgCkCr1', 'RH-88', 'RH-JSR', 'GT1', 'TgCkCr10', 'TgCkBr141']
    strains = ['BRC_TgH_18002_GUY-KOE', 'GUY-2004-ABE', 'BRC_TgH_18003_GUY-MAT', 'BRC_TgH_18009', 'BRC_TgH_18021', 'GUY-2003-MEL', 'BRC_TgH_18001_GUY-DOS']
    clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(filepath, tabpath)[0])
    results = calculateColor(translate(clusterTree, strains))
    write(results, outpath)
    print("end of ClusterPattern script")