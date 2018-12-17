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
        
'''int -> [int]
returns a list of colors in intrgb (like before) given a number'''
def generateColors(numColors):    
    #l and hInterval actually represents the number of intervals present, not the length of each.
    
    import colorsys
    max_h_interval = 8
    
    if numColors >= max_h_interval: hInterval = max_h_interval
    else: hInterval = numColors
    
    lInterval = (numColors-1)//max_h_interval + 1 
    
    colorList = []
    for x in range(1,lInterval+1):
        x = float(x)
        lum = x / (lInterval + 1)
        for y in range(1,hInterval+1):
            y = float(y)
            hue = (y + x%2/2) / hInterval
            print("Hue:{} Lum:{}".format(hue, lum))
            colorList.append(coordToInt(colorsys.hls_to_rgb(hue, lum, 1)))
    
    return iter(colorList[:numColors])
'''(f, f, f) -> int
helper for generateColors. Takes the rgb coords and return a int'''
def coordToInt(coord):    
    r = int(coord[0] * 255)
    g = int(coord[1] * 255)
    b = int(coord[2] * 255)
    
    hex = "0x{:02X}{:02X}{:02X}".format(r, g, b)
    print(hex)
    return int(hex, 0)

def createColorTable(patternTree, groups):
    import random
    
#     #quick fix to get the same colors; for 64 genomes
#    colorIter = iter([16737894, 1337452, 10040217, 6063411, 6710937, 6037279, 13421823, 16776960, 50637, 9737364, 16777215])
#     groups = ['GUYS', 'TgH', 'ME49', 'VEG', 'p89', 'TgCats', 'GAL-DOM1/2', 'GAL-DOM10', 'GT1', 'MISC', 'NONE']

    #for 8 genomes
#    colorIter = iter([255, 65280, 10066278, 16711680, 16777215])
#    groups = ['ME49', 'VEG', 'TYPEX', 'GT1', 'NONE']

    #yeast
#    colorIter = iter([255, 65280, 16711680, 16777215])
#    groups = ['S288C', 'AWR', 'UC5', 'NONE']
       
    results = {}
    colorIter = generateColors(len(groups))
    for group in groups:
        results[frozenset([group])] = next(colorIter)
    
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
            
#     printColorTable(results)
    return results

'''table is a dict of frozenset:value
just a helper function because i don't like their print'''
def printColorTable(table):
    maxKeyLength = max([len(x) for x in list(table.keys())])
    maxValueLength = max([len(str(x)) for x in list(table.values())])
#     print(maxKeyLength, maxValueLength)
    rowFormat = '{:<{klen}}\t{:<{vlen}}'
    result = "\n".join([rowFormat.format(x[0], hex(x[1]), klen=maxKeyLength, vlen=maxValueLength) for x in list(table.items())])
    return result
    
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
                output.write(bytes('#%d\n%s - %s: %d\n'%(count, selfName, itemName, block)).encode("utf-8"))
                count+=1


    
                
if __name__ == '__main__':
#Testing Only
    print(generateColors(10))