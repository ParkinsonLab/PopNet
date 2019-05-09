'''
Created on Aug 6, 2014

@author: javi

tries to find the parents of a strain
'''
import GroupComposition as gc


def findParents(strains, clusterMatrix, aggregateCountMatrix):
    
    
    #construct shortlist
    shortList = {}
    shortListLength = 5
    for strain in strains:
        shortList[strain] = []
    for strain, counts in aggregateCountMatrix.values()[0].items():
        if strain in strains:
            allowedStrainGroups = gc.getGroups(strain)
            allowedStrains = []
            for group, strainlist in allowedStrainGroups.items():
                for x in strainlist:
                    allowedStrains.append(x)
            countsList = []
            for strain2, count in counts.items():
                if strain2 is not strain and strain2 in allowedStrains:
                    countsList.append((count, strain2))
            shortList[strain] = [x[1] for x in sorted(countsList, reverse=True)][:shortListLength]
    
    #construct pair lists for each strain
    pairListMatrix = {}
    for strain, list in shortList.items():
        pairs = []
        pairListMatrix[strain] = {}
        for index, item1 in enumerate(list):
            for item2 in list[index:]:
                pairs.append((item1, item2))
        for pair in pairs:
            pairListMatrix[strain][pair] = [0,0,0] #edit this if you are going to use triplets and so on. 
    
    #check each pair
    totalLen = 0.
    for chrName, chr in clusterMatrix.items():
        totalLen += len(chr)
        for section in chr:
            for line in section:
                    for strain, pairs in pairListMatrix.items():
                        if strain in line:
                            for pair in pairs.keys():
                                for index, strain2 in enumerate(pair):
                                    if strain2 in line:
                                        pairListMatrix[strain][pair][2] += 1
                                        pairListMatrix[strain][pair][index] += 1
                                        break
    
    #return only top three for each strain:
    topLength = 3
    for strain, counts in pairListMatrix.items():
        topPicks = sorted([([z / totalLen for z in y], x) for x, y in counts.items()], reverse=True, key=lambda x:x[0][2])
        pairListMatrix[strain] = topPicks[:topLength]
       
    return pairListMatrix
                                
def printParents(pairListMatrix, outpath):
    with open(outpath, "w") as output:
        for strain, parentList in pairListMatrix.items():
            line = []
            line.append(strain)
            for count, pair in parentList:
                pattern = ",{:>20}" * len(pair) + ",{:>5}"*3
                args = list(pair) + count
                line.append(pattern.format(*args))
            line.append("")
            output.write("\n".join(line))
    
    
        
                    