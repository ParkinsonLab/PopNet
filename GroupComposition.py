'''
Created on Apr 28, 2014

@author: javi
'''
'''in fulfillment of John's request for a chromosome coloring scheme where each section
is colored according to which group is based on.. me looking at it. zzz.
Use the createColorTable and calculateColor from ClusterPattern to visualize.'''

'''Strain -> Number
number represents the group ID, as specified in this function'''
import numpy as np
def getGroups(strain):
    
    groups = {}
    groups['GUYS'] = ["BRC_TgH_18003_GUY-MAT", "BRC_TgH_18001_GUY-DOS", "BRC_TgH_18021", \
                  "BRC_TgH_18009", "BRC_TgH_18002_GUY-KOE", "GUY-2003-MEL", "GUY-2004-ABE", \
                  "RUB", "VAND"]
    groups['TgH'] = ["BRC_TgH_20005", "CASTELLS", "TgH26044", "BRC_TgH_21016"]
    groups['ME49'] = ["COUG", "GUY-2004-JAG1", "TgCat_PRC2", "ARI", "RAY", "PRU", \
                   "B41", "ME49", "B73"]
    groups['VEG'] = ["TgShUS28", "TgCkGy2", "ROD", "M7741", "VEG", "SOU", "G662M"]
    groups['p89'] = ["TgCatBr64", "p89", "TgCatBr3", "TgCatBr15"]
    groups['TgCats'] = ["TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
                   "TgCatBr44", "MAS", "TgCatBr1", "TgCatBr34"]
    groups['GAL-DOM1/2'] = ["CAST", "FOU", "GAB5-2007-GAL-DOM1", \
                   "GAB3-2007-GAL-DOM2", "BOF"]
    groups['GAL-DOM10'] = ["GAB2-2007-GAL-DOM2", "GAB5-2007-GAL-DOM6", "GAB1-2007-GAL-DOM10", "GAB3-2007-GAL-DOM9"]
    groups['GT1'] = ["TgDogCo17", "GT1", "RH-JSR", "RH-88", "TgCkCr1", "TgCkBr141"]
    notUsed  = ["TgRsCr1", "TgCtCo5", "TgCkCr10", "TgCatBr72", "TgCATBr9", "TgCatBr26" ]
    for name, group in groups.items():
        if strain in group:
            del groups[name]
    
    return groups

    
    '''ORIGINAL GROUPS
    groups['GUYS'] = ["BRC_TgH_18003_GUY-MAT", "BRC_TgH_18001_GUY-DOS", "BRC_TgH_18021", \
                  "BRC_TgH_18009", "BRC_TgH_18002_GUY-KOE", "GUY-2003-MEL", "GUY-2004-ABE", \
                  "RUB", "VAND"]
    groups['TgH'] = ["BRC_TgH_20005", "CASTELLS", "TgH26044", "BRC_TgH_21016"]
    groups['ME49'] = ["COUG", "GUY-2004-JAG1", "TgCat_PRC2", "ARI", "RAY", "PRU", \
                   "B41", "ME49", "B73", "SOU"]
    groups['VEG'] = ["TgShUS28", "TgCkGy2", "ROD", "M7741", "VEG"]
    groups['p89'] = ["G662M", "p89", "TgCatBr3", "TgCatBr15"]
    groups['TgCats'] = ["TgCatBr64", "TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
                   "TgCatBr44", "MAS", "TgCatBr1"]
    groups['GAL-DOM'] = ["TgCatBr72", "TgCATBr9", "TgCatBr26", "FOU", "GAB5-2007-GAL-DOM1", \
                   "GAB3-2007-GAL-DOM2", "GAB3-2007-GAL-DOM9", "GAB1-2007-GAL-DOM10", \
                   "GAB5-2007-GAL-DOM6", "GAB2-2007-GAL-DOM2", "BOF"]
    groups['GT1'] = ["CAST", "TgDogCo17", "GT1", "RH-JSR", "RH-88", "TgCkCr1", "TgCkBr141", \
                   "TgRsCr1", "TgCtCo5", "TgCkCr10"]
                   '''

'''(string, 2-Nested list, string) -> int
returns the index, in the first dimension, of the group that this strain is most
closely related to, in this particular block'''    
def matchToGroup(strain, groups, block):
    counts = []
    for name, group in groups.items():
        counts.append((countMatches(strain, group, block), name))
    return sorted(counts, reverse=True)[0][1]

'''(string, list, block) -> int
counts how many of the strains in the group clusters together with
the given strain'''
def countMatches(strain, group, block):
    for line in block:
        if strain in line:
            return len(set(group).intersection(line)) / float(len(group))
    return 0

def findComposition(strain, dataTree):
    results = {}
    groups = getGroups(strain)
    for name, chr in dataTree.items():
        results[name] = chrResults = []
        for block in chr:
            chrResults.append(frozenset([matchToGroup(strain, groups, block)]))
    return results

'''([Strains], dataTree) -> output
expanding on the single strain comparison, what if I want multiple strains from the same
group?'''
def multiComposition(strains, dataTree):
    from collections import Counter as counter
    resultsList = [findComposition(strain, dataTree) for strain in strains]
    compiledResults = {}
    for chr in sorted(resultsList[0].keys()):
        chrMatrix = np.array([strain[chr] for strain in resultsList])
        compiledMatrix = [] * chrMatrix.shape[1]
        for x in range(chrMatrix.shape[1]):
            #this line sorts the counter dictionary, and takes the first tuple (most common)'s first value (the color)
            compiledMatrix.append(sorted(counter(chrMatrix[:, x]).items(), key=lambda x: x[1], reverse=True)[0][0])
        compiledResults[chr] = compiledMatrix
    return compiledResults
    
if __name__ == '__main__':
    pass