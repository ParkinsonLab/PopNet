'''
Created on Apr 28, 2014

@author: javi
'''
'''in fulfillment of John's request for a chromosome coloring scheme where each section
is colored according to which group is based on.. me looking at it. zzz.
Use the createColorTable and calculateColor from ClusterPattern to visualize.'''

'''Strain -> Number
number represents the group ID, as specified in this function'''
def getGroups():
    
    groups = []
    groups.append(["BRC_TgH_18003_GUY-MAT", "BRC_TgH_18001_GUY-DOS", "BRC_TgH_18021", \
                  "BRC_TgH_18009", "BRC_TgH_18002_GUY-KOE", "GUY-2003-MEL", "GUY-2004-ABE", \
                  "RUB", "VAND"])
    groups.append(["BRC_TgH_20005", "CASTELLS", "TgH26044", "BRC_TgH_21016"])
    groups.append(["COUG", "GUY-2004-JAG1", "TgCat_PRC2", "ARI", "RAY", "PRU", \
                   "B41", "ME49", "B73", "SOU"])
    groups.append(["TgShUS28", "TgCkGy2", "ROD", "M7741", "VEG"])
    groups.append(["G662M", "p89", "TgCatBr3", "TgCatBr15"])
    groups.append(["TgCatBr64", "TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
                   "TgCatBr44", "MAS", "TgCatBr1"])
    groups.append(["TgCatBr72", "TgCATBr9", "TgCatBr26", "FOU", "GAB5-2007-GAL-DOM1", \
                   "GAB3-2007-GAL-DOM2", "GAB3-2007-GAL-DOM9", "GAB1-2007-GAL-DOM10", \
                   "GAB5-2007-GAL-DOM6", "GAB2-2007-GAL-DOM2", "BOF"])
    groups.append(["CAST", "TgDogCo17", "GT1", "RH-JSR", "RH-88", "TgCkCr1", "TgCkBr141", \
                   "TgRsCr1", "TgCtCo5", "TgCkCr10"])

    return groups

def getGroupNames():
    names = ['GUYs', 'TgH', 'ME49', 'VEG', 'p89', 'TgCats', 'GAL-DOM', 'GT1']
    return names
    
    
'''(string, 2-Nested list, string) -> int
returns the index, in the first dimension, of the group that this strain is most
closely related to, in this particular block'''    
def matchToGroup(strain, groups, names, block):
    counts = []
    for index, group in enumerate(groups):
        counts.append((countMatches(strain, group, block), names[index]))
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
    groups = getGroups()
    names = getGroupNames()
    for name, chr in dataTree.items():
        results[name] = chrResults = []
        for block in chr:
            chrResults.append(frozenset([matchToGroup(strain, groups, names, block)]))
    return results
        
    
if __name__ == '__main__':
    pass