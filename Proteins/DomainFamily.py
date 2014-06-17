'''
Created on Jun 13, 2014

@author: javi
'''

'''contains routines for choosing the best family from HMM searches using different
profiles. '''
import re


'''tree, tree -> tree
merges the information in target and source, which represent
search result using two profiles. Picks the best family of the two'''
def update(new, current):
    for proteinID, proteinInfo in new.items():
        if proteinID in current:
            currProtein = current[proteinID]
            for domain, score in proteinInfo.items():
                if domain not in currProtein or score[1] < currProtein[domain][1]:
                        current[proteinID][domain] = new[proteinID][domain]
        else:
            current[proteinID] = proteinInfo
    return current
        
'''list of trees -> tree

refer to HMMParser for the tree structure. 
Main running method for choosing the best family'''
def choose(treeList):
    base = treeList[0]
    others = treeList[1:]
    for other in others:
        base = update(other, base)
    return base

'''tree -> None(print)

prints a tree to file in a legible format'''
def printTree(tree, filepath):
    with open(filepath, "w") as output:
        count = 0
        for proteinID, proteinInfo in tree.items():
            output.write(">>{0}\n".format(proteinID))
            count += 1
            for domain, info in proteinInfo.items():
                family, score, length, cys, degen = info
                output.write("\tDomain {0:>3} Family {1} Length {3:>3} Cys {4:>2} Degen {5:>5} Score {2}\n".format(domain, family, score, length, cys, str(degen)))
        print("{0} has {1} proteins".format(re.match("^.+/(.+?)$", filepath).group(1), count))
            
# '''tree -> tree
# 
# prunes the current tree according to rules. Current rules are two of the three from the paper:
# 
# e-value < 1e-5
# at least 90 in length
# '''
# def prune(tree):
#     for proteinID, proteinInfo in tree.items():
#         for domainID, domainInfo in proteinInfo.items():
#             score = domainInfo[1]
#             length = domainInfo[2]
#             if score > 0.00001 or length < 90:
#                 del tree[proteinID][domainID]

