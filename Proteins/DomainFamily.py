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
            for domain, score in proteinInfo.items():
                current[proteinID] = updateDomains(score, current[proteinID])
        else:
            current[proteinID] = proteinInfo
    return current

'''tuple, dict -> dict
primary method for updating a domain. runs the other parts'''
def updateDomains(newDomain, currentProtein):
    decision = choosePosition(newDomain, currentProtein)
    if decision[1] == True:
        score = newDomain[1]
        oldScore = currentProtein[decision[0]][1]
        if score < oldScore: #smaller score is better
            currentProtein[decision[0]] = newDomain
    else:
        currentProtein = place(newDomain, currentProtein, decision[0])
    return currentProtein
        
    
'''tuple, dict, int -> dict
places the given domain at the appropriate index, and update all other indices'''
def place(newDomain, oldDomains, index):
    for domain in sorted(oldDomains.keys()[index-1:], reverse=True):
        oldDomains[domain+1] = oldDomains[domain]
    oldDomains[index] = newDomain
    return oldDomains

'''int, list of ints -> int
returns a or b, whichever one is closer to the query'''
def closest(query, list):
    results = [(abs(query-x), index) for index, x in enumerate(list)]
    return list[sorted(results)[0][1]]

'''tuple, dict -> (int, bool)
returns the index of where it should go, and whether its a replace or insert (True for replace)'''
def choosePosition(newDomain, currentProtein):
    oldCoords = [x[1][2] for x in sorted(currentProtein.items())]
    newCoords = newDomain[2]
    
    starts = [x[0] for x in oldCoords]
    ends = [x[1] for x in oldCoords]
    
    startClosest = closest(newCoords[0], starts+ends)
    endClosest = closest(newCoords[1], starts+ends)
    
    #prediction based on start
    if startClosest in starts:
        sPrediction = (starts.index(startClosest)+1, True)
    else:
        sPrediction = (ends.index(startClosest)+2, False)
            
    #prediction based on end
    if endClosest in ends:
        ePrediction = (ends.index(endClosest)+1, True)
    else:
        ePrediction = (starts.index(endClosest)+1, False)

    #pick one
    
    #ends
    if startClosest == endClosest:
        if startClosest in starts:
            return ePrediction
        elif endClosest in ends:
            return sPrediction
    #both predictions would be wrong
    elif newCoords[0] < startClosest and newCoords[1] > endClosest:
        return (sPrediction[0], True)
    elif sPrediction == ePrediction:
        return sPrediction
    elif startClosest in starts:
        return sPrediction
    elif endClosest in ends:
        return ePrediction
    else:
        raise Exception('unable to choose')


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
                output.write("\tDomain {0:>3} Family {1} Coords {3:>7} Cys {4:>2} Degen {5:>5} Score {2}\n".format(domain, family, score, "-".join([str(length[0]),str(length[1])]), cys, str(degen)))
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

'''trees -> None (output)

computes some stats about this tree'''
def summary(tree):
    pass