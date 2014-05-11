'''
Created on Apr 24, 2014

@author: javi
'''
from decimal import *
import numpy as np

getcontext().prec = 5

'''(dictionary) -> dictionary
Figure out what the lower bounds should be, and how to properly space the edge
widths to best reflect the observable connections'''
def pruneMatrix(matrix):
    #construct the np matrix
    sampleList = sorted(matrix.keys())
    npMatrix = np.zeros((len(sampleList), len(sampleList)))
    for index, source in enumerate(sampleList):
            for target in sampleList:
                npMatrix[sampleList.index(source)][sampleList.index(target)] = matrix[source][target]
                
    #identify lower bounds
    dropCutoff = analyzeMatrix(npMatrix)
    print(dropCutoff)
    filteredMatrix = npMatrix * (npMatrix >= dropCutoff)
    transformedMatrix = transformMatrix(filteredMatrix)
    
    return transformedMatrix
    
'''helper for prune. Calculates a few statistics about the matrix
and returns them

dropCutoff: half way between the average and the highest'''
def analyzeMatrix(npMatrix):
    cutoffStringency = 0.65
    
    #current method is to go half way between max and median. 
    actualValues = npMatrix * (npMatrix < npMatrix[0,0])
    
    cutoffs = np.zeros(np.shape(actualValues)[0])
    for x in np.arange(cutoffs.shape[0]):
        column = actualValues[:,x]
        clusterIterator = iter(intervalCluster(column))
        belowCutoff = 0.
        cutoff = 0
        
        try:
            while belowCutoff/np.shape(column)[0] < cutoffStringency:
                currCluster = next(clusterIterator)[0]
                cutoff = max(currCluster)
                belowCutoff += len(currCluster)
        except StopIteration:
            pass
        
        cutoffs[x] = cutoff
#         another old algorithm
#         cutoffs[x] = np.max(column) + np.sort(column)[len(column)/2] / 2
#     old algorithm, discarded. 
#     dropCutoff = (npMatrix.max() + sorted(actualValues.flat)[len(actualValues.flat)/2]) / 2
    return cutoffs

'''(ndarray) -> list of lists
short clustering algorithm. logic: if not within x% of any "centers", current
value becomes a new center
results sorted lowest to highest, both amongst and within clusters'''
def intervalCluster(npMatrix):
    anchors = []
    results = {}
    closeness = 0.2
    sortedMatrix = np.sort(npMatrix.flat)[::-1]
    
    for e in sortedMatrix[np.where(sortedMatrix>0)]:
        for anchor in anchors:
            if abs(1 - anchor/e) < closeness * (np.log10(e)+1):
                break
        else:
            anchors.append(e)
            results[e] = []
    anchors = np.array(anchors)
    
    for e in sortedMatrix:
        if e > 0:
            distances = abs(anchors - e)
            closest = np.where((distances - np.min(distances)) == 0)[0][0]
            results[anchors[closest]].append(e)
    
    return sorted([[sorted(values)] for anchor, values in results.items()])
        
        
    
'''same idea, really'''
def transformMatrix(npMatrix):
    #values are mapped log2 from zero to maxValue
    maxValue = 2
    scalingFactor = 10
    
    max = npMatrix.max()
    toFraction = npMatrix / max
    
#     retired, trying something new
#     highAverage = np.average(toFraction, weights=npMatrix>0)
#     
#     transformedMatrix = (np.power(toFraction / highAverage, scalingFactor)) * maxValue
    
    transformedMatrix = np.power(toFraction, np.log(toFraction) * -1 * scalingFactor) * maxValue
    
    return np.triu(transformedMatrix)
    
    
    
    
        
'''(dictionary) -> String
takes a single matrix (2x nested dictionary) and converts it to a series of nodes and edges
in accordance with the GML format, to be incorportated into a larger
GML file for cytoscape.'''
def fromMatrix(matrix, name):
    print(name + " cutoff is:")
    sampleList = sorted(matrix.keys())
    #nodes are two-tuples consisting of id and label
    nodes = [(index, name) for index, name in enumerate(sampleList)]
    edges = []
    prunedMatrix = pruneMatrix(matrix)
    
    #edges consists of three-tuples 
    for sourceIndex, source in enumerate(sampleList):
        for targetIndex, target in enumerate(sampleList):
            value = prunedMatrix[sourceIndex, targetIndex]
            if value > 0:
                edges.append((source, target, value))
    
    text = "\
    graph [\n\
{0}\n\
{1}\n\
{2}\n\
    ]\n".format(getHeader(name), "\n".join([getNodeText(node[0], node[1]) for node in nodes]), "\n".join([getEdgeText(sampleList.index(edge[0]), sampleList.index(edge[1]), edge[2]) for edge in edges]))
    return text



    
'''(Dictionary) -> None (write to file)
primary function of this encoder, to be called by the outside source.
Currently only accepts dataTree style input. Each parse call creates one file. 
'''
def parse(matrix, name, sampleList, outfile):
    with open(outfile, 'w') as output:
        output.write(fromMatrix(matrix, name))

'''for the sake of readability, the header info for the GML file will be
stored here'''
def getHeader(label):
    header = "\
        directed    0\n\
        label    \"{0}\"\n".format(label)
    return header

'''(Int, String) -> String
same idea for the node text.'''
def getNodeText(ID, label):
    text = "\
        node [\n\
            id    {0}\n\
            label    \"{1}\"\n\
            graphics [\n\
                w    1\n\
                h    1\n\
            ]\n\
        ]\n".format(ID, label)
    return text

'''(String, String, int) -> String
same.'''
def getEdgeText(source, target, width):
   
    if not width > 0 or source == target:
        return ""
    
    text = "\
        edge [\n\
            source    {0}\n\
            target    {1}\n\
            graphics [\n\
                width    {2}\n\
                type    \"line\"\n\
            ]\n\
        ]\n".format(source, target, width)
    return text

# '''(int) -> decimal
# given the raw # of connections, calculates an appropriate edge width. 
# setting inside'''
# def calibrateWidth(width):
#     min_width = 120
#     scaling_factor = Decimal(300)
#     
#     return max(0, Decimal(width - min_width) / scaling_factor)
    
if __name__ == '__main__':
    pass