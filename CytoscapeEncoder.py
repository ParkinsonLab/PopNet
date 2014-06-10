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
def fromMatrix(matrix, name, colorTable):
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
<?xml version=\"1.0\"?>\n\
    <graph label=\"{0}\"\n\
    xmlns:dc=\"http://purl.org/dc/elements/1.1/\"\n\
    xmlns:xlink=\"http://www.w3.org/1999/xlink\"\n\
    xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"\n\
    xmlns:cy=\"http://www.cytoscape.org\"\n\
    xmlns=\"http://www.cs.rpi.edu/XGMML\"\n\
    directed=\"0\" >\n\
{1}\n\
{2}\n\
    </graph>".format(name, "\n".join([getNodeText(node[0], node[1], colorTable) for node in nodes]), "\n".join([getEdgeText(sampleList.index(edge[0]), sampleList.index(edge[1]), edge[2], colorTable) for edge in edges]))
    return text


    
'''(Dictionary) -> None (write to file)
primary function of this encoder, to be called by the outside source.
Currently only accepts dataTree style input. Each parse call creates one file. 
'''
def parse(matrix, name, sampleList, outfile, colorTable):
    with open(outfile, 'w') as output:
        output.write(fromMatrix(matrix, name, colorTable))

# '''for the sake of readability, the header info for the GML file will be
# stored here'''
# def getHeader(label):
#     header = "".format(label)
#     return header

'''(Int, String) -> String
same idea for the node text.'''
def getNodeText(ID, label, colorTable):
    try:
        color = toHexColor(colorTable[label])
    except KeyError:
        color = "000000"
    
    text = "\
        <node label=\"{1}\" id=\"{0}\" >\n\
            <att name=\"name\" type=\"string\" value=\"{1}\"/>\n\
            <att name=\"group\" type=\"string\" value=\"{2}\"/>\n\
            <graphics h=\"10\" w=\"10\" fill=\"{2}\" />\n\
        </node>\n".format(ID, label, color)
    return text

'''(Int) -> hex
helper function for getNodeText, mostly. May have other uses'''
def toHexColor(num):
    binaryStr = "{:0>24b}".format(num)
    RED = int(binaryStr[0:8], base=2)
    GREEN = int(binaryStr[8:16], base=2)
    BLUE = int(binaryStr[16:24], base=2)
    
    pattern = "{:0>2X}"*3
    result = pattern.format(RED, GREEN, BLUE)
    return result

'''(String, String, int) -> String
same.'''
def getEdgeText(source, target, width, colorTable):
   
    if not width > 0 or source == target:
        return ""
    
    try:
        color = toHexColor(colorTable[target])
    except KeyError:
        color = "000000"
    
    text = "\
        <edge source=\"{0}\" target=\"{1}\" >\n\
            <graphics width=\"{2}\" fill=\"{3}\" />\n\
        </edge>\n".format(source, target, width, color)
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