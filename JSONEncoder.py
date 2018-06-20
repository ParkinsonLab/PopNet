'''
Created on Apr 24, 2014

@author: javi
'''
from decimal import *
import numpy as np
from CytoscapeEncoder import toHexColor

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
    
    #the axis needs to be changed for proper comparison against the full matrix!
    tempx = dropCutoff.shape[0]
    dropCutoff.shape = (tempx, 1)
    
    filteredMatrix = npMatrix * (npMatrix >= dropCutoff)
    np.fill_diagonal(filteredMatrix, 0)
# we're not filtering anymore
#     transformedMatrix = transformMatrix(filteredMatrix)
    transformedMatrix = transformMatrix(npMatrix)
    
#     #debug code delete later
#     key = 'EC9-8_S288C'
#     index = sampleList.index(key)
#     
#     print('{0} has cutoff off {1}'.format(key, str(dropCutoff[index])))
#     print('with the values \n' + str(npMatrix[index]))
#     print('result \n' + str(filteredMatrix[index]))
#     print('and the final one \n' + str(transformedMatrix[index]))
    
    return transformedMatrix
    
'''helper for prune. Calculates a few statistics about the matrix
and returns them

dropCutoff: half way between the average and the highest'''
def analyzeMatrix(npMatrix):
    #This value determines how many edges you get.
    cutoffStringency = 0.8
    
    #current method: cutoff is calculated from the interval cluster
    #then anything less than the cutoff is pruned. 
    actualValues = npMatrix * (npMatrix < npMatrix[0,0])

    cutoffs = np.zeros(np.shape(actualValues)[0])
    for x in np.arange(cutoffs.shape[0]):
        column = actualValues[:,x]
        clusterIterator = iter(intervalCluster(column))
        belowCutoff = 0.
        cutoff = 0
        previous = 0
        try:
            while belowCutoff/np.shape(column)[0] < cutoffStringency:
                currCluster = next(clusterIterator)[0]
                previous = cutoff
                cutoff = max(currCluster)
                belowCutoff += len(currCluster)
                cutoffs[x] = cutoff
        except StopIteration:
            cutoffs[x] = previous
        
        
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
    max = np.max(npMatrix.flat)
    anchors = []
    results = {}
    closeness = 0.2
    sortedMatrix = np.sort(npMatrix.flat)[::-1]
    
    for e in sortedMatrix[np.where(sortedMatrix>0)]:
        for anchor in anchors:
            #about this formula:
            #it basically measures if the value is [closeness] close to an anchor, to decide
            #if it should be a new anchor. The allowable range for anchors is increased for higher ranges,
            #hence the log part, but should always be kept to less than 30% to avoid poor clustering at
            #higher values. 
            if abs(1 - anchor/e) < closeness * (1+1/(1+np.exp(-1/(max-e)))):
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
    maxValue = 4.
    scalingFactor = 10.
    
    max = npMatrix.max()
    toFraction = npMatrix / max
    
#     retired, trying something new
#     highAverage = np.average(toFraction, weights=npMatrix>0)
#     
#     transformedMatrix = (np.power(toFraction / highAverage, scalingFactor)) * maxValue
    
#     transformedMatrix = np.power(toFraction, np.log(toFraction) * -1 * scalingFactor) * maxValue
    transformedMatrix = toFraction
    
    
    return transformedMatrix

        
'''(dictionary) -> String
takes a single matrix (2x nested dictionary) and converts it to a series of nodes and edges
in accordance with the GML format, to be incorportated into a larger
GML file for cytoscape.'''
def fromMatrix(matrix, name, colorTable, composition, groups):
    print(name + " cutoff is:")
    sampleList = sorted(matrix.keys())
    #nodes are two-tuples consisting of id and label
    nodes = [(index, nodeName) for index, nodeName in enumerate(sampleList)]
    edges = []
    prunedMatrix = pruneMatrix(matrix)
    
    #edges consists of three-tuples and excludes duplicated
    excluded = []
    for source in sorted(nodes):
        for target in sorted(nodes):
            if (source, target) not in excluded and not source == target:
                value = prunedMatrix[source[0], target[0]]
                if value > 0:
                    edge = getEdgeText(source, target, value)
                    edges.append(edge)    
                    excluded.append((target, source))
    
    text = {'names':sampleList, 
            'nodes':[getNodeText(node[0], node[1], composition[node[1]], groups[node[1]]) for node in nodes],
            'edges':edges, 
            'colorTable':{x: toHexColor(colorTable[x]) for x in colorTable}
            }
    return text


    
'''(Dictionary) -> None (write to file)
primary function of this encoder, to be called by the outside source.
Currently only accepts dataTree style input. Each parse call creates one file. 
'''
def parse(matrix, name, sampleList, outfile, colorTable, composition, groups):
    import json
    with open(outfile, 'w') as output:
        json.dump(fromMatrix(matrix, name, colorTable, composition, groups), output)

# '''for the sake of readability, the header info for the GML file will be
# stored here'''
# def getHeader(label):
#     header = "".format(label)
#     return header

'''(Int, String) -> String
same idea for the node text.'''
def getNodeText(ID, label, composition, group):      
    lengths = []
    ids = []
    for section in composition:
        lengths.append(str(section[1] - section[0] + 1))
        ids.append(str(section[2]))
#     #hacked for importing!
#     circosText = importCircos(composition)
    return {
        'name': label,
        'group': group,
        'id': ID,
        'ids': ids,
        'lengths': lengths
        }

'''(Int) -> hex
helper function for getNodeText, mostly. May have other uses'''
def toHexColor(num):
    binaryStr = "{:0>24b}".format(num)
    RED = int(binaryStr[0:8], base=2)
    GREEN = int(binaryStr[8:16], base=2)
    BLUE = int(binaryStr[16:24], base=2)
    
    pattern = "{:0>2X}"*3
    result = pattern.format(RED, GREEN, BLUE)
    return "#" + result

'''(String, String, int) -> String
same.'''
def getEdgeText(source, target, width):
    return {
                'source': source[1],
                'target': target[1],
                'width': width
                }

'''list -> tuple
takes a precalculated list of colors and puts them into the circos format'''
def importCircos(composition):
    valueList = []
    colorList = []
    for section in composition:
        valueList.append(str(section[1] - section[0] + 1))
        colorList.append(str(section[2]))
    
    colorString = ",".join(colorList)
    valueString = "\
            <att name=\"values\" type=\"list\">\n"
    for item in valueList:
        valueString += "\
                <att type=\"real\" value=\"{0}\"/>\n".format(item)
    valueString += "\
            </att>\n"
            
    return colorString, valueString
    
    
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