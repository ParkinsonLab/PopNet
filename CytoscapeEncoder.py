'''
Created on Apr 24, 2014

@author: javi
'''
from decimal import *
getcontext().prec = 5
'''(dictionary) -> String
takes a single matrix (2x nested dictionary) and converts it to a series of nodes and edges
in accordance with the GML format, to be incorportated into a larger
GML file for cytoscape.'''
def fromMatrix(matrix, name, sampleList):
    #nodes are two-tuples consisting of id and label
    nodes = [(index, name) for index, name in enumerate(sampleList)]
    edges = []
    
    #edges consists of three-tuples 
    for index, source in enumerate(sampleList):
        for target in sampleList[index:]:
            edges.append((source, target, matrix[source][target]))
    
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
        output.write(fromMatrix(matrix, name, sampleList))

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
                w    0.2\n\
                h    0.2\n\
            ]\n\
        ]\n".format(ID, label)
    return text

'''(String, String, int) -> String
same.'''
def getEdgeText(source, target, width):
    calibrated_width = calibrateWidth(width)
    
    if not calibrated_width > 0:
        return ""
    
    text = "\
        edge [\n\
            source    {0}\n\
            target    {1}\n\
            graphics [\n\
                width    {2}\n\
                type    \"line\"\n\
            ]\n\
        ]\n".format(source, target, calibrated_width)
    return text

'''(int) -> decimal
given the raw # of connections, calculates an appropriate edge width. 
setting inside'''
def calibrateWidth(width):
    min_width = 120
    scaling_factor = Decimal(300)
    
    return max(0, Decimal(width - min_width) / scaling_factor)
    
if __name__ == '__main__':
    pass