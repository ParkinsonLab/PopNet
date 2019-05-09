'''
Created on Apr 24, 2014

@author: javi
'''
from decimal import *
import numpy as np
import json

getcontext().prec = 5

'''(Dictionary) -> None (write to file)
primary function of this encoder, to be called by the outside source.
Currently only accepts dataTree style input. Each parse call creates one file. 
'''
def parse(matrix, color_table, composition, group_names, overall_clusters, sample_list, prefix):
    
    matrix = matrix.values / np.max(matrix.values)
    json_name = prefix + '.json'
    xgmml_name = prefix + '.xgmml'
    
    with open(json_name, 'w') as output:
        json.dump(jsonFromMatrix(matrix, color_table, composition, group_names, overall_clusters, sample_list), output)

    with open(xgmml_name, 'w') as output:
        output.write(xgmmlFromMatrix(matrix, color_table, composition, group_names, overall_clusters, sample_list, prefix))    
'''(dictionary) -> String
takes a single matrix (2x nested dictionary) and converts it to a series of nodes and edges
in accordance with the GML format, to be incorportated into a larger
GML file for cytoscape.'''
def jsonFromMatrix(matrix, color_table, composition, group_names, overall_clusters, sample_list):
    def getGroup(name):
        for i, g in enumerate(overall_clusters):
            if name in g:
                return group_names[i]
        raise ValueError('json from matrix got a name, {0}, not in group names'.format(name))
    
    group_names.append('None')
    nodes = [getNodeJson(idx, node, getGroup(node), composition[idx], group_names) for idx, node in enumerate(sample_list)]
    edges = []
    #edges consists of three-tuples and excludes duplicated
    for source in range(len(sample_list)):
        for target in range(source, len(sample_list)):
            if not source == target:
                value = matrix[source, target]
                if value > 0:
                    edge = getEdgeJson(sample_list[source], sample_list[target], value)
                    edges.append(edge)    
    
    
    sample_colors = {x: toHexColor(color_table['color'].loc[getGroup(x)]) for x in sample_list}
    group_colors = {x: toHexColor(color_table['color'].loc[x]) for x in group_names}
    colors = {}
    colors.update(sample_colors)
    colors.update(group_colors)

    text = {'names':sample_list, 
            'nodes':nodes,
            'edges':edges, 
            'colorTable': colors
            }

    return text

'''(dictionary) -> String
takes a single matrix (2x nested dictionary) and converts it to a series of nodes and edges
in accordance with the GML format, to be incorportated into a larger
GML file for cytoscape.'''
def xgmmlFromMatrix(matrix, color_table, composition, group_names, overall_clusters, sample_list, name):
    def getGroup(name):
        for i, g in enumerate(overall_clusters):
            if name in g:
                return group_names[i]
        raise ValueError('xgmml from matrix got a name, {0} not in group names'.format(name))

    nodes = [getNodeText(idx, node, color_table, composition[idx], getGroup(node)) for idx, node in enumerate(sample_list)]
    edges = []

    for source in range(len(sample_list)):
        for target in range(source, len(sample_list)):
            if not source == target:
                value = matrix[source, target]
                if value > 0:
                    edge = getEdgeText(source, target, value)
                    edges.append(edge)  
    
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
    </graph>".format(name, "\n".join(nodes), \
"\n".join(edges))

    return text


'''(Int, String) -> String
same idea for the node text.'''
def getNodeJson(ID, label, group, composition, group_names):      
    lengths = []
    ids = []
    for section in composition:
        lengths.append(str(section[0]))
        ids.append(group_names[section[1]])
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
    def clamp(x): 
        return int(max(0, min(x, 255)))
    r, g, b = num
    return "#{0:02x}{1:02x}{2:02x}".format(clamp(r * 255), clamp(g * 255), clamp(b * 255))

'''(String, String, int) -> String
same.'''
def getEdgeJson(source, target, width):
    return {
                'source': source,
                'target': target,
                'width': width
                }


def getNodeText(ID, label, color_table, composition, group):
    try:
        color = toHexColor(color_table['color'][group])
    except KeyError:
        print("KeyError!")
        color = "#000000"
        
    circosText = calculateCircos(composition, color_table)
#     #hacked for importing!
#     circosText = importCircos(composition)
    text = "\
        <node label=\"{1}\" id=\"{0}\" >\n\
            <att name=\"name\" type=\"string\" value=\"{1}\"/>\n\
            <att name=\"group\" type=\"string\" value=\"{2}\"/>\n\
            <att name=\"color\" type=\"string\" value=\"{6}\"/>\n\
            <att name=\"Gradient\" type=\"string\" value=\"circoschart: arcstart=90 firstarc=.6 arcwidth=.3 outlineWidth=0.0 colorlist=&quot;{3}&quot; showlabels=false  attributelist=&quot;{4}&quot;\"/>\n\
{5}\
            <graphics h=\"60\" w=\"60\" outline=\"{2}\" />\n\
        </node>\n".format(ID, label, group, circosText[0], "values", circosText[1], color)
    return text

def getEdgeText(source, target, width):

    if not width > 0 or source == target:
        return ""
        
    text = "\
        <edge source=\"{0}\" target=\"{1}\" >\n\
            <graphics width=\"{2}\" />\n\
            <att name=\"width\" type=\"float\" value=\"{2}\" />\n\
        </edge>\n".format(source, target, width)
    return text

'''list, dict -> tuple
creates the two parameter lists, colorList and valueList, for 
a node's circos chart'''
def calculateCircos(composition, color_table):
    valueList = []
    colorList = []
    for section in reversed(composition):
        valueList.append(str(section[0]))
        colorList.append(str(toHexColor(color_table['color'].loc[color_table.index[section[1]]])))

    colorString = ",".join(colorList)
    valueString = "\
            <att name=\"values\" type=\"list\">\n"
    for item in valueList:
        valueString += "\
                <att type=\"real\" value=\"{0}\"/>\n".format(item)
    valueString += "\
            </att>\n"
    return colorString, valueString
        
if __name__ == '__main__':
    pass