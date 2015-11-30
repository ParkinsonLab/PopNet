'''
Created on Jun 19, 2014

@author: javi
'''
'''Similar to cytoscape encoder, this file encodes the results of the SRS protein searches into
cytoscape format for ease of viewing'''
import os
import re 

'''String -> list
loads clusters into 2D list'''
def loadClusters(filepath):
    with open(filepath) as input:
        data = input.read()
    results = []
    
    lines = re.split("\n", data)[:-1]
    for line in lines:
        elements = line.split(" ")[1:]
        fixedelements = [re.split("[|]", element)[1] for element in elements]
        results.append(fixedelements)
    return results
        


'''String -> String
gets the color for that family, for the purpose of the node chart.'''
def getColor(family):
    familyDict = {'fam_1':'pink', 'fam_2':'red','fam_3':'orange','fam_4':'yellow','fam_5':'green','fam_6':'cyan','fam_7':'purple','fam_8':'brown', 'degen':'black'}
    return familyDict[family]


'''2D list -> 1D dictionary
gets all the nodes in the clusters, numbers them, and make
searchable by dictionary'''
def getNodes(clusters):
    count = 0
    nodes = {}
    for line in clusters:
        for item in line:
            nodes[item] = count
            count += 1
    return nodes


'''2D list, dictionary -> list
gets all the edges representing the clusters'''
def getEdges(clusters, nodes):
    edges = []
    for line in clusters:
        for index, e1 in enumerate(line):
            for e2 in line[index:]:
                edges.append((nodes[e1],nodes[e2]))
    return edges

'''String -> dictionary
gets the domain info for each node from the
hmmr files.'''                
def getNodeDomains(hmmrDirectory):
    files = [ x for x in os.listdir(hmmrDirectory) if os.path.isfile("/".join([hmmrDirectory,x]))]
    nodes = {}
    for file in files:
        with open("/".join([hmmrDirectory, file])) as input:
            data = input.read()
            proteins = re.findall("(?s)>>(.+?)\s+(.+?)\n(.+?)(?=\n>>|$)", data)
            for protein in proteins:
                strain = re.match("(.+?)_.*", protein[0]).group(1)
                name = protein[1]
                srs = re.search("SRS.*$", name)
                if srs:
                    name = srs.group(0)
                name = "{0}_{1}".format(strain, name)
                    
                domains = re.findall("Family\s+(.+?)\s", protein[2])
                degens = re.findall("Degen\s+(.+?)\s", protein[2])
                
                fixedDomains = []
                for domain, degen in zip(domains, degens):
                    if degen=='True':
                        fixedDomains.append("degen")
                    elif degen=='False':
                        fixedDomains.append(domain)
                    else:
                        raise Exception("bad degen value")
                
                nodes[protein[0]] = (name, fixedDomains)
    return nodes

# '''dict (nodes) -> dict
# aligns the domains to easily spot differences.'''
# def alignDomains(nodes):
#     



def encode(clusters, hmmrDirectory, name):
    print(name + " is being encoded")
    #nodes are two-tuples consisting of id and label
    nodes = getNodes(clusters)
    edges = getEdges(clusters, nodes)
    nodeComp = getNodeDomains(hmmrDirectory)
    
    #prep the texts
    nodeTexts = [getNodeText(index, nodeComp[node]) for node, index in sorted(nodes.items(), key=lambda x: x[1])]
    edgeTexts = [getEdgeText(edge) for edge in sorted(edges)]
    
    print("{0} nodes and {1} edges. Writing...".format(len(nodeTexts), len(edgeTexts)))
    
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
    </graph>".format(name, "\n".join(nodeTexts), "\n".join(edgeTexts))
    return text


def getNodeText(index, infoTuple):
    text = "\
        <node label=\"{1}\" id=\"{0}\" >\n\
            <att name=\"name\" type=\"string\" value=\"{1}\"/>\n\
            <att name=\"Gradient\" type=\"string\" value=\"stripechart: colorlist=&quot;{2}&quot;\"/>\n\
            <graphics h=\"10\" w=\"10\"/>\n\
        </node>\n".format(index, infoTuple[0], ",".join([getColor(x) for x in infoTuple[1]]))
    return text


'''(String, String, int) -> String
same.'''
def getEdgeText(edge):
   
    if edge[0] == edge[1]:
        return ""
    
    text = "\
        <edge source=\"{0}\" target=\"{1}\" >\n\
            <graphics width=\"1\"/>\n\
        </edge>\n".format(edge[0], edge[1])
    return text
