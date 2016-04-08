'''
Created on Jan 26, 2016

@author: javi

This module goes over the contents of each node and output a summary of its composition,
including the amount of foreign material in there. Reads from an xgmml file.
'''
import re
import GroupComposition as gc

#finds all the "features", i.e. all the block that are not self or 
def findFeatures(raw_text, colorTable):
    
    regex = '(?ims)<node.+/node>'
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    
    features_dict = {}
    
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        
        #gets rid of "self" and the blanks, leaving only "features"
        while '#000000' in gradient:
            gradient.remove('#000000')
        
        while colorTable[name] in gradient:
            gradient.remove(colorTable[name])
        
        features_dict[name] = gradient
    
    return features_dict
        
    
    

def parseXGMML(raw_text, groups):
    
    '''
    returns a nested dictionary {node_name:{group:value}}
    '''
    
    regex = '(?ims)<node.+/node>'
    rev_groups = gc.reverseGroups(groups)
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    values_dict = {}
            
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        
        values_dict[name] = {entry: 0 for entry in gradient}
        for entry, val in zip(gradient, values):
            values_dict[name][entry] += val
        
    return values_dict

#These methods work for a Single Block!
def getGradient(raw_text):
    regex = '(?ims)Gradient.+?(#.+?)&quot'
    return re.split(',', re.search(regex,raw_text).group(1))
def getValues(raw_text):
    regex = '(?ims)values(.+?)/att'
    regex2 = '[0-9]+'
    temp_text = re.search(regex, raw_text).group(1)
    return [int(x) for x in re.findall(regex2, temp_text)]
def getName(raw_text):
    regex = '(?i)name.+?value=\"(.+?)\"'
    return re.search(regex, raw_text).group(1)

#

def loadColorTable(color_path, reverse = False):
    '''
    the reverse function is much like the reverse groups.
    Because the whole reason for this is to support the make
    Summary function, which needs the reverse one.
    
    Normal = name: color
    '''
    
    with open(color_path) as input_type:
        data = input_type.read()
    
    results = {}
    
    for line in re.split('\n', data)[:-1]: #last line is empty.
        parts = re.split('\t', line)
        if reverse == True:
            results[parts[1]] = parts[0]
        else:
            results[parts[0]] = parts[1]
    
    return results
    
def makeSummary(input_path, group_path, color_path, out_path):
        
    groups = gc.loadGroups(group_path, '')
    del groups['NONE']
    rev_groups = gc.reverseGroups(groups)
    colors = loadColorTable(color_path)
    
    with open(input_path) as input_type:
        data = parseXGMML(input_type.read(), groups)
        
    with open(out_path, 'w') as output:
        
        group_names = sorted(groups.keys())
        output.write('\t'.join(['', ''] + group_names))
        output.write('\n')
        
        for group in group_names:
            for node in sorted(groups[group]):
                del data[node]['#000000']
                output.write('{0}\t{1}\t'.format(node, rev_groups[node]))
                total = float(sum(data[node].values()))
                for group in group_names:
                    try:
                        output.write('{0}\t'.format(str(data[node][colors[group]] / total * 100)))
                    except(KeyError):
                        output.write('0\t')
                output.write('\n')
    
    print('Node Summaries Completed')
        


if __name__ == '__main__':
#     a = "<graph label=/'Genome/' directed=/'0/'>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<edge source=/'1/'"
    
    directory = '/data/new/javi/yeast/pipeline'
    input_path = directory + '/matrix/cytoscape/cytoscapeGenome.xgmml'
    group_path = directory + '/matrix/groups.txt'
    color_path = directory + '/matrix/colors.txt'
    out_path = directory + '/matrix/node_summary.txt'
    makeSummary(input_path, group_path, color_path, out_path)
