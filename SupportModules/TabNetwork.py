'''
Created on Jul 14, 2016

@author: javi

This script aims to get the tab network (similar to tabNetwork.tsv) from the xgmml file.

Probably uses many methods from Node Summary
'''
import NodeSummary as ns
import GroupComposition as gc
import re 
import ChrNameSorter as cns
'''
similar to the NodeSummary one, but we care about different things here
'''
def parseXGMML(file_path, tab_path):
    
    with open(file_path) as input:
        raw_text = input.read()
    
    groups = gc.loadGroups(tab_path, '')
    
    
    regex = '(?ims)<node.+/node>'
    
    rev_groups = gc.reverseGroups(groups)
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    values_dict = {}
    meta_dict = {}
    
    for node in nodes:
        gradient = ns.getGradient(node)
        values = ns.getValues(node)
        name = ns.getName(node)
        group = ns.getGroup(node)
        color = ns.getColor(node)
        values_dict[name] = []
        for entry, val in zip(gradient, values):
            if not entry == '#000000':
                values_dict[name] += [entry] * val
        meta_dict[name] = (group, color)
        
    return values_dict, meta_dict

'''
Takes the whole tuple from parseKGMML
'''

def makeTabNetwork(parsed_result, guide):
    
    values, meta = parsed_result
    strains = list(values.keys())
    chr_list = sorted(guide.keys(), key = lambda x: cns.getValue(x, 'plasmodium'))
    result = ['\t'.join(['', ''] + strains), '\t'.join(['', ''] + [meta[x][0] for x in strains])]
    
    count = 0
    for chr in chr_list:
        for pos in guide[chr]:
            result.append('\t'.join([chr, str(pos)] + [values[x][count] for x in strains]))
            count+=1
    
    return '\n'.join(result)
    

def loadGuide(guide_path):
    chr_exp = '(?s)@(.+?)\n(.+?)(?=@|$)'
    block_size = 10000
    results = {}
    
    with open(guide_path, 'r') as input:
        data = input.read()
        
    for chr in re.findall(chr_exp, data):
        name = chr[0]
        count = chr[1].count('#')
        
        results[name] = []
        for x in range(count):
            results[name].append(x * block_size)
    
    return results
    
    
    
if __name__ == '__main__':
    directory = '/data/new/javi/plasmo/pipeline/matrix'
    xml_path = '/'.join([directory, 'cytoscape', 'cytoscapeGenome.xgmml'])
    tab_path = '/'.join([directory, 'persistentMatrix.tab'])
    guide_path = '/'.join([directory, 'persistentResult.txt'])
    out_path = '/'.join([directory, 'tabNetwork.tsv'])
    groups = gc.loadGroups(tab_path, '')
    guide = loadGuide(guide_path)
    
    with open(out_path, 'w') as output:
        output.write(makeTabNetwork(parseXGMML(xml_path, tab_path), guide))
        
    print('TabNetwork Finished.')