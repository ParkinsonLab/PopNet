'''
Created on Jan 26, 2016

@author: javi

This module goes over the contents of each node and output a summary of its composition,
including the amount of foreign material in there. Reads from an xgmml file.
'''
import re
import GroupComposition as gc
import numpy as np
import os
import re
import AutoGrouper as ag

#finds all the "features", i.e. all the block that are not self or 
def findFeatures(raw_text):
    
    regex = '(?ims)<node.+/node>'
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    
    features_dict = {}
    
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        
        #gets rid of the blanks, leaving only "features"
        to_remove = []
        for index, (color, value) in enumerate(zip(gradient, values)):
            if color is '#000000':
                to_remove.append(index)
        
        for index in to_remove:
            del gradient[index]
            del values[index]
        
        features_dict[name] = (gradient, values)
    
    return features_dict
        
def findFeatureDistribution(features_dict):
    
    values = []
    for name in features_dict:
        values += features_dict[name][1] #each slot is (gradient, values)
    
    values = np.array(values)
    mean = np.mean(values)
    std = np.std(values)
    
    return (mean, std)    

def parseXGMML(raw_text, groups):
    
    '''
    returns two nested dictionaries {node_name:{group:value}} and {node_name:(group, color)}
    '''
    
    regex = '(?ims)<node.+/node>'
    rev_groups = gc.reverseGroups(groups)
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    values_dict = {}
    meta_dict = {}
    
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        group = getGroup(node)
        color = getColor(node)
        values_dict[name] = {entry: 0 for entry in gradient}
        for entry, val in zip(gradient, values):
            values_dict[name][entry] += val
        meta_dict[name] = (group, color)
        
    return values_dict, meta_dict

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
def getGroup(raw_text):
    regex = '(?i)group.+?value=\"(.+?)\"'
    return re.search(regex, raw_text).group(1)
def getColor(raw_text):
    regex = '(?i)color.+?value=\"(.+?)\"'
    try:
        return re.search(regex, raw_text).group(1)
    except:
        return '#000000'
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

def makeColorTable(meta_dict, reverse = False):
    '''
    makes the color table from info gathered from the xgmml file
    normal = group:color
    '''
    
    color_table = {}
    
    for name in meta_dict:
        group, color = meta_dict[name]
        if group not in color_table:
            if not reverse:
                color_table[group] = color
            else:
                color_table[color] = group
    
    return color_table

def searchForFiles(directory):
    '''
    reads through the directory looking for numpy objects.
    returns a list of Complete paths to those files
    '''
    
    file_list = ['{0}/{1}'.format(directory, file) for file in os.listdir(directory) if file.endswith('.xgmml')]
    
    return file_list

def parseName(file_path):
    '''
    return the I, Pi, and Type of data for one matrix
    Takes complete files paths
    returns (section_length, i, pi, type)
    '''
    name_pattern = 'I([0-9]+?)PI([0-9]+?)S([0-9]+?)[.]'
    
    file_name = file_path.split('/')[-1]
    re_result = re.search(name_pattern, file_name)
    i = float(re_result.group(1)) / 10
    pi = float(re_result.group(2)) / 10
    sect_length = int(re_result.group(3))

    
    return i, pi, sect_length
    
def makeSummary(input_path, group_path, color_path, out_path):
        
    groups = gc.loadGroups(group_path, '')
    del groups['NONE']
    rev_groups = gc.reverseGroups(groups)
    colors = loadColorTable(color_path)
    
    with open(input_path) as input_type:
        data, meta_data = parseXGMML(input_type.read(), groups)
        
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
        
def findDistribution(directory, output_directory):
    '''
    main runner method for finding that distribution. Requires all xgmml files
    to be in the same folder, named as 'IXXPIXXSXX.xgmml'
    '''
    #path patterns. Append the folder path
    
    means_list = []
    std_list = []
    i_set = set([])
    pi_set = set([])
    sl_set = set([])
    
    files = searchForFiles(directory)
    for file in files:
        
        i, pi, section_length = parseName(file)
        i_set.add(i)
        pi_set.add(pi)
        sl_set.add(section_length)
        
        with open(file, 'r') as input:
            raw_text = input.read()
            
        features_dict = findFeatures(raw_text)
        results = findFeatureDistribution(features_dict) #(mean, std)
        means_list.append((i, pi, section_length, results[0]))
        std_list.append((i, pi, section_length, results[1]))
        
    means_arrdict = {sl: np.zeros((len(i_set), len(pi_set))) for sl in sl_set}
    std_arrdict = {sl: np.zeros((len(i_set), len(pi_set))) for sl in sl_set}
    i_set = sorted(i_set, reverse = True)
    pi_set = sorted(pi_set)
    
    for i, pi, sl, mean in means_list:
        means_arrdict[sl][i_set.index(i), pi_set.index(pi)] = mean
    
    for i, pi, sl, std in std_list:
        std_arrdict[sl][i_set.index(i), pi_set.index(pi)] = std
    
    for sl in sl_set:
        ag.graphHeatMap('{}/features_S{}.pdf'.format(output_directory, sl), [means_arrdict[sl], std_arrdict[sl]], [[min(pi_set), max(pi_set), min(i_set), max(i_set)]] * 2, ['Means', 'Standard Deviations'], 'Features Metric for SL {}'.format(sl))
    
    
    
    
    
if __name__ == '__main__':
#     a = "<graph label=/'Genome/' directed=/'0/'>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<edge source=/'1/'"
    
#     directory = '/data/new/javi/yeast/pipeline'
#     input_path = directory + '/matrix/cytoscape/cytoscapeGenome.xgmml'
#     group_path = directory + '/matrix/groups.txt'
#     color_path = directory + '/matrix/colors.txt'
#     out_path = directory + '/matrix/node_summary.txt'
#     makeSummary(input_path, group_path, color_path, out_path)

    directory = '/data/new/javi/yeast/networks'
    output_directory = '/data/new/javi/yeast/network_graphs'
    findDistribution(directory, output_directory)
    print('Run Complete')
