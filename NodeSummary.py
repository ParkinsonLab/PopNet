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
import matplotlib as mpl
import optiResultReader as orr
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import resource
import time


def averageFeatureLength(directory):
    sl = 10000
    results_dict = []
    
    files = searchForFiles(directory)
    for file in files:
        with open(file, 'r') as input:
            raw_text = input.read()
            
        features_dict = findSelfFeatures(raw_text)
        values = []
        for name in features_dict:
            values += features_dict[name][1] #each slot is (gradient, values)
    
        print('{0} has average length of {1}\n'.format(file, np.mean(values)))


#finds all the "features", i.e. all the block that are not self or 
def findFeatures(raw_text):
    
    regex = '(?ims)<node.+/node>'
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    
    features_dict = {}
    
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        col = getColor(node)
        
        #gets rid of the blanks, leaving only "features"
        to_remove = []
        for index, (color, value) in enumerate(zip(gradient, values)):
            if color == '#000000':
                to_remove.append(index)
        
        filtered_grad = [x for ind, x in enumerate(gradient) if ind not in to_remove]
        filtered_vals = [x for ind, x in enumerate(values) if ind not in to_remove]
        
        features_dict[name] = (filtered_grad, filtered_vals)
    
    return features_dict


def findSelfFeatures(raw_text):
    
    regex = '(?ims)<node.+/node>'
    nodes = re.split('</node>', re.search(regex, raw_text).group(0))[:-1] #last element is empty string
    
    features_dict = {}
    
    for node in nodes:
        gradient = getGradient(node)
        values = getValues(node)
        name = getName(node)
        col = getColor(node)
        
        #gets rid of the blanks, leaving only "features"
        to_remove = []
        for index, (color, value) in enumerate(zip(gradient, values)):
#             if color == col or color == '#000000': #NOT DELETING SELF
            if color == '#000000':
                to_remove.append(index)
        

        filtered_grad = [x for ind, x in enumerate(gradient) if ind not in to_remove]
        filtered_vals = [x for ind, x in enumerate(values) if ind not in to_remove]
        
        features_dict[name] = (filtered_grad, filtered_vals)
    
    return features_dict

# def findFeatureDistribution(features_dict):
#     
#     values = []
#     for name in features_dict:
#         values += features_dict[name][1] #each slot is (gradient, values)
#     
#     values = np.array(values)
#     mean = np.mean(values)
#     std = np.std(values)
#     
#     return (mean, std)    



def findFeatureDistribution(features_dict):
    #currently modded so that it's not percentage, but real numbers
    
    values = []
    for name in features_dict:
        values += features_dict[name][1] #each slot is (gradient, values)
    
    values = sorted(values)
    dist_dict = {}
    
#     #for counting features only
#     dist_dict[1] = len(values)
    
#     total = sum(values)
    total = 1
      
    count = 0
    pre = values[0] #UNCOMMENT AFTER
    dist_dict[pre] = 0
    for value in values:
        if value not in dist_dict:
            dist_dict[pre] = float(count) * pre / total
            count = 0
            dist_dict[value] = 0
        count += 1
        pre = value #UNCOMMENT AFTER
    dist_dict[pre] = float(count) * pre / total
    
    return dist_dict


def binResults(weighted_dict, bins):
    #first bin's always zero. not real. 
    results = {k:0. for k in bins[1:]}
    
    for k, v in weighted_dict.items():
        bin = min(bins[1:], key=lambda x:abs(x-k))
        results[bin] += v
        
    return results

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
            hex_str = '{:0<6}'.format(parts[0])
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

def makeLineGraph(coords_list, labels_list, title, subtitle_list, file_name):
    
    if len(coords_list) > 1:
        width = 2
        length = len(coords_list) // 2 + len(coords_list) % 2
    else:
        width = 1
        length = 1
    
    fig = mpl.figure.Figure()
    canvas = FigureCanvas(fig)
    fig.suptitle(title)
    
    for index, (coords, labels, subtitle) in enumerate(zip(coords_list, labels_list, subtitle_list)):
        ax = fig.add_subplot(length,width,index + 1)
        ax.set_title(subtitle)
        ax.set_xlabel(labels[0])
        ax.set_ylabel(labels[1])
        ax.plot(coords[0], coords[1])
        ax.set_xticklabels(ax.get_xticks(), rotation=45, ha='right')

    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()

def makeStackedBarGraph(coords_list, bar_labels_list, stack_label_list, title, axis_labels, file_name):    
    
    bar_width = 0.2
    step = 2000
    inds = np.arange(len(bar_labels_list))
    
    fig = mpl.figure.Figure(figsize=(15,10))
    canvas = FigureCanvas(fig)
    fig.suptitle(title)
    fig.subplots_adjust(wspace = 0.5, hspace = 0.5)
    
# for color bar
#     colorRange, colorTable = orr.generateColorRange(stack_label_list, step=step)
    
#     norm = Normalize(vmax = max(stack_label_list), vmin = min(stack_label_list), clip = True)
#     cmap_xvalues = norm(xrange(min(stack_label_list), max(stack_label_list), step))
#     cmap_xvalues[0] = 0.
#     cmap_xvalues[-1] = 1
#     cmap_list = [(val, c) for val, c in zip(cmap_xvalues, orr.generateColorRange(stack_label_list, raw=True, step=step)[1])]


    
#     cmap = LinearSegmentedColormap.from_list('mycmap', cmap_list)
#     sm = ScalarMappable(norm=norm, cmap = cmap)
#     sm.set_array(stack_label_list)
    
    colorRange, colorTable = orr.generateColors(stack_label_list, [0], step=step)
    
    ax = fig.add_subplot(2,2,1)
    ax.set_title(title, fontsize=28)
    ax.set_xlabel(axis_labels[0], fontsize=24)
    ax.set_ylabel(axis_labels[1], fontsize=24)

    pre = coords_list[0][1]
    bars = []
    bottom = np.zeros((len(coords_list[0][1]),))
    for i, (coords, stack_label) in enumerate(zip(coords_list, stack_label_list)):
#         if i == 0:
#             bar = ax.bar(inds, coords[1], color=colorTable[colorRange.index(coords[0])], align='center', linewidth=0)
#         else:
#             bar = ax.bar(inds, coords[1], color=colorTable[colorRange.index(coords[0])], align='center', bottom = bottom, linewidth=0)
        
        if i == 0:
            bar = ax.bar(inds, coords[1], color=colorTable[i], align='center', linewidth=0)
        else:
            bar = ax.bar(inds, coords[1], color=colorTable[i], align='center', bottom = bottom, linewidth=0)
            
        bars.append(bar)
        bottom = bottom + np.array(coords[1])    
    
    ax.set_xticks(inds)       
    ax.set_xticklabels(bar_labels_list, rotation=45, ha='right')
    ax.tick_params(axis='both', which='both', labelsize=18)
    ax.legend(bars, stack_label_list, loc=1, bbox_to_anchor=(1.2, 1), ncol=1, markerscale=0.2, fontsize=6)
#     fig.colorbar(sm)
    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()
    
    print('Graphing Resource usage: {}kb'.format(getattr(resource.getrusage(resource.RUSAGE_SELF), 'ru_maxrss') / 1000))
    
    
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
        
def findDistribution(directory, bins, title, axes, filepath):
    '''
    main runner method for finding that distribution. Requires all xgmml files
    to be in the same folder, named as 'IXXPIXXSXX.xgmml'
    
    currently changed the graph labels!
    '''
    
    #helpers
    def fillDict(inds, dict):
        for i in inds:
            if i not in dict:
                dict[i] = 0
        
    def transform(inds, dict_list):
        results = []
        for i in inds:
            tmp = []
            for k in dict_list:
                tmp.append(k[i])
            results.append((i, tmp))
        
        return results
    
    def combine(dict_list):
        result = {}
        for k in dict_list[0].keys():
            result[k] = sum([dict[k] for dict in dict_list]) / len(dict_list)
        return result
    #path patterns. Append the folder path
    
    results_dict = {}
    np_results_list = []
    i_set = set([])
    pi_set = set([])
    sl_set = set([])
    inds = set()
    
    #bins!
#     bins = [0,10000,20000,40000,60000,80000,120000,140000,160000,180000,200000]
#     bins = [0,20000,50000,100000,200000,500000,1000000,1500000,5000000]
    bins = [0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000]
    
    files = searchForFiles(directory)
    for file in files:
        
        i, pi, section_length = parseName(file)
        i_set.add(i)
        pi_set.add(pi)
        sl_set.add(section_length)
        
        with open(file, 'r') as input:
            raw_text = input.read()
#         #    Change this back!
#         features_dict = findFeatures(raw_text)
        features_dict = findSelfFeatures(raw_text)
        results = findFeatureDistribution(features_dict) #returns a dict {length: weight}
        
        
#         weighted_results = {k * section_length: results[k] for k in results} #NORMAL
        weighted_results = {k * 10000: results[k] * 10000 / 67 for k in results} #USE ONLY FOR THAT ONE GRAPH FOR YEAST. IF YOU DON'T RMB DON'T USE
#         weighted_results = {k: results[k] / 67 for k in results} #USE ONLY FOR THAT ONE GRAPH FOR YEAST. IF YOU DON'T RMB DON'T USE

        binned_results = binResults(weighted_results, bins)
        results = binned_results
        try:
            results_dict[section_length].append(results)
        except KeyError:
            results_dict[section_length] = []
            results_dict[section_length].append(results)
        
        for k, v in results.items():
            inds.add(k)   
    
    sl_set = sorted(sl_set)
    inds = sorted(inds)
    
    transformed_dict_list = []
    for sl in sl_set:
        dicts = results_dict[sl]
        for dict in dicts:
            fillDict(inds, dict)
        results_dict[sl] = combine(dicts)
    
    transformed_list = transform(inds, [results_dict[k] for k in sl_set])

# This is the label that is changed            
#     makeStackedBarGraph(transformed_list, sl_set, inds, 'Distribution of Features', ('Section Length', 'Features'),'{}/binned_features.pdf'.format(output_directory))
    makeStackedBarGraph(transformed_list, sl_set, inds, title, axes, filepath)

    
        
            
# The old line Style
#         with open(file, 'r') as input:
#             raw_text = input.read()
#             features_dict = findFeatures(raw_text)
#             results = findFeatureDistribution(features_dict) #(mean, std)
#             means_list.append((i, pi, section_length, results[0]))
#             std_list.append((i, pi, section_length, results[1]))
#         
#     means_arrdict = {sl: np.zeros((len(i_set), len(pi_set))) for sl in sl_set}
#     std_arrdict = {sl: np.zeros((len(i_set), len(pi_set))) for sl in sl_set}
#     i_set = sorted(i_set, reverse = True)
#     pi_set = sorted(pi_set)
#     
#     for i, pi, sl, mean in means_list:
#         means_arrdict[sl][i_set.index(i), pi_set.index(pi)] = mean
#     
#     for i, pi, sl, std in std_list:
#         std_arrdict[sl][i_set.index(i), pi_set.index(pi)] = std
#     
# #     for sl in sl_set:
# #         ag.graphHeatMap('{}/features_S{}.pdf'.format(output_directory, sl), [means_arrdict[sl], std_arrdict[sl]], [[min(pi_set), max(pi_set), min(i_set), max(i_set)]] * 2, ['Means', 'Standard Deviations'], 'Features Metric for SL {}'.format(sl))
#     
#     means_x = list(sorted(sl_set))
#     means_y = [np.average(means_arrdict[sl]) for sl in sorted(means_arrdict)]
#     means_coords = (means_x, means_y)
#     
#     std_x = means_x
#     std_y = [np.average(std_arrdict[sl]) for sl in sorted(std_arrdict)]
#     std_coords = (std_x, std_y)
#     
#     makeLineGraph([means_coords, std_coords], [('Section Length', 'Average Mean Size of Features'), ('Section Length', 'Average Standard Deviation in Feature Size')],'Feature Metric', ['Means', 'Standard Deviation'], '{}/features.pdf'.format(output_directory))
#     
    
if __name__ == '__main__':
# #    Test Code
#     a = "<graph label=/'Genome/' directed=/'0/'>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<node label=/'-_-_A/' id=/'0/'></node>\n\t<edge source=/'1/'"
      
#     directory = '/data/new/javi/yeast/pipeline/winvar2/SL10000'
#     input_path = directory + '/cytoscape/cytoscapeGenome.xgmml'
#     group_path = directory + '/groups.txt'
#     color_path = directory + '/colors.txt'
#     out_path = directory + '/node_summary.txt'
#     makeSummary(input_path, group_path, color_path, out_path)

# Currently used for find only nonself. Change in findDistribution and findFeaturesDistribution
    directory = '/data/new/javi/yeast/pipeline/penalty/networks2'
    output_directory = '/data/new/javi/yeast/pipeline/penalty/'
    start_time = time.time()
    
    #Various parameters for the graph
    title = 'Distribution of Features'
    x_axis = 'Max Gap Length'
    y_axis = 'Number of Base Pairs'
    filename = 'max_basic'
    
    #Vary the size of bins if too much of your data falls into one or a few of them.
    bins = [0,100000,200000,300000,400000,500000,600000,700000,800000,900000,1000000,2000000]
    
    #Don't edit
    path = '{0}/{1}.pdf'.format(output_directory, filename)
    axes = (x_axis, y_axis)
    
    findDistribution(directory, bins, title, axes, path)
    print('Graph took {} seconds to make.'.format(time.time() - start_time))
    print('Run Complete')

# # Currently used for find only self. Change in findDistribution and findFeaturesDistribution
#     directory = '/data/new/javi/assorted_networks'
#     output_directory = '/data/new/javi/'
#     start_time = time.time()
#     averageFeatureLength(directory)
#     print('Graph took {} seconds to make.'.format(time.time() - start_time))
#     print('Run Complete')
