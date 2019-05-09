'''
Created on May 3, 2016

@author: javi
'''

import re
import numpy as np
import matplotlib as mpl
import os
from os.path import isdir
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def readFile(file_path):
    
    data = np.genfromtxt(file_path, delimiter='\t', skip_header=1)[:,2:]
    
    return data

    
def makeBarGraph(coords_list, bar_labels_list, title, axis_labels, file_name):    
    
    bar_width = 0.1
    step = 2000
    inds = np.arange(len(bar_labels_list))
    
    fig = mpl.figure.Figure(figsize=(15,10))
    canvas = FigureCanvas(fig)
    fig.suptitle(title)
    fig.subplots_adjust(wspace = 0.5, hspace = 0.5)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(axis_labels[0])
    ax.set_ylabel(axis_labels[1])

    bar = ax.bar(inds, coords_list, align='center', linewidth=0, width=bar_width)
    
#     bar_labels_list = [0] + bar_labels_list
    ax.set_xticklabels(bar_labels_list, rotation=45, ha='right')
    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()

def rangeList(inc, n):
    #helper for multi bar graph.
    
    if n == 1:
        return [0]
    
    results = []
    min = -1. * inc * (n - 1) / 2
    
    for i in xrange(n):
        results.append(min + inc * i)
    
    return results
    
    
def makeMultiBarGraph(coords_dict, bar_labels_list, title, axis_labels, file_name):
    
    bar_width = 0.1
    inds = np.arange(len(bar_labels_list))
    inc = 0.15
    steps = rangeList(inc, len(coords_dict))
    color_list = ['red', 'blue', 'green', 'yellow', 'purple', 'aqua', 'black', 'maroon', 'olive', 'lime']

    
    fig = mpl.figure.Figure(figsize=(15,10))
    canvas = FigureCanvas(fig)
    fig.suptitle(title)
    fig.subplots_adjust(wspace = 0.5, hspace = 0.5)
    
    ax = fig.add_subplot(1,1,1)
    ax.set_xlabel(axis_labels[0], size=14)
    ax.set_ylabel(axis_labels[1], size=14)
    
    bars = []
    for step, color, (sl, coords) in zip(steps, color_list, sorted(coords_dict.items())):
        bar = ax.bar(inds + step, coords, align='center', linewidth=0, width=bar_width, label=str(sl), color=color)
        bars.append(bar)
#     bar_labels_list = [0] + bar_labels_list
    ax.set_xticks(inds)
    ax.set_xticklabels(bar_labels_list, rotation=45, ha='right')
    ax.tick_params(axis='both', which='both', labelsize=18)
    fig.legend(bars, sorted(coords_dict), 'upper right')
    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()
    
def multiDensity(directory, outpath):
    
    name_pattern = 'density_([0-9]+?).tsv'
    file_list = [x for x in os.listdir(directory) if re.search(name_pattern, x)]
#     bins = [0, 5, 10, 20, 50, 100, 200, 400, 600, 2000]
    bins = [0,1,2,3,4,5,10,20,50,100,1000]
    values_dict = {}
    for file in file_list:
        path = '/'.join([directory, file])
        sl = int(re.search(name_pattern, file).group(1))
        data = readFile(path)
        hist = np.histogram(data, bins=bins)[0]
        print(hist)
        hist = hist / float(data.shape[0]) #normalize
        
        print(data.shape)
        
        hist = np.insert(hist, 0, 0)
        values_dict[sl] = hist
    
    makeMultiBarGraph(values_dict, bins, 'SNP Densities', ('SNPs per Section', 'Fraction of Sections'), output_path)

def readClusters(file_path):
    '''
    for counting the distribution of clusters. reads persistentResult.txt
    '''
    
    with open(file_path, 'r') as input:
        data = input.read()
        sections = re.findall('(?s)#[0-9]+?\n(.+?)\n\n', data)
        data = None
    
    results = []
    for section in sections:
        lines = re.split('\n', section)[1:]
        for line in lines:
            size = len(line.split('\t'))
            results.append(size)
    
    results = np.array(results)
    
    return results

def multiClusterDist(directory, outpath):
    
    name_pattern = 'cluster_(.+?).txt'
    file_list = [x for x in os.listdir(directory) if re.search(name_pattern, x)]
    bins = [0, 1, 3, 4, 5, 6, 7, 8, 9, 10, 100]
    
    values_dict = {}
    for file in file_list:
        path = '/'.join([directory, file])
        try:
            sl = int(re.search(name_pattern, file).group(1))
        except:
            sl = re.search(name_pattern, file).group(1)
        data = readClusters(path)
        hist = np.histogram(data, bins=bins)[0]
        print(hist)
        hist = hist / float(data.shape[0]) #normalize
        
        print(data.shape)
        
        hist = np.insert(hist, 0, 0)
        values_dict[sl] = hist
    
    makeMultiBarGraph(values_dict, bins, 'Cluster Distribution', ('Number of Clusters', 'Fraction of Sections'), output_path)
    



    
if __name__ == '__main__':
#     filepath = '/data/new/javi/yeast/pipeline/winvar2/2K/density.txt'
#     output_path = '/data/new/javi/yeast/pipeline/winvar2/density_graph2K.pdf'
#     data = readFile(filepath)
#     print('max: {}, min: {}, average: {}'.format(np.max(data), np.min(data), np.average(data)))
# #     bins = np.logspace(0, np.log(np.max(data)), base = 10, num=10)
#     bins = [0, 5, 10, 20, 50, 100, 200, 400, 600, 2000]
#     hist, bin_edges = np.histogram(data, bins=bins)
#     makeBarGraph(np.insert(hist, 0, 0), bins, 'SNP Density at 2K', ('SNPs per Section', 'Number of Sections'), output_path)

    directory = '/data/new/javi/yeast/pipeline/winvar_negative/densities'
    output_path = '/data/new/javi/yeast/pipeline/winvar_negative/raw_density_graph.pdf'
    multiDensity(directory, output_path)
    
#     directory = '/data/new/javi/yeast/pipeline/winvar2/iClusters'
#     output_path = '/data/new/javi/yeast/pipeline/winvar2/icluster_graph.pdf'
#     multiClusterDist(directory, output_path)
    
    print('graph completed!')