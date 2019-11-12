'''
Created on Dec 16, 2014

@author: javi
'''
'''uses the whole genome MCLC data to group things automatically. by MCL'''

import time
import string
import multiprocessing as mp
import numpy as np
from functools import reduce
from sklearn.metrics import silhouette_score


import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

from IOTools import mcl, readTab, buildMatrix, writeGroup

def group(similarity_matrix, tab_path, out_path, params, autogroup=False):

    def getGroupName(n):
        '''
        get the nth groups name following a pattern
        '''
        if n > 26*26:
            raise ValueError('group::getGroupName encountered too many groups to name')
        chars = list(string.ascii_uppercase)
        m = len(chars)

        if n < m:
            return chars[n]
        else:
            return chars[int(n/m)] + chars[n % m]
        

    sample_list = readTab(tab_path)
    mci_string = buildMatrix(similarity_matrix.values, sample_list)

    if autogroup:
        params_to_use = optimize(similarity_matrix, tab_path, params)
    else:
        params_to_use = params
    
    cluster = mcl(mci_string, tab_path, params_to_use.getIVal(), params_to_use.getPiVal())
    names = [getGroupName(n) for n in range(len(cluster))]
    writeGroup(names, cluster, out_path)

    return names, cluster
            
def optimize(matrix, tab_path, params):
    '''
    find the best thing according silhouette index
    '''
    print('Beginning Optimization')    
    #list of (ival, pival) from small to large
    sample_list = readTab(tab_path)
    params_to_try = []
    i_list = stepList(params.getIMin(), params.getIMax(), params.getIStep(), reverse=True)
    pi_list = stepList(params.getPiMin(), params.getPiMax(), params.getPiStep())
    for x in i_list:
        for y in pi_list:
            params_to_try.append((x, y))
    
    pool = mp.Pool(None, getMetricsInit, (matrix, tab_path, sample_list))
    metrics = pool.starmap(getMetricsWorker, params_to_try, chunksize=10)
    n_clusters, silhouette = zip(*metrics)

    #generate graph
    output_name = 'Heatmaps.pdf'
    graph_title = 'S2 Clustering Metrics'

    cln_matrix = np.array(n_clusters).reshape(len(i_list), len(pi_list))
    sil_matrix = np.array(silhouette).reshape(len(i_list), len(pi_list))
    
    matricies = [cln_matrix, sil_matrix]
    extents = [[params.getPiMin(),params.getPiMax() + params.getPiStep(),params.getIMin(),params.getIMax() + params.getIStep()]]*2
    names = ['Number of Clusters', 'Silhouette Index']
    axis_labels = [('pI value', 'I value')] * 2
    graphHeatMap(output_name, matricies, extents, names, axis_labels, graph_title)

    #For silhouette, higher is better. return the best parameters
    best = np.argmax(silhouette)
    print(silhouette)
    
    ival, pival = params_to_try[best]
    params.setIVal(ival)
    params.setPiVal(pival)
    print('Optimization selected values I = {0}, pI = {1}'.format(ival, pival))
    print('Finished Optimization')
    return params

def getMetricsInit(_matrix, _tab_path, _sample_list):
    global matrix
    global tab_path
    global sample_list
    global mci_string
    matrix = _matrix.values
    tab_path = _tab_path 
    sample_list = _sample_list
    mci_string = buildMatrix(matrix, sample_list)

def getMetricsWorker(ival, pival):
    '''
    get the matrix and tab path from parent function
    then calculate some metrics

    returns (N, silhouette index)

    matrix, sample_list from parent

    currently we're doing silhouette index, higher is better
    '''
    def getLabels(y):
        label = y[0]
        names = y[1]
        return {name: label for name in names}

    global matrix 
    global tab_path 
    global sample_list
    global mci_string
    cluster = mcl(mci_string, tab_path, ival, pival)
    label_dict = reduce(lambda x, y: {**x, **getLabels(y)}, [(i, e) for i, e in enumerate(cluster)], {})
    labels = [label_dict[label] for label in sample_list]
    try:
        score = silhouette_score(matrix, labels, metric='precomputed')
    except ValueError:
        score = -1

    return (len(cluster), score)

def stepList(start, end, step, reverse=False):
    if start == end:
        return [start]
    list = []
    curr = start
    while curr <= end:
        list.append(round(curr, 1))
        curr += step
    
    if reverse:
        return [x for x in reversed(list)]
    else:
        return list
       
def graphHeatMap(file_name, matrices, extents, names, axis_labels, title):
    '''current style only gens 2x2 heat maps! change length, width to make
    different kinds'''

    length = 2
    width = 1
    axes = []
    
    fig = figure.Figure(figsize = (15, 15))
    canvas = FigureCanvas(fig)
    fig.subplots_adjust(wspace = 0.45, hspace = 0.45)
    fig.suptitle(title)
    
    for matrix, extent, name, ind, axis_label in zip(matrices, extents, names, range(len(matrices)), axis_labels):
#         ax = plt.subplot2grid((width, length), (ind // 2, ind % 2))
        ax = fig.add_subplot(width, length, ind + 1)
        ax.set_title(name, fontsize=24)
        ax.set_xlabel(axis_label[0], fontsize=24)
        ax.set_ylabel(axis_label[1], fontsize=24)
        img = ax.imshow(matrix, extent = extent, interpolation = 'none', vmin = np.nanmin(matrix), vmax = np.nanmax(matrix), cmap='bwr')
        ax.tick_params(axis='both', which='both', labelsize=18)
        fig.colorbar(img)
    
    
    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()
    