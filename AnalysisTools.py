'''
An aggregate script for independent anaysis functions
'''
import pandas
import numpy as np
import multiprocessing as mp
from colorsys import hls_to_rgb
from functools import reduce
from copy import deepcopy


def clustersToMatrix(clusters, sample_list):
    '''
    [[[string]]] -> [DataFrame]
    takes a bunch of cluster results and make them into similarity matrices
    '''
    #one base matrix to copy
    base_matrix = pandas.DataFrame(np.zeros((len(sample_list), len(sample_list))), index = sample_list, columns = sample_list)
    pool = mp.Pool(None, clusterWorkerInit, (base_matrix,))
    res = pool.map(clusterWorker, clusters)
    return res

def clusterWorkerInit(_base):
    global base 
    base = _base 

#decompose one cluster
def clusterWorker(cluster):
    '''
    does it for one cluster
    inherit base from init
    '''
    scaffold = deepcopy(base)
    for row in cluster:
        scaffold.loc[row, row] = 1
    
    return scaffold

def overallMatrix(matrices):

    return reduce(lambda x, y: x + y, matrices) 

def createColorTable(groups, overall_clusters, sample_list):
    '''
    [string], [[string]], [string] > df
    though really sample list is just a flat overall_clusters
    might as well keep it consistent. 
    takes the groups and the corresponding clusters and returns
    a df that is samples x 1, filled with what color each sample should be
    '''
    #figure out the hue and luminosity values
    n_groups = len(groups)
    n_hue_pref = 9
    hue_offset = 20
    lum_sets = [[0.5], [0.75, 0.25], [0.75, 0.50, 0.25]]

    n_lums = int(np.ceil(n_groups / n_hue_pref))
    print(n_lums)
    lums = lum_sets[n_lums - 1]
    
    n_hues = int(np.ceil(n_groups / n_lums))
    hue_max = 360-hue_offset * 2
    hue_gap = hue_max // n_hues
    hues = [(i * hue_gap + hue_offset)/360 for i in range(n_hues)]

    
    rgb_vals = [hls_to_rgb(hue, lum, 1) for lum in lums for hue in hues]
    #construct a table
    df = pandas.DataFrame('', columns = ['color'], index = groups)

    for group, rgb in zip(groups, rgb_vals[:len(groups)]):
        df['color'].loc[group] = rgb
    
    df.loc['None', 'color'] = (0,0,0)
    
    return df
