'''
Created on Dec 16, 2014

@author: javi
'''
'''uses the whole genome MCLC data to group things automatically. by MCL'''

import subprocess as subp
from subprocess import Popen
import re
import math
import numpy as np
import os
import random
import SnpSorter
import time
import multiprocessing as mp
from shutil import copyfile
import resource
import functools as ft

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.figure as figure
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas

def group(mclcTree, tabpath, outpath, mcipath, S2_iVal, S2_piVal):
    import MCLCounter as mclc
    sampleList = sorted(mclc.loadSampleList(tabpath))
    
    currString = buildMatrix(mclcTree, sampleList)

#     #MCL
#     S2_iVal = analyzeClm(currString, tabpath)
#     if S2_iVal == 0:
#         raise RuntimeError('No suitable iValue found!')
#         import sys
#         sys.exit()

    with open(mcipath, 'w') as moutpath:
        moutpath.write(currString)
           
    result = mcl(currString, tabpath, S2_iVal, S2_piVal)
        
     
    with open(outpath, 'w') as output:
        try:
            for line in result:
                gname = line[0]
                output.write("@{0:s}\n{1:s}\n".format(gname, "\t".join(line)))
        except:
            print('Overall Clustering Failed!')
            import sys
            sys.exit()
            
# DEPRECATED
# '''(String, list) -> [[string]]
# runs mcl for the matrix over a range of value, to find the best one.
# 
# Current Scheme:
# finds the largest PI/I that does not result in a singleton cluster
# '''
# def findOptimalPattern(currString, tabpath):
#     iMax = 8
#     piMax = 20
#     step = 0.1
#     patterns = []
#     
#     for x in reversed(np.arange(0, piMax, step)):
#         temp = mcl(currString, tabpath, iMax, x)
#         if not hasSingleton(temp): return temp
#         
#     for x in reversed(np.arange(0, iMax, step)):
#         temp = mcl(currString, tabpath, x, 0)
#         if not hasSingleton(temp): return temp
    

    
'''(dict, list, num, num) -> [[string]]
helper function for repeatedly running mcl over an array of I and PI values.
'''            
def mcl(currString, tabpath, iValue, piValue, raw=False):
    
    def mcl_pi():
        return subp.check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_nopi():
        return subp.check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_notab_pi():
        return subp.check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_notab_nopi():
        return subp.check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '1'])
    
    tempname = 'temp{}.mci'.format(random.randint(1, 1000000))
    
    with open(tempname, 'w') as temp:
        temp.write(currString)
    
    if tabpath is None and piValue > 0:
        result = mcl_notab_pi()
    elif tabpath is None:
        result = mcl_notab_nopi()
    else:
        if piValue > 0:
            result = mcl_pi()
        else:
            result = mcl_nopi()
            
    result = bytes.decode(result)
    os.remove(tempname)

    if raw:
        return result
    
    results = [re.split("\t", line) for line in re.split("\n",result.rstrip('\n'))]
    
    
    return results


#helper functions for clmInfo
def generateMclName(i, pi):
    
    name_pattern = 'out.{}'
    file_name = name_pattern.format('{i:0>2d}{pi:0>3d}'.format(i=int(i * 10), pi=int(pi * 10)))
    
    return file_name

def write_mcl(currString, i, pi, file_name):
    
    with open(file_name, 'w') as tmp:
        tmp.write(mcl(currString, None, i, pi, raw=True))

def setupClm(currString, i, pi_list):
    
    file_list = []
    for pi in pi_list:
        file_name = generateMclName(i, pi)
        write_mcl(currString, i, pi, file_name)
        file_list.append(file_name)
    
    return file_list

def setupClm_unpack(param_tuple):
    '''wrapper for calling by pool'''
    
    currString = param_tuple[0]
    i = param_tuple[1]
    pi_list = param_tuple[2]
    return setupClm(currString, i, pi_list)
###
    
# def clminfo(currString, tabpath, params):
#     '''runs clminfo in with a range of values, and finds the pattern with the least amount of change over time.'''
#     
#     #helper functions
#     def generateParams(currString, i_list, pi_list):
#         '''
#         generates the parameter list needed for running setupClm using a pool
#         '''
#         
#         return ((currString, i, pi_list) for i in i_list)
#     ###
#     
#     matrixName = 'in.mci'
#     with open(matrixName, 'w') as tmp:
#         tmp.write(currString)
#         
#     iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
#     piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
# 
#     pool = mp.Pool(4)
#     pool_params = generateParams(currString, iVals, piVals)
#     
#     files = pool.map(setupClm_unpack, pool_params)
#     pool.close()
#     pool.join()
#     
#     files = reduce(lambda x, y: x + y, files) #make 2d list flat
#     
#     result = subp.check_output(['clm', 'info', matrixName] + files, close_fds=True)
#     
#     #returns the clm dist (variance of information) but not using that atm.
# #    result = subp.check_output(['clm', 'dist', '-mode', 'sj', '--chain'] + files, close_fds=True)
#     
#     for file in files:
#         os.remove(file)
#     os.remove(matrixName)
#         
#     return result

def clminfo(currString, tabpath, params):
    '''runs clminfo in with a range of values, and finds the pattern with the least amount of change over time.'''
    
    #helper functions
    ###
    
    matrixName = 'in.mci'
    with open(matrixName, 'w') as tmp:
        tmp.write(currString)
        
    iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
    piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
    
    files = []
    for i in iVals:
        files += setupClm(currString, i, piVals)
#     files = reduce(lambda x, y: x + y, files) #make 2d list flat
    
    result = subp.check_output(['clm', 'info', matrixName] + files, close_fds=True)
    result = bytes.decode(result)
    #returns the clm dist (variance of information) but not using that atm.
#    result = subp.check_output(['clm', 'dist', '-mode', 'sj', '--chain'] + files, close_fds=True)
    
    for file in files:
        os.remove(file)
    os.remove(matrixName)
        
    return result

def analyzeClm(currString, tabpath, params):
     
    '''returns the optimal I value from 2 - 8, and pI value from 1 - 20 based on the stability, then based on efficiency.
    The maximum values are avoided in favor of any stable clustering patterns in the middle, but if it can't be helped then go for
    the highest possibility (I = 8, pi = 20). If it has to choo=se, it will overcluster rather than under. 
     
    (Actually not gonna do the pi thing)
     
    the out files are formatted as out.XXYY, where XX = I value * 10, and YY = piValue * 10. When the pi Value is not needed, it will default to 1.
    It seems that pi = 1 is like it's not there.'''
     
    start_time = time.time()
    effPattern = 'efficiency=(.+?)\s'
    clusPattern = 'clusters=(.+?)\s'
    namePattern= 'source=out\.(.+?)\s'
     
    data = clminfo(currString, tabpath, params)
    lines = re.split('===', data)
    #first position is efficiency, then the # of clusters, then the iValue, then the piValue used
    parsedLines = [(float(re.search(effPattern, line).group(1)), int(re.search(clusPattern, line).group(1)), float(re.search(namePattern, line).group(1)[:-3]) / 10, float(re.search(namePattern, line).group(1)[-3:]) / 10) for line in lines]
    clusNums = [x[1] for x in parsedLines]
     
     
    iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
    piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
    heat_cluster_matrix = np.zeros((len(iVals), len(piVals)))
    heat_eff_matrix = np.zeros_like(heat_cluster_matrix)
     
    for line in parsedLines:
        x = iVals.index(line[2])
        y = piVals.index(line[3])
        heat_cluster_matrix[x, y] = line[1]
        heat_eff_matrix[x, y] = line[0]
     
    if np.max(heat_eff_matrix) > 1:
        heat_eff_matrix = heat_eff_matrix * 0
         
    print('analyzeClm took {} seconds'.format(time.time() - start_time))
    return heat_cluster_matrix, heat_eff_matrix

# def analyzeClm(currString, tabpath, params):
#     
#     '''returns the optimal I value from 2 - 8, and pI value from 1 - 20 based on the stability, then based on efficiency.
#     The maximum values are avoided in favor of any stable clustering patterns in the middle, but if it can't be helped then go for
#     the highest possibility (I = 8, pi = 20). If it has to choose, it will overcluster rather than under. 
#     
#     (Actually not gonna do the pi thing)
#     
#     the out files are formatted as out.XXYY, where XX = I value * 10, and YY = piValue * 10. When the pi Value is not needed, it will default to 1.
#     It seems that pi = 1 is like it's not there.'''
#     
# #     start_time = time.time()
#     effPattern = 'efficiency=(.+?)\s'
#     clusPattern = 'clusters=(.+?)\s'
#     namePattern= 'source=out\.(.+?)\s'
#     
#     data = clminfo(currString, tabpath, params)
#     lines = re.split('===', data)
#     #first position is efficiency, then the # of clusters, then the iValue, then the piValue used
#     parsedLines = [(float(re.search(effPattern, line).group(1)), int(re.search(clusPattern, line).group(1)), float(re.search(namePattern, line).group(1)[:-3]) / 10, float(re.search(namePattern, line).group(1)[-3:]) / 10) for line in lines]
#     clusNums = [x[1] for x in parsedLines]
#     iVals = params.getIMin()
#     piVals = [params.getPiMin()]
#     heat_cluster_matrix = np.zeros((1,1))
#     heat_eff_matrix = np.zeros_like(heat_cluster_matrix)
#     heat_cluster_matrix[0,0] = parsedLines[0][1]
#     heat_eff_matrix[0,0] = parsedLines[0][0]
#     
#     if np.max(heat_eff_matrix) > 1:
#         heat_eff_matrix = heat_eff_matrix * 0
#     
# #     print('analyzeClm took {} seconds'.format(time.time() - start_time))
#     return heat_cluster_matrix, heat_eff_matrix
    
# #     diagnostics
#     import matplotlib
#     matplotlib.use('pdf')
#     import matplotlib.pyplot as plt
#     from matplotlib.backends.backend_pdf import PdfPages
#     x = [e[2] for e in parsedLines]
#     y = [e[1] for e in parsedLines]
#     z = [e[0] for e in parsedLines]
#     plt.subplot(3,1,1)
#     plt.plot(x, y)
#     plt.subplot(3,1,2)
#     plt.plot(x, z, 'r')
#     plt.show(block=False)
#     pdf = PdfPages('Autogroup.pdf')
#     pdf.savefig()
#     pdf.close()  
#     
#     
#     
#     return result
    
    
    
def buildMatrix(mclcTree, sampleList):
    
    
    def normalize(mclcTree):
        arr = np.array([list(x.values()) for x in mclcTree.values()])
        max = np.max(arr)
        min = np.min(arr)
        
        results = {}
        
        for k1, v1 in mclcTree.items():
            results[k1] = {}
            for k2, v2 in v1.items():
                results[k1][k2] = (v2 - min) / float(max - min)
        
        return results
    
    mclcTree = normalize(mclcTree)
    dimension = len(sampleList)    
    
    #Matrix Building
    currString = ""
    
    #The header portion of each matrix
    currString += "(mclheader\nmcltype matrix\ndimensions %dx%d\n)\n"%(dimension, dimension)
    currString += "(mcldoms\n"
    for index, key in enumerate(sampleList): #writes the doms string, as well as the tab file
        currString += "%d "%index                      
    currString += "$\n)\n"
        
    #The data portion
    currString += "(mclmatrix\nbegin\n"
    xcount = 0
    
    #This parts allows the matrix to divide by the first value,
    #thus normalizing all input_type to between zero and 1
    for xindex, x in enumerate(sampleList):
        currString += "%d\t"%xindex
        for yindex, y in enumerate(sampleList):
            #value is pre-normalized during matrix construction
            currString += "%s:%f "%(yindex, mclcTree[x][y])
        currString += "$\n"
    currString += ')'
    
    
    
    return currString
    

'''[[string]] -> bool
Helper: see if this pattern has singletons
'''
def hasSingleton(pattern):
    if len(pattern) == 1: return True
    for line in pattern:
        if len(line) == 1:
            return True
    return False

def notAllZeroes(matrix):
    '''checks a nested dictionary whether it's all zeroes'''
    for x in matrix.values():
        for y in x.values():
            if y > 0:
                return True
    return False
    
    
def interDistance(clusters, matrix):
    '''([[str]], {str{str:int}}) -> int
    computes the average of the average distances between all strains in each cluster to all other strains
    '''
    def add(x, y):
        return x + y
    def avg(list):
        return sum(list) / len(list)
    
    if len(clusters) <= 1:
        return float(0)
    
    if not notAllZeroes(matrix):
        return float(0)
    
    avg_dists = []
    for cluster in clusters:
        other_strains = ft.reduce(add, [x for x in clusters if x != cluster], [])
        avg_dists.append(avg([avg(list(map(lambda x: matrix[x][strain], other_strains))) for strain in cluster]))
    
    return avg(avg_dists)

def intraDistance(clusters, matrix):
    '''([[str]], {str{str:int}}) -> int
    computes the average of the average distances within each cluster
    '''
    def avg(list):
        return sum(list) / len(list)
    
    if not notAllZeroes(matrix):
        return float(0)
    
    avg_dists = []
    for cluster in clusters:
        if len(cluster) <= 1:
            continue
        clus_avgs = []
        for strain in cluster:
            other_strains = [x for x in cluster if x != strain]
            avg_dists.append(avg([avg(list(map(lambda x: matrix[x][strain], other_strains))) for strain in cluster]))
    
    return avg(avg_dists)
    
def loadMatrix(data, tab_file):
    '''(string(path), string(path)) -> {str{str:int}}
    loads the matrix from the first autogrouper
    '''
    def splitEntry(entry):
        return re.split(':', entry)
    def splitLine(data):
        return [re.split('\s', x) for x in re.split('\n', data.rstrip('\n'))]

    noheader = re.search('(?s)mclmatrix\nbegin\n(.+)\n', data).group(1)
    
    tab = loadTab(tab_file)
    
    matrix = {tab[int(line[0])] : {tab[int(x[0])] : float(x[1]) for x in map(splitEntry, line[1:-1])} for line in splitLine(noheader)}
#     #-1 because there's a semicolon at the end of each line.
#     matrix = np.array([[map(splitEntry, line[1:-1])]for line in splitLine(noheader)])
    
    return matrix

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

def loadTab(tabfile):
    '''(str)->list
    loads a tab file for mcl'''
    with open(tabfile, 'r') as tmp:
        data = tmp.read()
        
    tab = [re.match('^[0-9]+\s(.+)$', line).group(1) for line in re.split('\n', data) if line != '']
    
    return tab

#Helpers for analyze distance
def calcDistValues(currString, matrix, tab_file, i, pi_list):
    
    results = []
    for pi in pi_list:
        clusters = mcl(currString, tab_file, i, pi)
        inter_value = interDistance(clusters, matrix)
        intra_value = intraDistance(clusters, matrix)
        results.append((i, pi, inter_value, intra_value))
    
    return results
        
def calcDistValues_unpack(param_tuple):
    
    currString = param_tuple[0]
    matrix = param_tuple[1]
    tab_file = param_tuple[2]
    i = param_tuple[3]
    pi_list = param_tuple[4]
    
    return calcDistValues(currString, matrix, tab_file, i, pi_list)
# def analyzeDistance(currString, tab_file, params):
#     '''string(path), string(path) -> [num]
#     analyze the avg distance measure at a range of i values to see if there's a good one
#     '''
#     #helper functions
#        
#     def findMax(dictionary):
#         max = 0
#         for e in dictionary.values():
#             for e2 in e.values():
#                 if e2 > max:
#                     max = e2
#         return max
#          
#     def differential(values):
#         results = [(values[0][0], 0, 0)]
#         for ind, ele in enumerate(values[1:]):
#             results.append((ele[0], (ele[1] - values[ind][1]) / max(1,values[ind][1]), (ele[2] - values[ind][2]) / max(1,values[ind][2])))
#         return results
#      
#     def generateDistParams(currString, matrix, tab_file, i_list, pi_list):
#         '''generates params for the pool-based multiprocessing'''
#          
#         return [(currString, matrix, tab_file, i, pi_list) for i in i_list]
#     ###
#     start_time = time.time()     
#     matrix = loadMatrix(currString, tab_file)
#      
#     iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
#     piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
#  
#     param_list = generateDistParams(currString, matrix, tab_file, iVals, piVals)
#     pool = mp.Pool(4)
#      
#     results = pool.map(calcDistValues_unpack, param_list)
#     pool.close()
#     pool.join()
#      
#     results = reduce(lambda x, y: x + y, results) #make flat
#      
#     heat_inter_matrix = np.zeros((len(iVals), len(piVals)))
#     heat_intra_matrix = np.zeros_like(heat_inter_matrix)
#      
#     for result in results:
#         heat_inter_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[2]
#         heat_intra_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[3]
#      
# #     max = findMax(matrix)
#      
# #     heat_inter_matrix = max - heat_inter_matrix
# #     heat_intra_matrix = max - heat_intra_matrix
#      
#     result_matrix = heat_inter_matrix / heat_intra_matrix
#      
#     #gets rid of NaN values, as they don't equal themselves.
#     result_matrix[result_matrix != result_matrix] = 0
#     #gets rid of infs.
#     result_matrix[np.isinf(result_matrix)] = 0
#                                 
#     print('analyzeDist took {} seconds'.format(time.time() - start_time))
#     return result_matrix        
def analyzeDistance(currString, tab_file, params, separate=False):
    '''string(path), string(path) -> [num]
    analyze the avg distance measure at a range of i values to see if there's a good one
    '''
    #helper functions
       
    def findMax(dictionary):
        max = 0
        for e in dictionary.values():
            for e2 in e.values():
                if e2 > max:
                    max = e2
        return max
         
    def differential(values):
        results = [(values[0][0], 0, 0)]
        for ind, ele in enumerate(values[1:]):
            results.append((ele[0], (ele[1] - values[ind][1]) / max(1,values[ind][1]), (ele[2] - values[ind][2]) / max(1,values[ind][2])))
        return results
     
    def generateDistParams(currString, matrix, tab_file, i_list, pi_list):
        '''generates params for the pool-based multiprocessing'''
         
        return [(currString, matrix, tab_file, i, pi_list) for i in i_list]
    ###
    start_time = time.time()     
    matrix = loadMatrix(currString, tab_file)
     
    iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
    piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
    results = []
    
    for i in iVals:
        for pi in piVals:
            clusters = mcl(currString, tab_file, i, pi)
            inter_value = interDistance(clusters, matrix)
            intra_value = intraDistance(clusters, matrix)
            results.append((i, pi, inter_value, intra_value))
            
    heat_inter_matrix = np.zeros((len(iVals), len(piVals)))
    heat_intra_matrix = np.zeros_like(heat_inter_matrix)
     
    for result in results:
        heat_inter_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[2]
        heat_intra_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[3]
     
    max = findMax(matrix)
      
    heat_inter_matrix = 1 - heat_inter_matrix
    heat_intra_matrix = 1 - heat_intra_matrix
    
    #Standard Operation
    if not separate: 
        result_matrix = heat_inter_matrix / heat_intra_matrix
        result_matrix = heat_intra_matrix
    #gets rid of NaN values, as they don't equal themselves.
        result_matrix[result_matrix != result_matrix] = 0
    #gets rid of infs.
        result_matrix[np.isinf(result_matrix)] = 0
        return result_matrix
    else:
    #Separate Inter and Intra, for when you're not optimizing
        return heat_inter_matrix, heat_intra_matrix

                           
    print('analyzeDist took {} seconds'.format(time.time() - start_time))


    
    
# def analyzeDistance(currString, tab_file, params):
#     '''string(path), string(path) -> [num]
#     analyze the avg distance measure at a range of i values to see if there's a good one
#     '''
#     #helper functions
# #     start_time = time.time()
#     def findMax(dictionary):
#         max = 0
#         for e in dictionary.values():
#             for e2 in e.values():
#                 if e2 > max:
#                     max = e2
#         return max
#         
#     def differential(values):
#         results = [(values[0][0], 0, 0)]
#         for ind, ele in enumerate(values[1:]):
#             results.append((ele[0], (ele[1] - values[ind][1]) / max(1,values[ind][1]), (ele[2] - values[ind][2]) / max(1,values[ind][2])))
#         return results
# 
#     def generateDistParams(currString, matrix, tab_file, i_list, pi_list):
#         '''generates params for the pool-based multiprocessing'''
#         
#         return [(currString, matrix, tab_file, i, pi_list) for i in i_list]
#     ###
#         
#     matrix = loadMatrix(currString, tab_file)
#     
#     iVal = params.getIMin()
#     piVals = [params.getPiMin()]
# 
#     results = calcDistValues(currString, matrix, tab_file, iVal, piVals)
#     
# #     results = reduce(lambda x, y: x + y, results) #make flat
#     
#     heat_inter_matrix = np.zeros((1, 1))
#     heat_intra_matrix = np.zeros_like(heat_inter_matrix)
#     
#     for result in results:
#         heat_inter_matrix[0, 0] = result[2]
#         heat_intra_matrix[0, 0] = result[3]
#     
# #     max = findMax(matrix)
# #     
# #     heat_inter_matrix = max - heat_inter_matrix
# #     heat_intra_matrix = max - heat_intra_matrix
#     
#     result_matrix =  heat_intra_matrix / heat_inter_matrix
#     
#     #gets rid of NaN values, as they don't equal themselves.
#     result_matrix[result_matrix != result_matrix] = 0
#     #gets rid of infs.
#     result_matrix[np.isinf(result_matrix)] = 0
#                                
# #     print('analyzeDist iteration took {} seconds'.format(time.time() - start_time))
#     return result_matrix
    
    
#     #     diagnostics
#     import matplotlib
#     matplotlib.use('pdf')
#     import matplotlib.pyplot as plt
#     from matplotlib.backends.backend_pdf import PdfPages
#     x = [e[0] for e in results]
#     y = [e[1] for e in results]
#     z = [e[2] for e in results]
#     plot1 = plt.subplot(2,1,1)
#     plot1.set_title('InterDistance')
#     plot1.plot(x, y)
#     plot2 = plt.subplot(2,1,2)
#     plot2.set_title('IntraDistance')
#     plot2.plot(x, z, 'r')
# 
#     plt.show(block=False)
#     pdf = PdfPages('DistMeasure.pdf')
#     pdf.savefig()
#     pdf.close()  

    
    return results


def graphHeatMap(file_name, matrices, extents, names, axis_labels, title):
    '''current style only gens 2x2 heat maps! change length, width to make
    different kinds'''

    length = 2
    width = 2
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
    
def generateGraph(datapath, tabpath, params):
    
    output_name = 'Heatmaps.pdf'
    graph_title = 'S2 Clustering Metrics'
    
    with open(datapath) as f:
        data = f.read()
    
    cln_matrix, eff_matrix = analyzeClm(data, tabpath, params)
    inter_matrix, intra_matrix = analyzeDistance(data, tabpath, params, separate=True)
    
    matricies = [cln_matrix, eff_matrix, inter_matrix, intra_matrix]
    extents = [[params.getPiMin(),params.getPiMax() + params.getPiStep(),params.getIMin(),params.getIMax() + params.getIStep()]]*4
    names = ['Number of Clusters', 'Efficiency', 'Inter-Cluster Distance', 'Intra-Cluster Distance']
    axis_labels = [('pI value', 'I value')] * 4
    graphHeatMap(output_name, matricies, extents, names, axis_labels, graph_title)


#helpers for optimization
def findOptimalValue(matrix, params):
    '''for S1 and S2'''
    highest_location = np.where(matrix == np.nanmax(matrix))
    maxes = matrix.shape #the size of the whole thing
    optimal_i = params.getIMin() + params.getIStep() * (maxes[1] - highest_location[0][0] - 1)
    optimal_pi = params.getPiMin() + params.getPiStep() * highest_location[1][0]
    
    return optimal_i, optimal_pi

def loadPersistentGroups(file_path):
    '''for S3'''
    
    block_pattern = '(?s)#[0-9]+\n(.+?)\n\n'
    
    with open(file_path, 'r') as input_type:
        blocks = re.findall(block_pattern, input_type.read())
    
    temp_result = 0
    
    for block in blocks:
        lines = re.split('\n', block)
        temp_result += len(lines)
    
    return float(temp_result) / len(blocks)

#####    
# def findS1(matrix_dict, tab_path, params, section_length):
#     '''
#     Searches for the optimal I and pi values for phase 1 clustering.
#     Randomly samples 0.1% (flexible) of all sections.
#     Takes efficiency and inter/intra cluster dist into consideration.
#     '''   
#     
#     #helper functions
#     def S1Recursive(samples):
#         '''
#         the samples are just index, to be extrated from matrix_iterator
#         depends on analyzeClm and analyzeDist to return 3 matrices:
#         # of clusters, efficiency, and inter/intra distance
#         '''
#         sample = samples[0]
#         currString = SnpSorter.buildMatrix(matrix_iterator[sample], matrix_iterator[sample].keys())
#         clm_result = analyzeClm(currString, tab_path, params)
#         
#         clus_matrix = clm_result[0]
#         eff_matrix = clm_result[1]
#         
#         dist_matrix = analyzeDistance(currString, tab_path, params)
#         
#         if len(samples) == 1:
#             return clus_matrix, eff_matrix, dist_matrix
#         else:
#             prev_result = S1Recursive(samples[1:])
#             return clus_matrix + prev_result[0], eff_matrix + prev_result[1], dist_matrix + prev_result[2]
#     
#     #Attempts to run a number of sample matrices in order to estimate optimal I and pi vales
#     S1_start_time = time.time()
#      
#     matrix_iterator = []
#     for chr in matrix_dict:
#         matrix_iterator += matrix_dict[chr].values()
# 
# #Currently not sampling. Sampling procedure needs to be reworked somehow.
# #     n = len(matrix_iterator)
# #     m = n / 100
# #     samples = random.sample(range(n), m)
#     
# #     print('S1 test for i = {}, pi = {}'.format(str(params.get), str(pi)))
#     
#     samples = matrix_iterator
#     m = len(samples)
#     
#     clus_matrix, eff_matrix, dist_matrix = S1Recursive(samples)
#     
#     #get the average. the recursive func doesn't do this.
#     clus_matrix = clus_matrix / m
#     eff_matrix = eff_matrix / m
#     dist_matrix = dist_matrix / m
#     
#     
#     #Auto selection prioritize larger I and smaller pi
#     consensus_matrix = eff_matrix * dist_matrix
#     optimal_i, optimal_pi = findOptimalValue(consensus_matrix, params)
#     
#     #Graph results
#     output_folder = params.getOutputFolder()
#     graph_title = 'Phase 1 Clustering Parameter Optimization for Section Length {}'.format(section_length)
#     matricies = [clus_matrix, eff_matrix, dist_matrix]
#     extents = [[params.getPiMin(),params.getPiMax(),params.getIMin(),params.getIMax()]]*3
#     names = ['Number of Clusters', 'Efficiency', 'Inter/Intra-Cluster Distance']
# 
#     graphHeatMap(output_folder + "S1 Optimization for {}".format(section_length), matricies, extents, names, graph_title)
#     
#     print('S1 Optimization took {} seconds'.format(str(time.time() - S1_start_time)))  
#     return optimal_i, optimal_pi

def findS1(aggregateCount, tab_path, params):
    '''
    Searches for the optimal I and pi values for phase 1 clustering.
    Randomly samples 0.1% (flexible) of all sections.
    Takes efficiency and inter/intra cluster dist into consideration.
    
    This is actually way more like findS2 now. Just clusters the S2 matrix using a range of I and PI values
    and outputs the serialized matrix
    '''   
    
    #helper functions
#     def S1Recursive(matrix_iterator):
#         '''
#         the samples are just index, to be extrated from matrix_iterator
#         depends on analyzeClm and analyzeDist to return 3 matrices:
#         # of clusters, efficiency, and inter/intra distance
#         '''
# 
#         result_list = []
#         
#         while True:
#             try:
#                 matrix = matrix_iterator.next()
#             except StopIteration:
#                 break
#             
#             currString = SnpSorter.buildMatrix(matrix, matrix.keys())
#             clm_result = analyzeClm(currString, tab_path, params)
#             clus_matrix = clm_result[0]
#             eff_matrix = clm_result[1]
#             dist_matrix = analyzeDistance(currString, tab_path, params)
#             result_list.append((clus_matrix, eff_matrix, dist_matrix))
#             
#         print('S1 Resource usage: {}kb'.format(getattr(resource.getrusage(resource.RUSAGE_SELF), 'ru_maxrss') / 1000))
#             
#         if matrix_iterator.length() == 1:
#             return result_list
#         else:
#             reduced_result_list = reduce(lambda x, y: (x[0] + y[0], x[1] + y[1], x[2] + y[2]), result_list)
#             return reduced_result_list
    
#     class matrixDictIterator(object):
#         '''
#         have a generator instead of a new list to save memory
#         Not safe against empty input!
#         '''
#         def __init__(self, matrix_dict):
#             self.matrix_dict = matrix_dict
#             self.chr_iter = iter(matrix_dict.keys())
#             self.chr = self.chr_iter.next()
#             self.index = 0
#         
#         def __iter__(self):
#             return self
#         
#         def __next__(self):
#             return self.next()
#         
#         def next(self):
#             try:
#                 tmp = self.matrix_dict[self.chr][self.index]
#                 self.index += 1
#             except KeyError:
#                 try:
#                     self.chr = self.chr_iter.next()
#                     self.index = 0
#                     tmp = matrix_dict[self.chr][self.index]
#                 except StopIteration:
#                     raise StopIteration()  
#             return tmp
#                 
#         def length(self):
#             return sum([len(self.matrix_dict[x].keys()) for x in self.matrix_dict.keys()])
            
    #Attempts to run a number of sample matrices in order to estimate optimal I and pi vales
    S1_start_time = time.time()
     
    tab = loadTab(tab_path)
    currString = buildMatrix(aggregateCount, tab)
    
    dist_matrix = analyzeDistance(currString, tab_path, params)
    clus_matrix, eff_matrix = analyzeClm(currString, tab_path, params) #returns a tuple! (# of clusters, efficiency)
    
    clus_name_pattern = '{0}/I{1}PI{2}S{3}_CLUS'
    dist_name_pattern = '{0}/I{1}PI{2}S{3}_DIST'
    eff_name_pattern = '{0}/I{1}PI{2}S{3}_EFF'
    
    for m, n in zip([dist_matrix, clus_matrix, eff_matrix], [dist_name_pattern, clus_name_pattern, eff_name_pattern]):
        np.save(n.format(params.getOutputFolder(), int(params.getIVal() * 10), int(params.getPiVal() * 10), params.getSectionLength()), m)
    
                
def findS2(aggregateCount, tabfile, params, section_length):
    #finds the parameters for the second step (global) clustering
    S2_start_time = time.time()
    tab = loadTab(tabfile)
    currString = buildMatrix(aggregateCount, tab)
    
    dist_matrix = analyzeDistance(currString, tabfile, params)
    clus_matrix, eff_matrix = analyzeClm(currString, tabfile, params) #returns a tuple! (# of clusters, efficiency)
    
    clus_matrix = clus_matrix
    eff_matrix = eff_matrix
    dist_matrix = dist_matrix
    
    consensus_matrix = eff_matrix * dist_matrix
    optimal_i, optimal_pi = findOptimalValue(consensus_matrix, params)
    
    #Graph results
    output_folder = params.getOutputFolder()
    graph_title = 'Phase 2 Clustering Parameter Optimization for Section Length of {}'.format(section_length)
    matricies = [clus_matrix, eff_matrix, dist_matrix]
    extents = [[params.getPiMin(),params.getPiMax(),params.getIMin(),params.getIMax()]]*3
    names = ['Number of Clusters', 'Efficiency', 'Inter/Intra-Cluster Distance']
    
    graphHeatMap(output_folder + "S2 Optimization for {}".format(section_length), matricies, extents, names, graph_title)
    
    print(consensus_matrix)
    print('S2 Optimization took {} seconds'.format(str(time.time() - S2_start_time))) 
    return optimal_i, optimal_pi

 
def findS3(persistentGroups_path):
    #returns the load groups function results
    
    return loadPersistentGroups(persistentGroups_path)


    
    
if __name__ == "__main__":
#     workingDir = '/data/new/javi/plasmo/pipeline/matrix'
    workingDir = '/data/new/javi/yeast/pipeline/WinVar/matrix'
#     workingDir = '/data/new/javi/toxo/SNPSort20/matrix'
#     workingDir = '/data/new/javi/yeast/pipeline/matrix'
    os.chdir(workingDir)
    datapath = workingDir + '/groups.txt.mci'
    tabpath = workingDir + '/persistentMatrix.tab'
    generateGraph(datapath, tabpath, 'HeatMaps.pdf', 'Yeast-10K')
    

    print('Autogrouper Finished.')

