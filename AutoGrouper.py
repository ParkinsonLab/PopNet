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

iMin = 2
iMax = 8
piMin = 0
piMax = 10


istep = 0.2
pistep = 0.5


def group(mclcTree, tabpath, outpath, mcipath, iVal, piVal):
    import MCLCounter as mclc
    sampleList = sorted(mclc.loadSampleList(tabpath))
    
    currString = buildMatrix(mclcTree, sampleList)

#     #MCL
#     iVal = analyzeClm(currString, tabpath)
#     if iVal == 0:
#         raise RuntimeError('No suitable iValue found!')
#         import sys
#         sys.exit()

    with open(mcipath, 'w') as moutpath:
        moutpath.write(currString)
           
    result = mcl(currString, tabpath, iVal, piVal)
        
     
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
        return subp.check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '8'])
    
    def mcl_nopi():
        return subp.check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '8'])
    
    def mcl_notab_pi():
        return  subp.check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '8'])
    
    def mcl_notab_nopi():
        return subp.check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '8'])
    
    tempname = 'temp.mci'
    with open(tempname, 'w') as temp:
        temp.write(currString)
    
    if tabpath is None and piValue > 0:
        result = mcl_notab_pi()
    elif tabpath is None:
        result = mcl_notab_nopi()
    elif piValue > 0:
        result = mcl_pi()
    else:
        result = mcl_nopi()

    if raw:
        return result
    
    results = [re.split("\t", line) for line in re.split("\n",result.rstrip('\n'))]
    return results


def clminfo(currString, tabpath):
    '''runs clminfo in with a range of values, and finds the pattern with the least amount of change over time.'''
    
    
    namePattern = 'out.{}'
    matrixName = 'in.mci'
    files = []
    
    with open(matrixName, 'w') as tmp:
        tmp.write(currString)
        
    iVals = stepList(iMin, iMax, istep, reverse=True)
    piVals = stepList(piMin, piMax, pistep)
    
    for iVal in iVals:
        for piVal in piVals:
            filename = namePattern.format('{i:0>2d}{pi:0>3d}'.format(i=int(iVal * 10), pi=int(piVal * 10)))
            files.append(filename)
            with open(filename, 'w') as tmp:
                tmp.write(mcl(currString, None, iVal, piVal, raw=True))

        
    result = subp.check_output(['clm', 'info', matrixName] + files, close_fds=True)
    
    #returns the clm dist (variance of information) but not using that atm.
#    dresult = subp.check_output(['clm', 'dist', '-mode', 'sj', '--chain'] + files, close_fds=True)
    
    for file in files:
        os.remove(file)
        
    return result


def analyzeClm(currString, tabpath):
    
    '''returns the optimal I value from 2 - 8, and pI value from 1 - 20 based on the stability, then based on efficiency.
    The maximum values are avoided in favor of any stable clustering patterns in the middle, but if it can't be helped then go for
    the highest possibility (I = 8, pi = 20). If it has to choose, it will overcluster rather than under. 
    
    (Actually not gonna do the pi thing)
    
    the out files are formatted as out.XXYY, where XX = I value * 10, and YY = piValue * 10. When the pi Value is not needed, it will default to 1.
    It seems that pi = 1 is like it's not there.'''
    
    effPattern = 'efficiency=(.+?)\s'
    clusPattern = 'clusters=(.+?)\s'
    namePattern= 'source=out\.(.+?)\s'
    
    data = clminfo(currString, tabpath)
    lines = re.split('===', data)
    #first position is efficiency, then the # of clusters, then the iValue, then the piValue used
    parsedLines = [(float(re.search(effPattern, line).group(1)), int(re.search(clusPattern, line).group(1)), float(re.search(namePattern, line).group(1)[:-3]) / 10, float(re.search(namePattern, line).group(1)[-3:]) / 10) for line in lines]
    clusNums = [x[1] for x in parsedLines]
    
    
    iVals = stepList(iMin,iMax,istep, reverse=True)
    piVals = stepList(piMin,piMax,pistep)
    heat_cluster_matrix = np.zeros((len(iVals), len(piVals)))
    heat_eff_matrix = np.zeros_like(heat_cluster_matrix)
    
    for line in parsedLines:
        x = iVals.index(line[2])
        y = piVals.index(line[3])
        heat_cluster_matrix[x, y] = line[1]
        heat_eff_matrix[x, y] = line[0]
    
    return heat_cluster_matrix, heat_eff_matrix
    
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
        arr = np.array([x.values() for x in mclcTree.values()])
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
    #thus normalizing all input to between zero and 1
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


def interDistance(clusters, matrix):
    '''([[str]], {str{str:int}}) -> int
    computes the average of the average distances between all strains in each cluster to all other strains
    '''
    def add(x, y):
        return x + y
    def avg(list):
        return sum(list) / len(list)
    
    if len(clusters) == 1:
        return float('nan')
    
    avg_dists = []
    for cluster in clusters:
        other_strains = reduce(add, [x for x in clusters if x != cluster], [])
        avg_dists.append(avg([avg(map(lambda x: matrix[x][strain], other_strains)) for strain in cluster]))
    
    return avg(avg_dists)

def intraDistance(clusters, matrix):
    '''([[str]], {str{str:int}}) -> int
    computes the average of the average distances within each cluster
    '''
    def avg(list):
        return sum(list) / len(list)
    
    avg_dists = []
    for cluster in clusters:
        if len(cluster) == 1:
            continue
        clus_avgs = []
        for strain in cluster:
            other_strains = [x for x in cluster if x != strain]
            clus_avgs.append(avg(map(lambda x: matrix[x][strain], other_strains)))
        avg_dists.append(avg(clus_avgs))
    
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
    
    with open(tab_file, 'r') as tmp:
        data = tmp.read()
    tab = [re.match('^[0-9]+\s(.+)$', line).group(1) for line in re.split('\n', data) if line != '']
    
    matrix = {tab[int(line[0])] : {tab[int(x[0])] : float(x[1]) for x in map(splitEntry, line[1:-1])} for line in splitLine(noheader)}
#     #-1 because there's a semicolon at the end of each line.
#     matrix = np.array([[map(splitEntry, line[1:-1])]for line in splitLine(noheader)])
    
    return matrix

def stepList(start, end, step, reverse=False):
    list = []
    curr = start
    while curr <= end:
        list.append(round(curr, 1))
        curr += step
    
    if reverse:
        return [x for x in reversed(list)]
    else:
        return list
    
def analyzeDistance(currString, tab_file):
    '''string(path), string(path) -> [num]
    analyze the avg distance measure at a range of i values to see if there's a good one
    '''
    
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

        
    matrix = loadMatrix(currString, tab_file)
    
    results = []
    iVals = stepList(iMin, iMax, istep, reverse=True)
    piVals = stepList(piMin, piMax, pistep)
    for i in iVals:
        for pi in piVals:
            clusters = mcl(currString, tab_file, i, pi)
            inter_value = interDistance(clusters, matrix)
            intra_value = intraDistance(clusters, matrix)
            #results order: i, pi, inter, intra
            results.append((i, pi, inter_value, intra_value))
    
    heat_inter_matrix = np.zeros((len(iVals), len(piVals)))
    heat_intra_matrix = np.zeros_like(heat_inter_matrix)
    for result in results:
        heat_inter_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[2]
        heat_intra_matrix[iVals.index(result[0]), piVals.index(result[1])] = result[3]
    
    max = findMax(matrix)
    
    heat_inter_matrix = max - heat_inter_matrix
    heat_intra_matrix = max - heat_intra_matrix
    
    
    return heat_inter_matrix, heat_intra_matrix

    
    
    
    
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


def graphHeatMap(filename, matricies, extents, names, title):
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    plt.summer()
    plt.subplots_adjust(wspace = 0.3, hspace = 0.3)
    plt.suptitle(title)
    length = len(matricies) / 2
    width = 2
    axes = []
    
    for matrix, extent, name, ind in zip(matricies, extents, names, range(len(matricies))):
        ax = plt.subplot2grid((width, length), (ind // 2, ind % 2))
        ax.set_title(name)
        img = ax.imshow(matrix, extent = extent, interpolation = 'none')
        plt.colorbar(img)
    
    
    pdf = PdfPages(filename)
    pdf.savefig()
    pdf.close()
    
def generateGraph(datapath, tabpath, output_name, graph_title):
    
    with open(datapath) as f:
        data = f.read()
    
    cln_matrix, eff_matrix = analyzeClm(data, tabpath)
    inter_matrix, intra_matrix = analyzeDistance(data, tabpath)
    
    matricies = [cln_matrix, eff_matrix, inter_matrix, intra_matrix]
    extents = [[piMin,piMax,iMin,iMax]]*4
    names = ['Number of Clusters', 'Efficiency', 'Inter-cluster Distance', 'Intra-cluster Distance']

    graphHeatMap(output_name, matricies, extents, names, graph_title)
        
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

