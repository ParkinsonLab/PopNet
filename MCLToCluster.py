'''
Created on May 11, 2015

@author: javi

This script provides functions that can convert an existing MCL file or string into a
cluster file readable by Cluster3.0. The tab file is necessary for the names. 
'''

import re
import numpy as np
import os
import ChrTranslator as ct
import AutoGrouper as ag

prefixDict = {
              '': 'NON',
              'Mali': 'AF_MA',
              'BurkinaFaso': 'AF_BF',
              'Gambia': 'AF_GA',
              'Ghana': 'AF_GH',
              'Cambodia(Pursat)': 'AS_CU',
              'Cambodia(Tasanh)' : 'AS_CT',
              'Cambodia(Pailin)' : 'AS_CA',
              'Cambodia(Ratanakiri)' : 'AS_CR',
              'Thailand' : 'AS_TH',
              'Vietnam' : 'AS_VI'
              }


def toCluster(elements, sampleList):
    '''(str, list) -> str
    constructs the table used for clustering'''
        
    clusString = ''
    #top row
    clusString += '\t'.join([''] + sampleList) #empty element in the front

    for i, line in enumerate(elements):
        clusString += '\n'
        clusString += '\t'.join([sampleList[i]] + map(str, line))
    
    return clusString

def formatMCL(mclList, sampleList):
    '''([str], [str]) -> [[[num]]]
    Add together a list of mcl strings into one big matrix'''
    
    def getValue(e):
        '''get the value from one element in the mcl file like 1:xxxx'''
        return float(re.split(':', e)[1])
    
    bodyPattern = '(?s)begin\n(.+?)\n[)]'
    
    
    results = np.zeros((len(sampleList), len(sampleList)))
    for matrix in mclList:
        mclString = re.search(bodyPattern, matrix).group(1)
        elements = [map(getValue, re.split('\s', line)[1:-1]) for line in re.split('\n', mclString)]
        results += np.array(elements)
        
    return results / results[0,0]

def findMatrices(locs, data):
    '''((str, [int]), str) -> [str]
    returns the matrices at the locations specified by locs.
    The input dict should contain the chr number in the first position and
    a list corresponding to the stretch of interest'''
    
    basePattern = '(?s)@.*?_{chrName}\n.*?#{matrixID}\n(.*?\n[)]\n)(?=#|$)'
    
    chrNum, matrices = locs
    
    matrixList = []
    for matrixID in matrices:
        pattern = basePattern.format(chrName = ct.reverseTranslate(chrNum), matrixID = str(matrixID))
        matrixList.append(re.search(pattern, data).group(1))
    
    return matrixList

def loadLocationTable(filename):
    delimiter = ','
    
    with open(filename, 'r') as infile:
        data = infile.read()
    
    results = {}
    lines = re.split('\n', data)
    for line in lines:
        try:
            elements = re.split(delimiter, line)
            results[elements[0]] = elements[1]
        except:
            pass
    
    return results
    
def addInfoToOutput(output, loctable, grouptable, sampleList):
    '''(str, dict, [str]) -> str
    adds the country code to sample names according to the table'''
    def addInfo(matchobj):

        try:
            return '_'.join([prefixDict[loctable[name]],'GP{}'.format(grouptable[name]), name])
        except KeyError:
            return '_'.join([prefixDict[''], 'GP{}'.format(grouptable[name]), name])
    
    
    for name in sampleList:
        output = re.sub(name, addInfo, output)
    return output
    
        
if __name__ == '__main__':
    folder = '/data/new/javi/plasmo/pipeline/matrix/'
#     filename = 'groups.txt.mci'
    filename = 'persistentMatrix.txt'
    tabname = 'persistentMatrix.tab'
    groupsname = 'groups.txt'
    locname = 'locTable.csv'
    outname = 'section6_cluster.tsv'
    locs = ('7', range(125, 135))
    
    
    os.chdir(folder)
    with open(filename, 'r') as infile:
        data = infile.read()
    

    sampleList = ag.loadTab(tabname)

    import GroupComposition as gc
    groups = gc.loadGroups(groupsname, None)
    revgroups = gc.reverseGroups(groups)
     
    locTable = loadLocationTable(locname)
    
    #for multiple matrices  
    matrices = formatMCL(findMatrices(locs, data), sampleList)

#     #for single matrix
#     matrices = formatMCL([data], sampleList)
    
    with open(outname, 'w') as outfile:
        output = addInfoToOutput(toCluster(matrices, sampleList), locTable, revgroups, sampleList)
        outfile.write(output)

    print('Conversion script completed.')