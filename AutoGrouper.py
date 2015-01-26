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

def group(mclcTree, tabpath, outpath):
    import MCLCounter as mclc
    sampleList = sorted(mclc.loadSampleList(tabpath))
    
    currString = buildMatrix(mclcTree, sampleList, outpath)

    #MCL   
    
    
    result = findOptimalPattern(currString, tabpath)
        
     
    with open(outpath, 'w') as output:
        try:
            for line in result:
                gname = line[0]
                output.write("@{0:s}\n{1:s}\n".format(gname, "\t".join(line)))
        except:
            print('Overall Clustering Failed!')

'''(String, list) -> [[string]]
runs mcl for the matrix over a range of value, to find the best one.

Current Scheme:
finds the largest PI/I that does not result in a singleton cluster
'''
def findOptimalPattern(currString, tabpath):
    iMax = 8
    piMax = 20
    step = 0.1
    patterns = []
    
    for x in reversed(np.arange(0, piMax, step)):
        temp = mcl(currString, tabpath, iMax, x)
        if not hasSingleton(temp): return temp
        
    for x in reversed(np.arange(0, iMax, step)):
        temp = mcl(currString, tabpath, x, 0)
        if not hasSingleton(temp): return temp
    

    
'''(dict, list, num, num) -> [[string]]
helper function for repeatedly running mcl over an array of I and PI values.
'''            
def mcl(currString, tabpath, iValue, piValue):
            
    p1 = Popen(["echo", currString], stdout = subp.PIPE)
    p2 = Popen(["mcl", "-", "-use-tab", tabpath, "-I", str(iValue), "-pi", str(piValue), "-o", "-", "-q", "x", "-V", "all"], stdin = p1.stdout, stdout = subp.PIPE, close_fds=True)
    
    results = [re.split("\t", line) for line in re.split("\n",p2.stdout.read())[:-1]]

    return results


def buildMatrix(mclcTree, sampleList, outpath):
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
    
    with open(outpath + ".mci", 'w') as moutpath:
        moutpath.write(currString)
    
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
