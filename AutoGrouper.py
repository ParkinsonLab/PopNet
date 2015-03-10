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
    iVal = analyzeClm(currString, tabpath)
    if iVal == 0:
        raise RuntimeError('No suitable iValue found!')
        import sys
        sys.exit()
       
    result = mcl(currString, tabpath, iVal, 1)
        
     
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
    
    if len(currString) > 100:
        tempname = 'temp.mci'
        with open(tempname, 'w') as temp:
            temp.write(currString)
        result = subp.check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '8'])
    else:
        p1 = subp.Popen(["echo", currString], stdout = subp.PIPE)
        result = subp.check_output(["mcl", "-", "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '8'], stdin = p1.stdout, close_fds=True) #The -I option is the inflation value. Play around with it. 
        p1.stdout.close()

            
#     p1 = Popen(["echo", currString], stdout = subp.PIPE)
#     p2 = Popen(["mcl", "-", "-use-tab", tabpath, "-I", str(iValue), "-pi", str(piValue), "-o", "-", "-q", "x", "-V", "all"], stdin = p1.stdout, stdout = subp.PIPE, close_fds=True)
#     
    results = [re.split("\t", line) for line in re.split("\n",result)[:-1]]
    
    if raw:
        return result
    
    return results


def clminfo(currString, tabpath):
    '''runs clminfo in with a range of values, and finds the pattern with the least amount of change over time.'''
    
    step = 0.2
    imin = 2
    imax = 8
    pimin = 1
    pimax = 20
    
    namePattern = 'out.{}'
    matrixName = 'in.mci'
    files = []
    
    
    with open(matrixName, 'w') as f:
        f.write(currString)
        
    iVal = imin
    while iVal <= imax:
        filename = namePattern.format(str(int(iVal * 10 )))
        files.append(filename)
        subp.Popen(['mcl', matrixName, '-I', str(iVal), '-o', filename, "-V", "all", '-te', '8'], close_fds=True)
        iVal += step     
    
    
#     iVal = imin
#     while iVal <= imax:
#         filename = namePattern.format(str(int(iVal * 10 )) + '10')
#         files.append(filename)
#         subp.Popen(['mcl', matrixName, '-I', str(iVal), '-o', filename, "-V", "all", '-te'. '8'], close_fds=True)
#         iVal += step
#         
    piVal = pimin
    while piVal <= pimax:
        filename = namePattern.format(str(int(80 + piVal * 10)))
        files.append(filename)
        subp.Popen(['mcl', matrixName, '-I', '8', '-pi', str(piVal), '-o', filename, "-V", "all", '-te', '8'], close_fds=True)
        piVal += step
        
    result = subp.check_output(['clm', 'info', matrixName] + files, close_fds=True)
    dresult = subp.check_output(['clm', 'dist', '-mode', 'sj', '--chain'] + files, close_fds=True)
        
    return result, dresult


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
    distpattern = 'd=(.+?)\s'
    
    data, distdata = clminfo(currString, tabpath)
    print(data, distdata)
    lines = re.split('===', data)
    dists = [re.search(distpattern, l).group(1) for l in re.split('\n', distdata.rstrip('\n'))]
    #first position is efficiency, then the # of clusters, then the iValue used
    parsedLines = [(float(re.search(effPattern, line).group(1)), int(re.search(clusPattern, line).group(1)), int(re.search(namePattern, line).group(1))) for line in lines]
    clusNums = [x[1] for x in parsedLines]
    
    lowest = parsedLines[0][1]
    highest = parsedLines[-1][1]
    
    filteredLines = [x for x in parsedLines if (highest > x[1] > lowest)]
    
    if len(filteredLines) > 0:
        sfl = sorted(filteredLines, key=lambda x: clusNums.count(x[1]), reverse=True)
        candidates = sorted(sfl[:clusNums.count(sfl[0][1])], key = lambda x: (x[0], x[1]), reverse=True)
    
        if len(candidates) > 1 and clusNums.count(candidates[0][1]) > 1:
            result = candidates[0][2] / 10
#             iresult = (candidates[0][2] // 100) / 10
#             piresult = (candidates[0][2] % 100) / 10
            cn = candidates[0][1]
        else:
            iresult = 8
#             piresult = 19
            cn = candidates[0][1]
    else:
        result = 8
#         piresult = 19
        cn = parsedLines[-1][1]
    
    print('Autogrouper: {:d} was chosen as I value, for {:d} clusters'.format(result, cn))
    
#     diagnostics
    import matplotlib
    matplotlib.use('pdf')
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    x = [e[2] for e in parsedLines]
    y = [e[1] for e in parsedLines]
    z = [e[0] for e in parsedLines]
    plt.subplot(3,1,1)
    plt.plot(x, y)
    plt.subplot(3,1,2)
    plt.plot(x, z, 'r')
    plt.subplot(3,1,3)
    plt.plot(x[:-1], dists, 'g')
    plt.show(block=False)
    pdf = PdfPages('Autogroup.pdf')
    pdf.savefig()
    pdf.close()  
       
    return result
    
    
    
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

if __name__ == "__main__":
#     workingDir = '/data/new/javi/plasmo/pipeline/matrix'
    workingDir = '/data/new/javi/toxo/SNPSort20/matrix'
    datapath = workingDir + '/groups.txt.mci'
    tabpath = 'persistentMatrix.tab'
    import os
    os.chdir(workingDir)
    with open(datapath) as f:
        data = f.read()
    print(analyzeClm(data, tabpath))

