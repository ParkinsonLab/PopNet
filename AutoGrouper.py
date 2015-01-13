'''
Created on Dec 16, 2014

@author: javi
'''
'''uses the whole genome MCLC data to group things automatically. by MCL'''

import subprocess as subp
from subprocess import Popen
import re

def group(mclcTree, tabpath, outpath):
    import MCLCounter as mclc
    sampleList = sorted(mclc.loadSampleList(tabpath))
    dimension = len(sampleList)
    iValue = 8
    piValue = 2
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
            
    p1 = Popen(["echo", currString], stdout = subp.PIPE)
    p2 = Popen(["mcl", "-", "-use-tab", tabpath, "-I", str(iValue), "-pi", str(piValue), "-o", "-", "-q", "x", "-V", "all"], stdin = p1.stdout, stdout = subp.PIPE, close_fds=True)
    result = p2.stdout.read()
        
    with open(outpath, 'w') as output:
        try:
            for line in re.split("\n",result)[:-1]:
                gname = re.match("([\S].+?)(?:\s|$)", line).group(1)
                output.write("@{0:s}\n{1:s}\n".format(gname, line))
        except:
            print('Overall Clustering Failed!')
    with open(outpath + ".mci", 'w') as moutpath:
        moutpath.write(currString)
            
