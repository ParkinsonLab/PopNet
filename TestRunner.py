'''
Created on Apr 4, 2016

@author: javi

Partial runner for isolating parts for development
'''


import SnpSorter as snps
import TabularLoader as tl
import NexusEncoder as nex
import MCLCounter as mclc
import Chr_NexusEncoder as chrnex
import os
import similarity as sim
import CytoscapeEncoder as ce
import GroupComposition as gc
import ClusterPattern as cp
import ParentFinder as pf
import AutoGrouper as ag
import SupportModules.NodeSummary as ns
import sys
from os import listdir
from os.path import isfile, join, isdir
import re
import time

start_time = time.time()



#Directory containing input
#     baseDirectory = '/data/new/javi/plasmo/pheno_select'
baseDirectory = '/data/new/javi/toxo/Test'
#     baseDirectory = '/scratch/j/jparkins/xescape/plasmo/pipeline'

#Do Not Change
outputDirectory = baseDirectory + '/matrix/'
cytoscapeDirectory = outputDirectory + '/cytoscape'
rawresultpath = outputDirectory + "results.txt"
perresultpath = outputDirectory + "persistentResult.txt"
tabpath = outputDirectory + "persistentMatrix.tab"
outpath = cytoscapeDirectory + "/cytoscape{0}.xgmml"
tab_networkpath = cytoscapeDirectory + "/tabNetwork.tsv"
matrixDirectory = cytoscapeDirectory + "/countMatrices"
matrixoutpath = matrixDirectory + "/{0}.txt"
densitypath = outputDirectory + "density.txt"
countpath = outputDirectory + "counted.txt"  
grouppath = outputDirectory + "groups.txt"
groupmcipath = grouppath + ".mci"
colorout_path = outputDirectory + "colors.txt"
params_path = outputDirectory + "params.txt"

#Settings
mode = 'toxoplasma' #options: toxoplasma, yeast, plasmodium
input = 'tabular' #options: nucmer, tabular
filename = 'Toxo20.txt[test]' #for tabular file only

#debug mode calculated all parameters
debug = True
debug_level = 2 #1 = only S1, 2 = S1 + S2, 3 = all 3
blength_min = 2000
blength_max = 10000
blength_step = 1000

S1iVal_min = 1
S1iVal_max = 10
S1iVal_step = 2

S1piVal_min = 1
S1piVal_max = 10
S1piVal_step = 2

S2iVal_min = 1
S2iVal_max = 10
S2iVal_step = 0.5

S2piVal_min = 1
S2piVal_max = 10
S2piVal_step = 0.5
    
#set these according to debug info
blength = 8000

S1iVal = 8
S2piVal = 19

S2iVal = 5
S2piVal = 1.5

reference = "ME49"
graph_filename = 'YHeatMaps.pdf'
graph_title = 'Yeast-part'

os.chdir(baseDirectory)
griggpath = baseDirectory + '/' + filename
data = tl.load(griggpath, reference)
#     excludepath = outputDirectory + '/exclude.txt'
#     data = tl.load(griggpath, reference, excludepath)
dataTree = data[0]
sampleList = sorted(data[1])

print('filling data tree')
dataTree = snps.fillDataTree(dataTree, sampleList, reference)
  
print("generating matrix")     
if not isdir(outputDirectory):    
    os.mkdir(outputDirectory)
matrix = snps.calculateMatrix(dataTree, sampleList, blength)


S1iVal, S1piVal = ag.findS1(matrix, S1iVal_max, S1iVal_min, S1iVal_step, S1piVal_max, S1piVal_min, S1piVal_step)

print('run time is {} seconds'.format(time.time() - start_time))




