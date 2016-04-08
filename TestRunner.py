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
import ParamWrapper as pw

if __name__ == '__main__':
    start_time = time.time()
    
    
    
    #Directory containing input_type
    #     base_directory = '/data/new/javi/plasmo/pheno_select'
    base_directory = '/data/new/javi/toxo/Test'
    #     base_directory = '/scratch/j/jparkins/xescape/plasmo/pipeline'
    
    #Do Not Change
    outputDirectory = base_directory + '/matrix/'
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
    
    S1_directory = outputDirectory + "/S1_optimization/"
    S2_directory = outputDirectory + "/S2_optimization/"
    
    #Settings
    organism = 'toxoplasma' #options: toxoplasma, yeast, plasmodium
    input_type = 'tabular' #options: nucmer, tabular
    file_name = 'Toxo20.txt' #for tabular file only
    
    #optimize organism calculated all parameters
    optimize = True
    optimize_level = 2 #1 = only S1, 2 = S1 + S2, 3 = all 3
    section_length_min = 2000
    section_length_max = 10000
    section_length_step = 1000
    
    S1_iVal_min = 2
    S1_iVal_max = 6
    S1_iVal_step = 2
    
    S1_piVal_min = 2
    S1_piVal_max = 6
    S1_piVal_step = 2
    
    S2_iVal_min = 2
    S2_iVal_max = 6
    S2_iVal_step = 2
    
    S2_piVal_min = 2
    S2_piVal_max = 6
    S2_piVal_step = 2
        
    #set these according to optimize info
    section_length = 8000
    
    S1_iVal = 8
    S2_piVal = 19
    
    S2_iVal = 5
    S2_piVal = 1.5
    
    reference = "ME49"
    graph_filename = 'YHeatMaps.pdf'
    graph_title = 'Yeast-part'
    
    #######
    if optimize == True and optimize_level >= 1:
        S1params = pw.ParamWrapper()
        
        S1params.setIMax(S1_iVal_max)
        S1params.setIMin(S1_iVal_min)
        S1params.setIStep(S1_iVal_step)
        S1params.setPiMax(S1_piVal_max)
        S1params.setPiMin(S1_piVal_min)
        S1params.setPiStep(S1_piVal_step)
        S1params.setOutputFolder(S1_directory)
        if not isdir(S1_directory):
            os.mkdir(S1_directory)
    
    if optimize == True and optimize_level >= 2:
        S2params = pw.ParamWrapper()
        
        S2params.setIMax(S2_iVal_max)
        S2params.setIMin(S2_iVal_min)
        S2params.setIStep(S2_iVal_step)
        S2params.setPiMax(S2_piVal_max)
        S2params.setPiMin(S2_piVal_min)
        S2params.setPiStep(S2_piVal_step)
        S2params.setOutputFolder(S2_directory)
        if not isdir(S2_directory):
            os.mkdir(S2_directory)
    
    #######
    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
    
            
    os.chdir(base_directory)
    griggpath = base_directory + '/' + file_name
    data = tl.load(griggpath, reference)
    #     excludepath = outputDirectory + '/exclude.txt'
    #     data = tl.load(griggpath, reference, excludepath)
    dataTree = data[0]
    sampleList = sorted(data[1])
    
    print('filling data tree')
    dataTree = snps.fillDataTree(dataTree, sampleList, reference)
    
    
    loop = True 
    if optimize == True and optimize_level >= 3:
        section_length = section_length_min
        bin_list = []
        
        
    while loop:
        print("generating matrix")     
        if not isdir(outputDirectory):    
            os.mkdir(outputDirectory)
        matrix = snps.calculateMatrix(dataTree, sampleList, section_length)
        
        start_time = time.time()
        S1_iVal, S1_piVal = ag.findS1(matrix, tabpath, S1params, section_length)
        print('S1 optimization took {} seconds'.format(str(time.time() - start_time)))
        reloaded_dataTree = mclc.loadClusters(perresultpath, tabpath)
        sampleList = reloaded_dataTree[1]
        dataMatrix = mclc.toMatrix(reloaded_dataTree[0])
        counted = mclc.count(reloaded_dataTree, countpath)
        aggregateCount = mclc.aggregate(counted[0]).values()[0]

#         i, pi = ag.findS2(aggregateCount, tabpath, S2params)
#         print('i val is' + str(i[0]))
#         print('pi val is' + str(pi[0]))
    
    
        if optimize == True and optimize_level >= 3:
            bin_list.append((section_length, ag.findS3(perresultpath)))
            section_length += section_length_step
            if section_length > section_length_max:
                loop = False
                print(bin_list)
        else:
            loop = False

    
    
    print('run time is {} seconds'.format(time.time() - start_time))




