'''
Created on Dec 15, 2014
 
@author: javi
'''

'''list of data-specific changes you ought to make:
1. directory
2. reference name
4. VCFnameconverter
6. possibly the inflation values
7. possible penalty score in contigsComposition
'''
'''Current Mode: Griggs'''
 
if __name__ == '__main__':
    
    import SnpSorter as snps
    import GriggsLoader as gl
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
    import sys
    from os import listdir
    from os.path import isfile, join, isdir
    import re


    baseDirectory = '/data/new/javi/toxo/SNPSort20'
#     baseDirectory = '/scratch/j/jparkins/xescape/SNPSort'
    outputDirectory = baseDirectory + '/matrix/'
    cytoscapeDirectory = outputDirectory + '/cytoscape'
    rawresultpath = outputDirectory + "results.txt"
    perresultpath = outputDirectory + "persistentResult.txt"
    tabpath = outputDirectory + "persistentMatrix.tab"
    outpath = cytoscapeDirectory + "/cytoscape{0}.xgmml"
    matrixDirectory = cytoscapeDirectory + "/countMatrices"
    matrixoutpath = matrixDirectory + "/{0}.txt"
    densitypath = outputDirectory + "density.txt"
    countpath = outputDirectory + "counted.txt"  
    grouppath = outputDirectory + "groups.txt"
    reference = 'ME49'
    mode = 'toxo'
    
    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
    
    #Grigg data
    
    os.chdir(baseDirectory)
    griggpath = baseDirectory + '/OrderedSNPV8.txt'
    data = gl.load(griggpath, reference)
#     excludepath = outputDirectory + '/exclude.txt'
#     data = gl.load(griggpath, reference, excludepath)
    dataTree = data[0]
    sampleList = sorted(data[1])

################################################################################################
    
    #Analysis
    with open(rawresultpath, "w") as results:
        snps.record(dataTree, results, sampleList)
 
    print("calculating density")
   
    snps.snpDensity(dataTree,densitypath,sampleList)
     
    print("generating matrix")     
    if not isdir(outputDirectory):    
        os.mkdir(outputDirectory)
    matrix = snps.calculateMatrix(dataTree, sampleList)
    snps.recordMatrix(matrix, sampleList)
        
    print('encoding to nexus')     
    treetuple = (dataTree, sampleList)
    nex.nexusOutput(gl.aggregateForNexus(treetuple))
    
#         #for reanalyzing Only!
#     os.chdir(outputMatrix)
#     print("remcl")
#     snps.remcl()
#     print("Reanalyzing Matrix")
#     snps.analyzeMatrix(perresultpath)
#     #reloaded      
#     dataTree = mclc.loadClusters(perresultpath, tabpath)
    
    dataTree = mclc.loadClusters(perresultpath, tabpath)
    sampleList = dataTree[1]
    dataMatrix = mclc.toMatrix(dataTree[0])
    counted = mclc.count(dataTree, countpath)
    aggregateCount = mclc.aggregate(counted[0]).values()[0]
#    ag.group(aggregateCount, tabpath, grouppath)
    
    #to load groups        
    groups = gc.loadGroups(grouppath, "")
    expandedGroups = gc.expandGroups(groups)
    strainList = [x[0] for x in expandedGroups]
    groupColors = cp.createColorTable({"everything" : groups.keys()}, groups.keys())
#    for 64 genomes
#     density = gc.loadDensity(densityPath) 
      
    #8genomes doesn't use density! The loadMultiDensity function is currently hacked, and sets everything to 100!!!!
    density = gc.loadMultiDensity(densitypath)   
            
            
    composition = gc.cytoscapeComposition(strainList, dataMatrix, density, grouppath)  
    colorTable = {}
    for strain in expandedGroups:
        colorTable[strain[0]] = groupColors[frozenset([strain[1]])]
    for group in groups.keys():
        colorTable[group] = groupColors[frozenset([group])]
               
    #enables the use of black spacer tiles
    colorTable['SPACER'] = 0
    print(colorTable)
       
#    for the whole thing

     
    mclc.printMatrix(aggregateCount, matrixoutpath.format("aggregate"))
     
    aggregateComp = gc.aggregate(composition)
    ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp)
    print("Runner Completed")
