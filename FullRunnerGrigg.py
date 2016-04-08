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
    import SupportModules.ColorSummary as cs

#     base_directory = '/data/new/javi/toxo/WinVar'
#     base_directory = '/scratch/j/jparkins/xescape/SNPSort'
    base_directory = '/data/new/javi/toxo/SNPSort20'
    
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
    color_outpath = outputDirectory + "colors.txt"
    
    #Settings
    reference = None
    file_name = 'OrderedSNPV8.txt'
#     file_name = 'SimulatedToxo.txt'
    organism = 'toxoplasma'
    section_length = 10000
    autogroup = True
    S2_iVal = 4
    S2_piVal = 1.5

    graph_filename = 'HeatMaps.pdf'
    graph_title = 'Toxo-10K'
    
    
    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
    
    #Grigg data
    os.chdir(base_directory)
    griggpath = base_directory + '/' + file_name
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
   
    snps.snpDensity(dataTree,densitypath,sampleList, section_length)
     
    print("generating matrix")     
    if not isdir(outputDirectory):    
        os.mkdir(outputDirectory)
    matrix = snps.calculateMatrix(dataTree, sampleList, section_length)
    snps.recordMatrix(matrix, sampleList)
    
    os.chdir(outputDirectory)
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
    if autogroup:
            ag.group(aggregateCount, tabpath, grouppath, groupmcipath, S2_iVal, S2_piVal)
            ag.generateGraph(groupmcipath, tabpath, graph_filename, graph_title)
        
    #to load groups        
    groups = gc.loadGroups(grouppath, "")
    expandedGroups = gc.expandGroups(groups)
    strainList = [x[0] for x in expandedGroups]
    groupColors = cp.createColorTable({"everything" : groups.keys()}, groups.keys())
#    for 64 genomes
#     density = gc.loadDensity(densityPath) 
      
    #8genomes doesn't use density! The loadMultiDensity function is currently hacked, and sets everything to 100!!!!
    density = gc.loadMultiDensity(densitypath)   
            
            
    composition = gc.cytoscapeComposition(strainList, dataMatrix, grouppath)  
    colorTable = {}
    for strain in expandedGroups:
        colorTable[strain[0]] = groupColors[frozenset([strain[1]])]
    for group in groups.keys():
        colorTable[group] = groupColors[frozenset([group])]
               
    #enables the use of black spacer tiles
    colorTable['SPACER'] = 0
    cs.outputColors(groups, colorTable, color_outpath)
       
#    for the whole thing

     
    mclc.printMatrix(aggregateCount, matrixoutpath.format("aggregate"))
     
    aggregateComp = gc.aggregate(composition, organism)
    gc.tab_output(composition, sampleList, colorTable, section_length, tab_networkpath)
    rev_groups = gc.reverseGroups(groups)
    ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp, rev_groups)
    print("Runner Completed")
