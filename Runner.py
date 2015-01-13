'''
Created on Feb 26, 2014

@author: javi
'''

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

if __name__ == '__main__':
    print("Runner Script Start")
#     #Job: Encode whole genome information into MCL format.
#     filePath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
#     tabPath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
#     workingDirectory = "/data/javi/Toxo/64Genomes"
#     outPath = workingDirectory + "/GenomicMCL/out.clusters"
#     MCLoutpath = workingDirectory + "/GenomicMCL/Genome.mci"
#     MCLtabpath = workingDirectory + "/GenomicMCL/Genome.tab"
#     os.chdir(workingDirectory)
#     loadedClusterTuple = mclc.loadClusters(filePath, tabPath)
#     countedTuple = mclc.count(loadedClusterTuple, outPath)
#     aggregatedClusters = mclc.aggregate(countedTuple[0])['Genome'] #Temporary, until I standardize
#     print(aggregatedClusters)
#     nex.toMCL(aggregatedClusters, MCLoutpath)
#     nex.toMCLTab(aggregatedClusters, MCLtabpath)
#      
#      #Job: Similarity Cytoscape
#     #yeast
#     filepath = "/data/new/javi/yeast/results/matrix/persistentResult.txt"
#     tabpath = "/data/new/javi/yeast/results/matrix/persistentMatrix.tab"
#     outpath = "/data/new/javi/yeast/results/matrix/cytoscape/fixed_simcytoscape{0}.xgmml"
#     matrixoutpath = "/data/new/javi/yeast/results/matrix/countMatrices/fixed_sim{0}.txt"
#     densityPath = "/data/new/javi/yeast/results/matrix/density.txt"
#     countPath = "/data/new/javi/yeast/results/matrix/counted.txt"
#     simpath = "/data/new/javi/yeast/results/matrix/sim"  
#      
#        
#     groupColors = {}
#     groupColors['S288C'] = 255
#     groupColors['AWR'] = 65280
#     groupColors['UC5'] = 267386880
#     groupColors['NONE'] = 0
#        
# #     #8 genomes. Groups and Colorlist have also been modified!
# #     filepath = "/data/javi/Toxo/BaseData/tempz/matrix/persistentResult.txt"
# #     tabpath = "/data/javi/Toxo/BaseData/tempz/matrix/persistentMatrix.tab"
# #     simpath = "/data/javi/Toxo/BaseData/tempz/matrix/simFiles"
# #     outpath = "/data/javi/Toxo/BaseData/tempz/matrix/simcytoscape/cytoscape{0}.xgmml"
# #     matrixoutpath = "/data/javi/Toxo/BaseData/tempz/matrix/simcountMatrices/{0}.txt"
# #     densityPath = "/data/javi/Toxo/BaseData/tempz/matrix/density.txt"
# #     countPath = "/data/javi/Toxo/BaseData/tempz/matrix/counted.txt"
# #      groupColors = {}
# #     groupColors['ME49'] = 255
# #     groupColors['VEG'] = 65280
# #     groupColors['TYPEX'] = 0
# #     groupColors['GT1'] = 267386880
# #     groupColors['NONE'] = 0
#        
#     dataTree = mclc.loadClusters(filepath, tabpath)
#     dataMatrix = mclc.toMatrix(dataTree[0])
#     counted = mclc.count(dataTree, countPath)
#     
# #    for 64 genomes
# #     density = gc.loadDensity(densityPath) 
#      
#     #8genomes doesn't use density! The loadMultiDensity function is currently hacked, and sets everything to 100!!!!
#     groups = gc.getGroups("")
#     expandedGroups = gc.expandGroups(groups)
#     strainList = [x[0] for x in expandedGroups]
#     density = gc.loadMultiDensity(densityPath)
#              
#        
#         
#     colorTable = {}
#     for strain in expandedGroups:
#         colorTable[strain[0]] = groupColors[strain[1]]
#     
#     #enables the use of black spacer tiles
#     simColors = sim.loadToColors(simpath)
#         
#         
# #    for the whole thing
#     aggregateCount = mclc.aggregate(counted[0]).values()[0]
#     aggregateComp = simColors
#     ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp)
#       
      
#     #Job: plasmodium to cytoscape
#     
#     #Job: Encoding to Cytoscape File.
#     #     #hybrid


#     baseDirectory = '/data/javi/Toxo/64Genomes/Filtered'
#     outputDirectory = baseDirectory + '/matrix'
#     cytoscapeDirectory = outputDirectory + '/cytoscape'
#     
#     filepath = outputDirectory + "persistentResult.txt"
#     tabpath = outputDirectory + "persistentMatrix.tab"
#     outpath = cytoscapeDirectory + "/cytoscape{0}.xgmml"
#     matrixoutpath = cytoscapeDirectory + "/countMatrices/{0}.txt"
#     densitypath = outputDirectory + "density.txt"
#     countpath = outputDirectory + "counted.txt"  
#     grouppath = outputDirectory + "groups.txt"  

   
#          
#     dataTree = mclc.loadClusters(filepath, tabpath)
#     dataMatrix = mclc.toMatrix(dataTree[0])
#     counted = mclc.count(dataTree, countpath)
#       
# #     #to load groups        
# #     groups = gc.getGroups("")
#     groups = gc.loadGroups(grouppath, "")
#     expandedGroups = gc.expandGroups(groups)
#     strainList = [x[0] for x in expandedGroups]
#     groupColors = cp.createColorTable({"everything" : groups.keys()}, groups.keys())
# #    for 64 genomes
# #     density = gc.loadDensity(densityPath) 
#      
#     #8genomes doesn't use density! The loadMultiDensity function is currently hacked, and sets everything to 100!!!!
#     density = gc.loadMultiDensity(densitypath)   
#            
#            
#     composition = gc.cytoscapeComposition(strainList, dataMatrix, density, grouppath)  
#     colorTable = {}
#     for strain in expandedGroups:
#         colorTable[strain[0]] = groupColors[frozenset([strain[1]])]
#     for group in groups.keys():
#         colorTable[group] = groupColors[frozenset([group])]
#               
#     #enables the use of black spacer tiles
#     colorTable['SPACER'] = 0
#     print(colorTable)
# ##    for by chr
# #     for name, chr in counted[0].items():
# #         mclc.printMatrix(chr, matrixoutpath.format(name))
# #         ce.parse(chr, name, counted[1], outpath.format(name), colorTable, composition[name])
#       
# #    for the whole thing
#     aggregateCount = mclc.aggregate(counted[0]).values()[0]
#     
#     mclc.printMatrix(aggregateCount, matrixoutpath.format("aggregate"))
#     
#     aggregateComp = gc.aggregate(composition)
#     ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp)
  
 
 
 
#     #Job: Running the Group Composition Script
#     filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
#     tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
#     outpath = "/data/javi/Toxo/64Genomes/Filtered/GroupComposition.txt"
#     dataTree = mclc.toMatrix(mclc.loadClusters(filepath, tabpath)[0])
# #     groups = ["TgCatBr64", "TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
# #                    "TgCatBr44", "MAS", "TgCatBr1"]
# #     composition = gc.multiComposition(groups, dataTree)
#     composition = gc.findComposition("GT1", dataTree)
#     colored = cp.calculateColor(composition)
#     cp.write(colored, outpath)




#         #Job: findind parents
#     filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
# #     filepath = "/home/javi/testzone/Griggs Stuff/persistentResults.txt"
#     tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
#     outpath = "/data/javi/Toxo/64Genomes/Filtered/parents.csv"
# 
#     dataTree = mclc.loadClusters(filepath, tabpath)
#     dataMatrix = mclc.toMatrix(dataTree[0])
#     counted = mclc.count(dataTree, "/data/javi/Toxo/64Genomes/Filtered/counted.txt")
#     aggregateCountMatrix = mclc.aggregate(counted[0])
#     
#     groups = gc.getGroups("")
#     expandedGroups = gc.expandGroups(groups)
#     strainList = [x[0] for x in expandedGroups]
#     
#     pairListMatrix = pf.findParents(strainList, dataMatrix, aggregateCountMatrix)
#     pf.printParents(pairListMatrix, outpath)
    
    #Job: Single strain group composition view in from 64 Strains to Javaviewer.
    baseDirectory = '/data/javi/Toxo/64Genomes/Filtered'
    outputDirectory = baseDirectory + '/matrix'
    filepath = outputDirectory + "persistentResult.txt"
    tabpath = outputDirectory + "persistentMatrix.tab"
    outpath = outputDirectory + "composition.sim"
    densitypath = outputDirectory + "density.txt"
    countpath = outputDirectory + "counted.txt"  
    grouppath = outputDirectory + "groups.txt"  
    
    
    strain = 'ME49'
    dataTree = mclc.toMatrix(mclc.loadClusters(filepath, tabpath)[0])
    density = gc.loadMultiDensity(densitypath)
    groups = gc.loadGroups(grouppath, '')
    expandedGroups = gc.expandGroups(groups)
    strainList = [x[0] for x in expandedGroups]
    
    colorTable = {}
    groupColors = cp.createColorTable({"everything" : groups.keys()}, groups.keys())
    for group in groups.keys():
        colorTable[group] = groupColors[frozenset([group])]

    composition = gc.findContigsComposition(strain, dataTree, density, groups)
    for chrName, chr in composition.items():
        templist = []
        for position in chr:
            templist.append(groupColors[position])
        composition[chrName] = templist
    
    cp.write(composition, outpath)
    
    print("Runner Script End")