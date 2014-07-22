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

if __name__ == '__main__':
#     print("Runner Script Start")
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
#     #Job: Similarity Comparisons
#      
#     #Job: Encoding to Cytoscape File
#     filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
#     tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
#     outpath = "/data/javi/Toxo/64Genomes/Filtered/cytoscape/cytoscape{0}.xgmml"
#     matrixoutpath = "/data/javi/Toxo/64Genomes/Filtered/countMatrices/{0}.txt"
#     dataTree = mclc.loadClusters(filepath, tabpath)
#     dataMatrix = mclc.toMatrix(dataTree[0])
#     counted = mclc.count(dataTree, "/data/javi/Toxo/64Genomes/Filtered/counted.txt")
#       
#     groups = gc.getGroups("")
#     expandedGroups = gc.expandGroups(groups)
#     strainList = [x[0] for x in expandedGroups]
#     groupColors = cp.createColorTable({"everything" : groups.keys()})
#     composition = gc.cytoscapeComposition(strainList, dataMatrix)
#      
#     colorTable = {}
#     for strain in expandedGroups:
#         colorTable[strain[0]] = groupColors[frozenset([strain[1]])]
#     for group in groups.keys():
#         colorTable[group] = groupColors[frozenset([group])]
#     print(colorTable)
#   
#     for name, chr in counted[0].items():
#         mclc.printMatrix(chr, matrixoutpath.format(name))
#         ce.parse(chr, name, counted[1], outpath.format(name), colorTable, composition[name])

    #Job: Running the Group Composition Script
    filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
    tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
    outpath = "/data/javi/Toxo/64Genomes/Filtered/GroupComposition.txt"
    dataTree = mclc.toMatrix(mclc.loadClusters(filepath, tabpath)[0])
#     groups = ["TgCatBr64", "TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
#                    "TgCatBr44", "MAS", "TgCatBr1"]
#     composition = gc.multiComposition(groups, dataTree)
    composition = gc.findComposition("GT1", dataTree)
    colored = cp.calculateColor(composition)
    cp.write(colored, outpath)


    print("Runner Script End")