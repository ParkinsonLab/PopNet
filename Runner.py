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
    print("Runner Script Start")
#     #Job: Encode whole genome information into MCL format.
#     filePath = "/data/javi/Toxo/64Genomes/matrixBackup/persistentResult.txt"
#     tabPath = "/data/javi/Toxo/64Genomes/matrixBackup/persistentMatrix.tab"
#     workingDirectory = "/data/javi/Toxo/64Genomes"
#     outPath = workingDirectory + "/GenomicMCL/out.clusters"
#     MCLoutpath = workingDirectory + "/GenomicMCL/Genome.mci"
#     MCLtabpath = workingDirectory + "/GenomicMCL/Genome.tab"
#     os.chdir(workingDirectory)
#     loadedClusterTuple = MCLCounter.loadClusters(filePath, tabPath)
#     countedTuple = MCLCounter.count(loadedClusterTuple, outPath)
#     aggregatedClusters = MCLCounter.aggregate(countedTuple[0])['Genome'] #Temporary, until I standardize
#     print aggregatedClusters
#     NexusEncoder.toMCL(aggregatedClusters, MCLoutpath)
#     NexusEncoder.toMCLTab(aggregatedClusters, MCLtabpath)
    
    #Job: Similarity Comparisons
    
#     #Job: Encoding to Cytoscape File
#     filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
#     tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
#     outpath = "/data/javi/Toxo/64Genomes/Filtered/cytoscape{0}.gml"
#     
#     dataTree = mclc.loadClusters(filepath, tabpath)
#     counted = mclc.count(dataTree, "/data/javi/Toxo/64Genomes/Filtered/counted.txt")
#     for name, chr in counted[0].items():
#         ce.parse(chr, name, counted[1], outpath.format(name))

    #Job: Running the Group Composition Script
    filepath = "/data/javi/Toxo/64Genomes/Filtered/persistentResult.txt"
    tabpath = "/data/javi/Toxo/64Genomes/Filtered/persistentMatrix.tab"
    outpath = "/data/javi/Toxo/64Genomes/Filtered/GroupComposition.txt"
    dataTree = mclc.toMatrix(mclc.loadClusters(filepath, tabpath)[0])
    
    composition = gc.findComposition('ME49', dataTree)
    colored = cp.calculateColor(composition)
    cp.write(colored, outpath)
    print("Runner Script End")