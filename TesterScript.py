'''
Created on Sep 18, 2013

@author: javi
'''
import os
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from subprocess import call
import time
import random
import csv

if __name__ == '__main__':
    
#     datatype = np.dtype([('name', tuple), ('data', list)])
#      
#     c = np.array([('name1', [1,2,3]),
#                   ('name2', [1,2,3]),
#                   ('name3', [1,2,3]),
#                   ('name4', [1,2,3]),
#                   ('name5', [1,2,3])],dtype=datatype )
#     z = [[[1, 10],
#           [1, 11],
#           [2, 10],
#           [20, 100],
#           [41, 105]],
#          [[1, 12],
#           [8, 13],
#           [2, 9],
#           [22, 101],
#           [19, 100]]]
#     a = np.array(z)
#     tolerance = 15
#     element = [20, 100]
#      
#     mask1 = np.abs(a[:,:,0] - element[0]) < tolerance
#     mask2 = np.abs(a[:,:,1] - element[1]) < tolerance
#     des = np.where(mask1 * mask2)        
#     print(mask1)
#     print(mask2)
#     print(des)

#     import ClusterSimulations as cs
#     import MCLCounter
#     filepath = r"F:\Documents\ProjectData\64Genomes\Counting\persistentResult.txt"
#     tabpath = r"F:\Documents\ProjectData\64Genomes\Counting\persistentMatrix.tab"
#     outpath = r"F:\Documents\ProjectData\64Genomes\Counting\keyblocks-simple.txt"
#      
#     clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(filepath, tabpath)[0])['@TGME49_chrXI']
#     sampleList = ['TgCat_PRC2', 'COUG']
#     print(cs.strainPairComparison(clusterTree, sampleList))
#     sampleList = ['G662M', 'COUG']
#     print(cs.strainPairComparison(clusterTree, sampleList))
#     sampleList = ['G662M', 'TgCat_PRC2']
#     print(cs.strainPairComparison(clusterTree, sampleList))
    
    
    path = "/home/javi/testzone/basic/hi.txt"
    with open(path, "w") as output:
        output.write("line1\nline2    aaaaa\n line3    bbbbb\n line4    ccccc")
    with open(path, 'r') as input:
        input.readline()
        print(input.read())





