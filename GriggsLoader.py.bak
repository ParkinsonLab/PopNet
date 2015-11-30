'''
Created on Feb 16, 2014

@author: javi
'''
import re
import copy
import os
import NexusEncoder


def load(path):
    print "loading.."
    with open(path) as source:
        nameList = filter(None, re.split("\t",source.readline().replace("\r\n", "").replace("\n", ""))) #gets names, split, filter for preceeding tabs. Name will be in order.
        rawData = re.split("\n", source.read())
        tree = {}
        chrBranch = {}
        posBranch = {}
    
            
        for line in rawData[:-1]:
            lineSplit = re.split("\t", line.replace("\n", ""))
            if lineSplit[0] not in tree:
                chrBranch = {}
                tree[lineSplit[0]] = chrBranch
            
            if int(lineSplit[1]) not in chrBranch:
                posBranch = {}
                chrBranch[int(lineSplit[1])] = posBranch
            
            for snp, samp in zip(lineSplit[3:], nameList):
#                print "%s\t%s"%(samp, snp)
                posBranch[samp] = snp.replace("\r", "")
                
    print "loading done"
    return (tree, nameList)

def aggregateForNexus(treeTuple):
    print "begin aggregate"
    results = []
    baseBlock = {}
    for name in treeTuple[1]:
        baseBlock[name] = ""
    block = copy.deepcopy(baseBlock)
    
    for chrs in treeTuple[0].items():
        for pos in sorted(chrs[1].items()):
            for name, info in pos[1].items():
                block[name] += info
                last = name
            if len(block[last]) >= 33:
                results.append(block)
                block = copy.deepcopy(baseBlock)
        results.append(block)
    print "aggregate done" 
    return {"Genome" : results}

if __name__ == '__main__':    
    os.chdir("/data/javi/Toxo/64Genomes")
    #OrderedSNPV6 
    treeTup =  load("/data/javi/Toxo/64Genomes/OrderedSNPV6.txt")  
    NexusEncoder.nexusOutput(aggregateForNexus(treeTup))
    
    with open("density.txt", "w") as densityOutput:
    	import SnpSorter
    	SnpSorter.snpDensity(treeTup[0], densityOutput, treeTup[1])
    
    print "end of GriggsLoader script"
            