'''
Created on Jun 13, 2014

@author: javi
'''
import Proteins.HMMParser as hmp
import Proteins.DomainFamily as dmf
import os
import re


if __name__ == '__main__':
    #merging multiple hmm families. This script specifically deals with the file/folder structure.
    directory = "/data/javi/Proteins/results"
    resultPath = "/data/javi/Proteins/HMMResolved"
    try:
        os.mkdir(resultPath)
    except:
        pass
    
    print("reading data")
    data = {}
    for dir in os.walk(directory).next()[1]:
        family = dir
        path = "/".join([directory, dir])
        files = [ x for x in os.listdir(path) if os.path.isfile("/".join([path,x]))]
        for file in files:
            strain = re.split("\.",file)[0]
            if strain not in data:
                data[strain] = []
            data[strain].append(hmp.read("/".join([path, file]), family))
    
    print("analyzing")
    results = {}
    for strain, info in data.items():
        print("resolving strain: " + strain)
        results[strain] = dmf.choose(info)
        dmf.printTree(results[strain], "{0}/{1}.hmmr".format(resultPath, strain))
    
    print("HMM reorganiation completed.")