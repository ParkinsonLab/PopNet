'''
Created on Feb 1, 2016

@author: javi

outputs a structure-compatible file. 
'''
import ChrTranslator as ct
import os
import sys
import SnpSorter as snps
import re

convertDict = {'A' : '0', 'T' : '1', 'C' : '2', 'G' : '3'}

def rearrangeTree(dataTree, sampleList):
    '''rearrange the tree such that it is
    ordered by sample, as opposed to by chr, and converted to numbers'''
    global convertDict
       
    newTree = {}
    for sample in sampleList:
        newTree[sample] = {}
        for chr in dataTree:
            newTree[sample][chr] = {}
            
    for chr in dataTree:
        for pos in dataTree[chr]:
            for sample in dataTree[chr][pos]:
                newTree[sample][chr][pos] = convertDict[dataTree[chr][pos][sample].upper()]
    
    return newTree

def record(dataTree, outpath, sampleList, mode):
    '''records structure formatted files. One line per sample with one header line'''
    newTree = rearrangeTree(dataTree, sampleList)
    
    chrs = sorted(dataTree.keys(), key = lambda x: ct.translate(x, mode = mode))
    
    posNames = ['\t'.join(sorted([chr + str(position) for position in dataTree[chr]])) for chr in chrs] #each element is a string including all the positions in a chr
    headerString = '\t'.join(posNames)
    
    placeholder = '\t'.join(['x' * 6]) #fastStructure needs something in diploid, so each line is just copied once. Also needs first 6 columns to be nothing ??? 
    
    with open(outpath, 'w') as output:
        output.write(headerString)
        output.write('\n')
    
        for sample in sorted(newTree):
            sampleStringPreformat = ['\t'.join([newTree[sample][chr][position] for position in sorted(dataTree[chr])]) for chr in chrs]
            sampleStringPreformat.insert(0, sample)
            output.write(placeholder + '\t')
            output.write('\t'.join(sampleStringPreformat))
            output.write('\n')
            output.write(placeholder + '\t')
            output.write('\t'.join(sampleStringPreformat))
            output.write('\n')



def formatToStructure(directory, outpath, mode):
    '''main runner method, directory contains all necessary files
    currently supports: toxoplasma, yeast, plasmodium'''
    
    #load the data like in FullRunner
    if mode == 'toxoplasma':
            #Grigg data
        import GriggsLoader as gl
        filename = 'SortedSNPs.txt'
        reference = None
        os.chdir(directory)
        griggpath = directory + '/' + filename
        data = gl.load(griggpath, reference)
#         excludepath = outputDirectory + '/exclude.txt'
#         data = gl.load(griggpath, reference, excludepath)
        dataTree = data[0]
        sampleList = sorted(data[1])
    else:
        pattern = '^(.+?)[\.].*'
      
        os.chdir(directory) 
        try:
            onlyfiles = [ f for f in os.listdir(directory) if (os.path.isfile(os.path.join(directory,f)) and (f.endswith(".snps") or f.endswith(".vcf"))) ]
        except Exception:
            print("\n%s is not a valid file." % f)
            sys.exit()
           
        dataTree = {}   #actually a dictionary
        sampleList = []
              
        if reference not in sampleList: sampleList.append(reference)
               
        for f in onlyfiles:
            print("\nProcessing %s ..." % f)           
            sampleName = re.match(pattern, f).group(1).upper()
            if sampleName not in sampleList: 
                if f.endswith("vcf"):
                    with open("%s_coverage.min"%(re.split("\.", f)[0]), "r") as minCoverage, open(f, "r") as data:
                        dataTree = snps.addData(data, sampleName, dataTree, minCoverage, reference, mode)
                else:
                    with open(f, "r") as data:
                        dataTree = snps.addData(data, sampleName, dataTree, None, reference, mode)
                sampleList.append(sampleName)
            else:
                print("Duplicate for {0}".format(sampleName))
        sampleList = sorted(sampleList)
    
    record(dataTree, outpath, sampleList, mode)
    
if __name__ == '__main__':
    directory = '/data/new/javi/toxo/structure-test'
    mode = 'toxoplasma'
    outpath = directory + '/structure_data.txt'
    print('Initiating...')
    formatToStructure(directory, outpath, mode)
    print('Structure Output Complete.')