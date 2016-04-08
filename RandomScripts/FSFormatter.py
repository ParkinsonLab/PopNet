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


def rearrangeTree(dataTree, sampleList):
    '''rearrange the tree such that it is
    ordered by sample, as opposed to by chr, and converted to numbers'''
    global convertDict
       
    newTree = {}
    for chr in dataTree:
        newTree[chr] = {}
        for sample in sampleList:
            newTree[chr][sample] = {}
            
    for chr in dataTree:
        for pos in dataTree[chr]:
            for sample in dataTree[chr][pos]:
                newTree[chr][sample][pos] = dataTree[chr][pos][sample]
    
    return newTree

def recordFS(dataTree, outpath, sampleList, organism):
    '''records structure formatted files. One line per sample with one header line'''
    suffix = '.phase'
    sampleListFileName = 'samplelist.ids'
    fileListFileName = 'filelist.txt'
    sampleList = sorted(sampleList)
    newTree = rearrangeTree(dataTree, sampleList)
    chrs = sorted(dataTree.keys(), key = lambda x: ct.translate(x, organism = organism))
    
        
    
    fileList = []
    for chr in chrs:
        file_name = chr + suffix
        fileList.append(file_name)
        headerString = '\n'.join([str(len(sampleList)), str(len(dataTree[chr]))])
        posNames = ' '.join(['P'] + sorted([str(position) for position in dataTree[chr]])) #each element is a string including all the positions in a chr
        with open('/'.join([outpath, file_name]), 'w') as output:
            output.write(headerString + '\n')
            output.write(posNames + '\n')
            for sample in sampleList:
                sampleString = ''.join([newTree[chr][sample][position] for position in sorted(dataTree[chr])])
                output.write(sampleString)
                output.write('\n')
                
    with open('/'.join([outpath, sampleListFileName]), 'w') as output:
        output.write('\n'.join(sampleList))
    
    with open('/'.join([outpath, fileListFileName]), 'w') as output:
        output.write(' '.join(fileList))



def formatToFS(directory, outpath, organism):
    '''main runner method, directory contains all necessary files
    currently supports: toxoplasma, yeast, plasmodium
    outpath is a directory where all the times are to be stored'''
    
    #load the data like in FullRunner
    if organism == 'toxoplasma':
            #Grigg data
        import GriggsLoader as gl
        file_name = 'SortedSNPs.txt'
        reference = None
        os.chdir(directory)
        griggpath = directory + '/' + file_name
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
                        dataTree = snps.addData(data, sampleName, dataTree, minCoverage, reference, organism)
                else:
                    with open(f, "r") as data:
                        dataTree = snps.addData(data, sampleName, dataTree, None, reference, organism)
                sampleList.append(sampleName)
            else:
                print("Duplicate for {0}".format(sampleName))
        sampleList = sorted(sampleList)
    
    recordFS(dataTree, outpath, sampleList, organism)
    
if __name__ == '__main__':
    directory = '/data/new/javi/toxo/structure-toxo'
    organism = 'toxoplasma'
    outpath = directory
    print('Initiating...')
    formatToFS(directory, outpath, organism)
    print('fineStructure Output Complete.')