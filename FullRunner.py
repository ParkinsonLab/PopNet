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
'''Current Mode: Plasmodium'''
 
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


    baseDirectory = '/data/new/javi/plasmo/pheno_select'

#     baseDirectory = '/scratch/j/jparkins/xescape/plasmo/pipeline'

    
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
    groupmcipath = grouppath + ".mci"
    
    #Settings
    mode = 'plasmodium'
    blength = 4000
    autogroup = True
    iVal = 6
    piVal = 2
    
    graph_filename = 'PHeatMaps.pdf'
    graph_title = 'Plasmo-part'
    
    
    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
    
#     #Grigg data
#     os.chdir(baseDirectory)
#     griggpath = baseDirectory + '/OrderedSNPV8.txt'
#     data = gl.load(griggpath, reference)
# #     excludepath = outputDirectory + '/exclude.txt'
# #     data = gl.load(griggpath, excludepath)
#     dataTree = data[0]
#     sampleList = sorted(data[1])

################################################################################################
    
    #Sequence data
    #Locates all vcf and snps files in the directory and adds them to the dataTree one by one.   
       
    #Set the sample name pattern here!
    #        pattern = '^(.+?)[_].*' #for the old plasmodium stuff
    pattern = '^(.+?)[\.].*'
      
    os.chdir(baseDirectory) 
    try:
        onlyfiles = [ f for f in listdir(baseDirectory) if (isfile(join(baseDirectory,f)) and (f.endswith(".snps") or f.endswith(".vcf"))) ]
    except Exception:
        print("\n%s is not a valid file." % f)
        sys.exit()
       
    dataTree = {}   #actually a dictionary
    sampleList = []
          
    reference = "S288C"
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
#####################################################################################################   
     
    #Analysis
    with open(rawresultpath, "w") as results:
        snps.record(dataTree, results, sampleList)
  
    print("calculating density")
    
    snps.snpDensity(dataTree,densitypath,sampleList, blength)
              
    import DriftDetection as dd
    dataTree = dd.scan(dataTree)
         
    print('filling data tree')
    dataTree = snps.fillDataTree(dataTree, sampleList, reference)
      
    print("generating matrix")     
    if not isdir(outputDirectory):    
        os.mkdir(outputDirectory)
    matrix = snps.calculateMatrix(dataTree, sampleList, blength)
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
        ag.group(aggregateCount, tabpath, grouppath, groupmcipath, iVal, piVal)
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
    print(colorTable)
       
#    for the whole thing

     
    mclc.printMatrix(aggregateCount, matrixoutpath.format("aggregate"))
     
    aggregateComp = gc.aggregate(composition, mode)
    ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp)
    print("Runner Completed")