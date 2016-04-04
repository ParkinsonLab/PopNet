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

 
if __name__ == '__main__':
    
    import SnpSorter as snps
    import TabularLoader as tl
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
    import ParamFinder as pmf
    import SupportModules.NodeSummary as ns
    import sys
    from os import listdir
    from os.path import isfile, join, isdir
    import re

#Directory containing input
#     baseDirectory = '/data/new/javi/plasmo/pheno_select'
    baseDirectory = '/data/new/javi/yeast/partial'
#     baseDirectory = '/scratch/j/jparkins/xescape/plasmo/pipeline'

    #Do Not Change
    outputDirectory = baseDirectory + '/matrix/'
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
    colorout_path = outputDirectory + "colors.txt"
    params_path = outputDirectory + "params.txt"
       
    #debug mode calculated all parameters
    debug = True
    debug_level = 2 #1 = only S1, 2 = S1 + S2, 3 = all 3
    blength_min = 2000
    blength_max = 10000
    blength_step = 1000
    
    S1iVal_min = 1
    S1iVal_max = 10
    S1iVal_step = 0.5
    
    S1piVal_min = 1
    S1piVal_max = 10
    S1piVal_step = 0.5
    
    S2iVal_min = 1
    S2iVal_max = 10
    S2iVal_step = 0.5
    
    S2piVal_min = 1
    S2piVal_max = 10
    S2piVal_step = 0.5
    
    #Settings
    mode = 'yeast' #options: toxoplasma, yeast, plasmodium
    input = 'nucmer' #options: nucmer, tabular
    filename = '' #for tabular file only
    
    #set these according to debug info
    blength = 8000
    
    S1iVal = 8
    S2piVal = 19
    
    S2iVal = 5
    S2piVal = 1.5
    
    reference = "S288C"
    graph_filename = 'YHeatMaps.pdf'
    graph_title = 'Yeast-part'
    
    

################################################################################################

    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
            
            
    if input == 'tabular':
        #Grigg data
        os.chdir(baseDirectory)
        griggpath = baseDirectory + '/' + filename
        data = tl.load(griggpath, reference)
    #     excludepath = outputDirectory + '/exclude.txt'
    #     data = tl.load(griggpath, reference, excludepath)
        dataTree = data[0]
        sampleList = sorted(data[1])
        
    elif input == 'nucmer':
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
        
    with open(rawresultpath, "w") as results:
        snps.record(dataTree, results, sampleList)
        
#####################################################################################################   
    
    S1_set = False
    S2_set = False
    blength_set = False
            
    #Analysis
    def analyze():
    
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
        
        if debug == True and debug_level >= 1 and S1_set == False:
            S1iVal, S1piVal = ag.findS1()
            S1_set = True
        
        snps.recordMatrix(matrix, sampleList)
        
        os.chdir(outputDirectory)
        print('encoding to nexus')     
        treetuple = (dataTree, sampleList)
        nex.nexusOutput(tl.aggregateForNexus(treetuple))
            
        dataTree = mclc.loadClusters(perresultpath, tabpath)
        sampleList = dataTree[1]
        dataMatrix = mclc.toMatrix(dataTree[0])
        counted = mclc.count(dataTree, countpath)
        aggregateCount = mclc.aggregate(counted[0]).values()[0]
        
        if debug == True and debug_level >= 2:
            ag.group(aggregateCount, tabpath, grouppath, groupmcipath, S2iVal, S2piVal)
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
        gc.tab_output(composition, sampleList, colorTable, blength, tab_networkpath, colorout_path)
        rev_groups = gc.reverseGroups(groups)
        ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp, rev_groups)
        
                    
        if debug == True and debug_level >= 3:
            S1eff = ag.checkS1Eff()
            features = ns.findFeatures()
            return S1eff, features
        
    for temp_blength in range(blength_min, blength_max, blength_step):
        blength = temp_blength
        analyze()
    
    print("Runner Completed")