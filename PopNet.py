'''
Created on Dec 15, 2014
 
@author: javi

Main Runner Class
'''
  
def sanitize_input(config, var_list, debug_var_list):
    
    #check various inputs to make sure they are fine

    if not os.path.isdir(config.get('Settings', 'base_directory')):
        raise ValueError("base directory doesn't exist")
    
    if not config.get('Settings', 'organism') in ['toxoplasma', 'yeast', 'plasmodium', 'strep']:
        raise ValueError("Unsupported organism: {}".format(config.get('Settings', 'organism')))
    
    if not config.get('Settings', 'input_type') in ['nucmer', 'tabular']:
        raise ValueError("Unsupported input_type type: {}".format(config.get('Settings', 'input_type')))
    

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
    import NodeSummary as ns
    import sys
    from os import listdir
    from os.path import isfile, join, isdir
    import re
    import time
    import ParamWrapper as pw
    import configparser


#config loading    
    var_list = ['base_directory', 'organism', 'input_type', 'file_name', 'section_length', 'S1_iVal', 'S1_piVal', 'S2_iVal', 'S2_piVal', \
                'reference', 'optimize']
    
    debug_var_list = ['optimize_level', 'section_length_min', 'section_length_max', 'section_length_step', 'S1_iVal_max', 'S1_iVal_min', 'S1_iVal_step', \
                      'S1_piVal_max', 'S1_piVal_min', 'S1_piVal_step', 'S2_iVal_max', 'S2_iVal_min', 'S2_iVal_step', 'S2_iVal_max', \
                      'S2_iVal_min', 'S2_iVal_step'] 
    
    config_file_path = sys.argv[1]
    config = configparser.SafeConfigParser({'autogroup': True})
    config.read(config_file_path)
    
    
    #Settings
    base_directory = config.get('Settings', 'base_directory')
    try:
        outputDirectory = config.get('Settings', 'output_directory')
    except:
        outputDirectory = base_directory + '/matrix/'
    
    organism = config.get('Settings', 'organism')
    #options: toxoplasma, yeast, plasmodium
    
    input_type = config.get('Settings', 'input_type')
    #options: nucmer, tabular
    
    if input_type == 'tabular':
        file_name = config.get('Settings', 'file_name') #for tabular file only
    else:
        file_name = 'Not Used.'
    #set these according to optimize info
    section_length = config.getint('Settings', 'section_length')
    
    S1_iVal = config.getfloat('Settings', 'S1_iVal')
    S1_piVal = config.getfloat('Settings', 'S1_piVal')

    S2_iVal = config.getfloat('Settings', 'S2_iVal')
    S2_piVal = config.getfloat('Settings', 'S2_piVal')
    
    reference = config.get('Settings', 'reference')
    optimize = config.getboolean('Settings', 'optimize')
    autogroup = config.getboolean('Settings', 'autogroup')

    
################################################################################################
#Debug Info
#optimize organism calculated all parameters
    if optimize == True:
        optimize_level = config.getint('Settings', 'optimize_level') #1 = only S1, 2 = S1 + S2, 3 = all 3
        section_length_min = config.getint('Optimization', 'section_length_min')
        section_length_max = config.getint('Optimization', 'section_length_max')
        section_length_step = config.getint('Optimization', 'section_length_step')
        
        S1_iVal_min = config.getfloat('Optimization', 'S1_iVal_min')
        S1_iVal_max = config.getfloat('Optimization', 'S1_iVal_max')
        S1_iVal_step = config.getfloat('Optimization', 'S1_iVal_step')
        
        S1_piVal_min = config.getfloat('Optimization', 'S1_piVal_min')
        S1_piVal_max = config.getfloat('Optimization', 'S1_piVal_max')
        S1_piVal_step = config.getfloat('Optimization', 'S1_piVal_step')
        
        S2_iVal_min = config.getfloat('Optimization', 'S2_iVal_min')
        S2_iVal_max = config.getfloat('Optimization', 'S2_iVal_max')
        S2_iVal_step = config.getfloat('Optimization', 'S2_iVal_step')
        
        S2_piVal_min = config.getfloat('Optimization', 'S2_piVal_min')
        S2_piVal_max = config.getfloat('Optimization', 'S2_piVal_max')
        S2_piVal_step = config.getfloat('Optimization', 'S2_piVal_step')
    
#Directory containing input_type
#     base_directory = '/data/new/javi/plasmo/pheno_select'
#     base_directory = '/data/new/javi/yeast/partial'
#     base_directory = '/scratch/j/jparkins/xescape/plasmo/pipeline'

    #Do Not Change
    
    cytoscapeDirectory = outputDirectory + '/cytoscape'
    rawresultpath = outputDirectory + "/results.txt"
    perresultpath = outputDirectory + "/persistentResult.txt"
    permatrixpath = outputDirectory + "/persistentMatrix.txt"
    tabpath = outputDirectory + "/persistentMatrix.tab"
    outpath = cytoscapeDirectory + "/cytoscape{0}.xgmml"
    tab_networkpath = cytoscapeDirectory + "/tabNetwork.tsv"
    matrixDirectory = cytoscapeDirectory + "/countMatrices"
    matrixoutpath = matrixDirectory + "/{0}.txt"
    densitypath = outputDirectory + "/density.txt"
    countpath = outputDirectory + "/counted.txt"  
    grouppath = outputDirectory + "/groups.txt"
    groupmcipath = grouppath + ".mci"
    colorout_path = outputDirectory + "/colors.txt"
    logpath = outputDirectory + "/log.txt"
    matrix_density_path = outputDirectory + '/matrix_density.tsv'
    S1_directory = outputDirectory + "/S1_optimization/"
    S2_directory = outputDirectory + "/S2_optimization/"   
################################################################################################
#Setup
    start_time = time.time()
    
    for folder in [outputDirectory, cytoscapeDirectory, matrixDirectory]:
        if not isdir(folder):  
            os.mkdir(folder)
            
    if optimize == True and optimize_level >= 1:
        S1params = pw.ParamWrapper()
        
        S1params.setIMax(S1_iVal_max)
        S1params.setIMin(S1_iVal_min)
        S1params.setIStep(S1_iVal_step)
        S1params.setPiMax(S1_piVal_max)
        S1params.setPiMin(S1_piVal_min)
        S1params.setPiStep(S1_piVal_step)
        S1params.setOutputFolder(S1_directory)
        if not isdir(S1_directory):
            os.mkdir(S1_directory)
    
    if optimize == True and optimize_level >= 2:
        S2params = pw.ParamWrapper()
        
        S2params.setIMax(S2_iVal_max)
        S2params.setIMin(S2_iVal_min)
        S2params.setIStep(S2_iVal_step)
        S2params.setPiMax(S2_piVal_max)
        S2params.setPiMin(S2_piVal_min)
        S2params.setPiStep(S2_piVal_step)
        S2params.setOutputFolder(S2_directory)
        if not isdir(S2_directory):
            os.mkdir(S2_directory)
    
    default_params = pw.ParamWrapper()
    
    default_params.setIMax(20)
    default_params.setIMin(2)
    default_params.setIStep(0.5)
    default_params.setPiMax(20)
    default_params.setPiMin(1)
    default_params.setPiStep(0.5)
    
            
    S1_set = False
    S2_set = False
    blength_set = False
    log = []
    sanitize_input(config, var_list, debug_var_list)
    
    log.append('Logging Parameters Used:')
    for name in sorted(var_list):
        log.append('{0}\t{1}'.format(name, str(globals()[name])))
    if optimize:
        for name in sorted(debug_var_list):
            log.append('{0}\t{1}'.format(name, str(globals()[name])))
    
    print('Config Loading Successful')
    log.append('Config Loading Successful\n')
#################################################################################################
#Input Processing
            
    if input_type == 'tabular':
        #Grigg data
        os.chdir(base_directory)
        griggpath = base_directory + '/' + file_name
        data = tl.load(griggpath, reference)
    #     excludepath = outputDirectory + '/exclude.txt'
    #     data = tl.load(griggpath, reference, excludepath)
        dataTree = data[0]
        sampleList = sorted(data[1])
        grey_dict = data[2] #see if there are any grey regions.
        grey_dict = tl.convertGreyDict(grey_dict, section_length)
        #format = {sample:[(chr, position)]}
        
    elif input_type == 'nucmer':
        #Sequence data
        #Locates all vcf and snps files in the directory and adds them to the dataTree one by one.  
    
        #Set the sample name pattern here!
        #        pattern = '^(.+?)[_].*' #for the old plasmodium stuff
        pattern = '^(.+?)[\.].*'
          
        os.chdir(base_directory) 
        try:
            onlyfiles = [ f for f in listdir(base_directory) if (isfile(join(base_directory,f)) and (f.endswith(".snps") or f.endswith(".vcf"))) ]
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
        grey_dict = {n:[] for n in sampleList}
        
    with open(rawresultpath, "w") as results:
        snps.record(dataTree, results, sampleList)
        
#################################################################################################   
#Analysis

    
            
    #Analysis
    loop = True 
    if optimize == True and optimize_level >= 3:
        section_length = section_length_min
        blength_list = []
    
    os.chdir(outputDirectory)
      
    print("calculating density")    
    snps.snpDensity(dataTree,densitypath,sampleList, section_length)
              
    import DriftDetection as dd
    dataTree = dd.scan(dataTree)
         
    print('filling data tree')
    dataTree = snps.fillDataTree(dataTree, sampleList, reference)
        
    while loop: 
        print("generating matrix")     
        if not isdir(outputDirectory):    
            os.mkdir(outputDirectory)
        matrix = snps.calculateMatrix(dataTree, sampleList, section_length)
        snps.outputMatrixDensity(matrix, matrix_density_path)
        
        
        snps.normalizeMatrix(matrix)
        snps.recordTab(sampleList, tabpath)
        
        if optimize == True and optimize_level >= 1 and S1_set == False:
            print('Step 1 Optimization')
            S1_iVal, S1_piVal = ag.findS1(matrix, tabpath, S1params, section_length)
            log.append('Step 1 Optimization for section_length of {0}:\ni = {1}\npi = {2}'.format(str(section_length), str(S1_iVal), str(S1_piVal)))
            log.append('Search range was i between {0} and {1}, pi between {2} and {3}\n'.format(str(S1_iVal_max), str(S1_iVal_min), str(S1_piVal_max), str(S1_piVal_min)))
         
        print('Cluster step 1')
        snps.recordMatrix(matrix, sampleList, tabpath, permatrixpath, perresultpath, S1_iVal, S1_piVal)
        
        
        print('encoding to nexus')     
        treetuple = (dataTree, sampleList)
        nex.nexusOutput(tl.aggregateForNexus(treetuple))
            
        reloaded_dataTree = mclc.loadClusters(perresultpath, tabpath)
        sampleList = reloaded_dataTree[1]
        dataMatrix = mclc.toMatrix(reloaded_dataTree[0])
        counted = mclc.count(reloaded_dataTree, countpath)
        aggregateCount = list(mclc.aggregate(counted[0]).values())[0]
        
        if optimize == True and optimize_level >= 2:
            print('Step 2 Optimization')
            S2_iVal, S2_piVal = ag.findS2(aggregateCount, tabpath, S2params, section_length)
            log.append('Step 2 Optimization for section_length of {0}:\ni = {1}\npi = {2}'.format(str(section_length), str(S2_iVal), str(S2_piVal)))
            log.append('Search range was i between {0} and {1}, pi between {2} and {3}\n'.format(str(S2_iVal_max), str(S2_iVal_min), str(S2_piVal_max), str(S2_piVal_min)))
        
        if autogroup == True:   
            ag.group(aggregateCount, tabpath, grouppath, groupmcipath, S2_iVal, S2_piVal)
            ag.generateGraph(groupmcipath, tabpath, default_params)
        
        #to load groups        
        groups = gc.loadGroups(grouppath, "")
        expandedGroups = gc.expandGroups(groups)
        strainList = [x[0] for x in expandedGroups]
        groupColors = cp.createColorTable({"everything" : groups.keys()}, groups.keys())
    #    for 64 genomes
    #     density = gc.loadDensity(densityPath) 
          
        #8genomes doesn't use density! The loadMultiDensity function is currently hacked, and sets everything to 100!!!!
        density = gc.loadMultiDensity(densitypath)   
                
                
        composition = gc.cytoscapeComposition(strainList, dataMatrix, grouppath, grey_dict)  
        colorTable = {}
        for strain in expandedGroups:
            colorTable[strain[0]] = groupColors[frozenset([strain[1]])]
        for group in groups.keys():
            colorTable[group] = groupColors[frozenset([group])]
                   
        #enables the use of black spacer tiles
        colorTable['SPACER'] = 0
        colorTable['GREY'] = 13421772
        print(colorTable)
        with open(colorout_path, 'w') as color_out:
            color_out.write(cp.printColorTable(colorTable))
           
    #    for the whole thing
        mclc.printMatrix(aggregateCount, matrixoutpath.format("aggregate"))
         
        aggregateComp = gc.aggregate(composition, organism)
        gc.tab_output(composition, sampleList, colorTable, section_length, tab_networkpath)
        rev_groups = gc.reverseGroups(groups)
        ce.parse(aggregateCount, "Genome", counted[1], outpath.format("Genome"), colorTable, aggregateComp, rev_groups)
        
        os.chdir(base_directory)
        if optimize == True and optimize_level >= 3:
            blength_list.append((section_length, ag.findS3(perresultpath)))   
            section_length += section_length_step
            if section_length > section_length_max:
                loop = False 
                log.append('Step 3 Optimization results:\nblength\tAverage number of clusters')
                for entry in blength_list:
                    log.append('{0}\t{1}'.format(entry[0], entry[1]))
                print(blength_list)
        else:
            loop = False
        
    with open(logpath, 'w') as log_output:
        log_output.write('\n'.join(log))
    
    print("Runner Completed")
    print('Run time was {0} seconds'.format(time.time() - start_time))