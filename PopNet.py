'''
Created on Dec 15, 2014
 
@author: javi

Main Runner Class
'''
#module imports
import os
import sys
import re
import time
import configparser
import logging
from pathlib import Path

#package imports
from PrimaryCluster import primaryCluster
import SecondaryCluster as sc

from LoadInput import loadToPandas
import Exporter as exporter

#import CytoscapeEncoder as exporter
import ChromosomePainter as gc

import ParamWrapper as pw
import IOTools as io 
import AnalysisTools as at 
  
def sanitize_input(config):
    
    #check various inputs to make sure they are fine

    if not os.path.isdir(config.get('Settings', 'base_directory')):
        raise ValueError("base directory doesn't exist")

    
def configLogger(_path):
    
    logger = logging.getLogger()
    logger.setLevel(logging.INFO)
    
    formatter = logging.Formatter('%(asctime)s:\n%(message)s\n')
    
    fh = logging.FileHandler(_path)
    fh.setLevel(logging.INFO)
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    
    sh = logging.StreamHandler()
    sh.setLevel(logging.INFO)
    sh.setFormatter(formatter)
    logger.addHandler(sh)
    
    return logger

def main(config_file_path):    
#config loading    
    var_list = ['base_directory', 'organism', 'input_type', 'file_name', 'section_length', 'S1_iVal', 'S1_piVal', 'S2_iVal', 'S2_piVal', \
                'reference', 'optimize']
    
    config_file_path = sys.argv[1]
    config = configparser.SafeConfigParser()
    config.read(config_file_path)
    
    #setup
    start_time = time.time()

    #Settings
    try:
        output_directory = Path(config.get('Settings', 'output_directory'))
        input_path = Path(config.get('Settings', 'input_path'))
        prefix = input_path.parts[-1].split('.')[0]

        #set these according to config info
        s1_params = pw.ParamWrapper()
        s2_params = pw.ParamWrapper()

        s1_params.setSectionLength(config.getint('Settings', 'section_length')) 
        s1_params.setIVal(config.getfloat('Settings', 'S1_iVal'))
        s1_params.setPiVal(config.getfloat('Settings', 'S1_piVal'))
        
        s2_params.setIVal(config.getfloat('Settings', 'S2_iVal'))
        s2_params.setPiVal(config.getfloat('Settings', 'S2_piVal'))

        #set stuff for autogroup
        s2_params.setIMax(10)
        s2_params.setIMin(2)
        s2_params.setIStep(0.5)
        s2_params.setPiMax(10)
        s2_params.setPiMin(1)
        s2_params.setPiStep(0.5)
        
        reference = config.get('Settings', 'reference')
        autogroup = config.getboolean('Settings', 'autogroup')

    except:
        raise RuntimeError('Error reading configureation file')

    

    #output paths 
    if not output_directory.is_dir():
        output_directory.mkdir()
    os.chdir(output_directory)

    cytoscape_path = Path("{0}.xgmml".format(prefix))
    json_path = Path("{0}.json".format(prefix))


    tab_network_path = Path("chromosome_paintings.tsv")
    matrixout_path = Path("overall_similarity.tsv")
    heatmaps_path = Path("heatmaps.pdf")

    density_path = Path("density.txt")
    group_path = Path("groups.txt")
    tab_path = Path("tab.txt")

    colorout_path = Path("colors.txt")
    log_path = Path("log.txt")
    
    #other variables
    # accepted_values = ['A', 'T', 'C', 'G', 'N']  
    logger = configLogger(log_path)
    
    #let's log some params used later
    # msg = '\n'.join(['{0}\t{1}'.format(name, str(globals()[name])) for name in sorted(var_list)])
    # logger.log(msg)

    logger.info('config loaded')

#Input Processing

    #tabular data from GTAK loaded to pandas
    hdf_path = input_path.parent / '{0}.h5'.format(prefix)
    df, sample_list = loadToPandas(hdf_path, input_path, reference)
    os.chdir(output_directory)
        

#Analysis
        
    logger.info('Start Primary Clustering')
    io.writeTab(sample_list, tab_path) 
    clusters, chr_breaks = primaryCluster(df, sample_list, s1_params, logger)

    # diagnostic
    chr_names = list(df.index.get_level_values(0))
    io.writePrimaryClusters(chr_names, chr_breaks, clusters, Path('pclusters.txt'))

    matrices = at.clustersToMatrix(clusters, sample_list)   
    overall_matrix = at.overallMatrix(matrices)

    logger.info('Start Secondary Clustering')
    group_names, overall_clusters = sc.group(overall_matrix, tab_path, group_path, s2_params)
            
    color_table = at.createColorTable(group_names, overall_clusters, sample_list)
    color_table.to_csv(colorout_path)
  
    logger.info('calculating composition')
    condensed_matrices = gc.condenseToGroupMatrix(matrices, group_names, overall_clusters, sample_list)        
    composition = gc.getChromosomePaintings(condensed_matrices, chr_breaks, overall_clusters, group_names, sample_list)
       
#    for the whole thing
    logger.info('writing output')
    # io.writeTabularPainting(composition, tab_network_path)
    # io.writeOverallMatrix(overall_matrix, sample_list)
    exporter.parse(overall_matrix, color_table, composition, group_names, overall_clusters, sample_list, prefix)
    print("PopNet Completed")
    print('Run time was {0} seconds'.format(time.time() - start_time))
    
if __name__ == '__main__':
    main(sys.argv[1])