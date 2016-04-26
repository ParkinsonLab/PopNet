'''
Created on Apr 14, 2016

@author: javi
'''
import AutoGrouper as ag
import numpy as np
import re
import os

def readFile(file_path):
    '''
    reads one file and return the average of the matrix
    '''
    
    matrix = np.load(file_path)
    return np.average(matrix)
    

def parseName(file_path):
    '''
    return the I, Pi, and Type of data for one matrix
    Takes complete files paths
    returns (section_length, i, pi, type)
    '''
    name_pattern = 'I([0-9]+?)PI([0-9]+?)S([0-9]+?)_(.+?)[.]'
    
    file_name = file_path.split('/')[-1]
    re_result = re.search(name_pattern, file_name)
    i = float(re_result.group(1)) / 10
    pi = float(re_result.group(2)) / 10
    sect_length = int(re_result.group(3))
    type = re_result.group(4)
    
    return i, pi, sect_length, type

def searchForFiles(directory):
    '''
    reads through the directory looking for numpy objects.
    returns a list of Complete paths to those files
    '''
    
    file_list = ['{0}/{1}'.format(directory, file) for file in os.listdir(directory) if file.endswith('.npy')]
    
    return file_list


def readOptiResults(input_directory, output_directory):
    '''
    main method for getting some optimization results
    '''
    type_list = ['CLUS', 'EFF', 'DIST']
    results_list = []
    i_list = set()
    pi_list = set()
    section_length_list = set()
    numpy_files = searchForFiles(input_directory)
        
    for file in numpy_files:
        i, pi, section_length, type = parseName(file)
        avg = readFile(file)
        results_list.append((section_length, i, pi, type, avg))
        i_list.add(i)
        pi_list.add(pi)
        section_length_list.add(section_length)
        
    i_list = sorted(i_list, reverse=True)
    pi_list = sorted(pi_list)
    section_length_list = sorted(section_length_list)
    
    if not len(i_list) * len(pi_list) * len(type_list) * len(section_length_list) == len(numpy_files):
        raise ValueError('Incorrect Number of Files')
    
    matrix_dict = {section_length: {type: np.zeros((len(i_list), len(pi_list))) for type in type_list} for section_length in section_length_list}

    for result in results_list:
        matrix_dict[result[0]][result[3]][i_list.index(result[1]), pi_list.index(result[2])] = result[4]
        
    #graphing
    title_dict = {'CLUS': 'Average Number of Clusters Formed', 'EFF':'Average Clustering Efficiency', 'DIST':'Average Inter/Intra Cluster Distance'}
    for sect_length in matrix_dict:
        file_name = '{0}/SL{1}.pdf'.format(output_directory, sect_length)
        sub_titles = [title_dict[type] for type in type_list]
        matrices = [matrix_dict[sect_length][type] for type in type_list]
        extents = [[min(pi_list), max(pi_list), min(i_list), max(i_list)] for type in type_list]
        graph_title = 'Optimization for Section Length {}'.format(sect_length)
        ag.graphHeatMap(file_name, matrices, extents, sub_titles, graph_title)

    
if __name__ == '__main__':
    input_directory = '/data/new/javi/yeast/matrices'
    output_directory = '/data/new/javi/yeast/opti_graphs2'
    
    readOptiResults(input_directory, output_directory)
    print('Graphs Completed!')