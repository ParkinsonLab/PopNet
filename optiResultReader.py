'''
Created on Apr 14, 2016

@author: javi
'''
import AutoGrouper as ag
import numpy as np
import re
import os
import matplotlib as mpl
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas


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


def generateColorRange(x_list, **kwargs):
    '''
    sort the lists before use
    '''
    
    import colorsys

    try:
        step = kwargs['step']
    except:
        step = 1
    
    x_range = range(min(x_list), max(x_list)+ 1, step)
    length = len(x_range)
    
    try:
        raw = kwargs['raw']
        RGB_arr = np.zeros((length,), dtype=object)
    except:
        raw = False
        RGB_arr = np.zeros((length,), dtype='|S7')
    
#     HSV_arr = np.zeros(len(x_list), len(y_list), dtype=(('Hue', '>f4'), ('Sat', '>f4'), ('Lumi', '>f4')))
     
    for x_ind, x in enumerate(x_range):
        hue = float(x_ind) / length * 0.9 #avoid going back to the reds.
        rgb = colorsys.hsv_to_rgb(hue, 1, 0.5)
        if raw:
            RGB_arr[x_ind] = rgb
        else:
            RGB_arr[x_ind] = '#{:02X}{:02X}{:02X}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    
    return x_range, RGB_arr

def generateColors(x_list, y_list, **kwargs):
    '''
    sort the lists before use
    '''
    
    import colorsys
    
    if len(y_list) == 1:
        twod = False
    else:
        twod = True

    try:
        raw = kwargs['raw']
        RGB_arr = np.zeros((len(x_list), len(y_list)), dtype=object)
    except:
        raw = False
        RGB_arr = np.zeros((len(x_list), len(y_list)), dtype='|U7')
        
#     HSV_arr = np.zeros(len(x_list), len(y_list), dtype=(('Hue', '>f4'), ('Sat', '>f4'), ('Lumi', '>f4')))
    for x_ind, x in enumerate(x_list):
        hue = float(x_ind) / len(x_list) * 0.9 #avoid going back to the reds.
        for y_ind, y in enumerate(y_list):
            if len(y_list) <= 1:
                sat = 1
            else:
                sat = float(y_ind) / len(y_list) * 0.8 + 0.2
            rgb = colorsys.hsv_to_rgb(hue, sat, 0.5)
            if raw:
                RGB_arr[x_ind, y_ind] = rgb
            else:
                RGB_arr[x_ind, y_ind] = '#{:02X}{:02X}{:02X}'.format(int(rgb[0] * 255), int(rgb[1] * 255), int(rgb[2] * 255))
    
    if twod == False:
        return x_list, RGB_arr[:,0]
    else:
        return x_list, RGB_arr
            
def makeMultiLineGraph(coords_list, labels_list, title, subtitle_list, file_name):
    '''
    coords needs [[(x, y, label)]]
    '''
    
    if len(coords_list) > 1:
        width = 2
        length = len(coords_list) // 2 + len(coords_list) % 2
    else:
        width = 1
        length = 1
    
    fig = mpl.figure.Figure(figsize=(15,15))
    canvas = FigureCanvas(fig)
    fig.subplots_adjust(wspace = 0.8, hspace = 1.5)
    fig.suptitle(title)
    lines = []
    line_labels = []
    for index, (coords, labels, subtitle) in enumerate(zip(coords_list, labels_list, subtitle_list)):
        ax = fig.add_subplot(width,length,index + 1)
        ax.set_title(subtitle, size=28)
        ax.set_xlabel(labels[0], size=24)
        ax.set_ylabel(labels[1], size =24)
        for coord in coords:
            line, = ax.plot(coord[0], coord[1], label=coord[2], color=coord[3], linewidth=2)
#             ax.text(coord[0][-1] + 0.1, coord[1][-1], coord[2], size=2)
            if index == 0:
                lines.append(line)
                line_labels.append(coord[2])
#             ax.set_xticklabels([int(x) for x in ax.get_xticks()], rotation=45, ha='right')
        ax.tick_params(axis='both', which='both', labelsize=18)
        
    fig.legend(lines, line_labels, loc='lower right')
            
    pdf = PdfPages(file_name)
    pdf.savefig(fig)
    pdf.close()

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
    title_dict = {'CLUS': 'Number of Clusters', 'EFF':'Clustering Efficiency', 'DIST':'Inter/Intra Distance'}
    labels_dict = {'CLUS': 'Average Number of Clusters Formed', 'EFF':'Average Clustering Efficiency', 'DIST':'Average Inter/Intra Cluster Distance'}
    
    
    lines_list = []
    subtitle_list = []
    labels_list = []
    line_index_list = sorted([0, len(i_list) - 1, len(i_list) / 2])
#     line_index_list = range(len(i_list))
#     color_matrix = generateColors(line_index_list, section_length_list)
    color_matrix = generateColors(section_length_list, line_index_list)
    for type in type_list:
        line_sub_list = []
        subtitle = title_dict[type]
#         for sl in section_length_list:
        for sl in [10000]:
            for ind in line_index_list:
#             for ind in [7]:
                x = pi_list
                y = matrix_dict[sl][type][ind]
                if type == 'DIST':
                    y = 1. / matrix_dict[sl][type][ind]
                label = 'I{0}SL{1}'.format(int(i_list[ind]), sl)
                color = color_matrix[1][line_index_list.index(ind), section_length_list.index(sl)]
#                 color = color_matrix[1][section_length_list.index(sl),line_index_list.index(ind)]
                line_sub_list.append((x, y, label, color))
        lines_list.append(line_sub_list)
        subtitle_list.append(subtitle)
        labels_list.append(('pi Value', labels_dict[type]))
        
    file_name = '{}/Opti.pdf'.format(output_directory)
    graph_title = 'Optimization'
    
    makeMultiLineGraph(lines_list, labels_list, graph_title, subtitle_list, file_name)
#     for sect_length in matrix_dict:
#         file_name = '{0}/SL{1}.pdf'.format(output_directory, sect_length)
#         sub_titles = [title_dict[type] for type in type_list]
#         matrices = [matrix_dict[sect_length][type] for type in type_list]
#         extents = [[min(pi_list), max(pi_list), min(i_list), max(i_list)] for type in type_list]
#         graph_title = 'Optimization for Section Length {}'.format(sect_length)
#         ag.graphHeatMap(file_name, matrices, extents, sub_titles, graph_title)

    
if __name__ == '__main__':
    input_directory = '/data/new/javi/yeast/matrices'
    output_directory = '/data/new/javi/yeast/opti_graphs2'
    
    readOptiResults(input_directory, output_directory)
    print('Graphs Completed!')