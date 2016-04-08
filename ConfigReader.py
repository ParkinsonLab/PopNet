'''
Created on Apr 8, 2016

@author: javi

Reads config from a file! So I can put this stuff onto a cluster.
'''

import re

def clean(string):
    #strips whitespace from strings
    chars = '\t\n\s '
    return string.strip(chars)

def readConfigFile(file_path, var_list, optional_var_list):
   
#################################
#Default Values

    base_directory = '/data/new/javi/toxo/test_optimization'
    
    #Settings
    organism = 'toxoplasma' #options: toxoplasma, yeast, plasmodium
    input_type = 'tabular' #options: nucmer, tabular
    file_name = 'Toxo20.txt' #for tabular file only
    
    #set these according to optimize info
    section_length = 8000
    
    S1_iVal = 8
    S1_piVal = 19
    
    S2_iVal = 5
    S2_piVal = 1.5
    
    reference = "Me49"
    graph_filename = 'THeatMaps.pdf'
    graph_title = 'Toxo-part'
    
    
################################################################################################
#Debug Info
#optimize organism calculated all parameters
    optimize = True
    optimize_level = 3 #1 = only S1, 2 = S1 + S2, 3 = all 3
    section_length_min = 4000
    section_length_max = 10000
    section_length_step = 2000
    
    S1_iVal_min = 2
    S1_iVal_max = 8
    S1_iVal_step = 1
    
    S1_piVal_min = 2
    S1_piVal_max = 8
    S1_piVal_step = 1
    
    S2_iVal_min = 2
    S2_iVal_max = 8
    S2_iVal_step = 1
    
    S2_piVal_min = 2
    S2_piVal_max = 8
    S2_piVal_step = 1

#################################################################################################
    
    kw_pattern = "{0}[\s]=(.+?)\n"
    
    
    
    with open(file_path, 'r') as input_type:
        data = input_type.read()
    
    var_dict = {}    
    for var in var_list:
        try:
            val = clean(re.search(kw_pattern.format(var)).group(1))
            var_dict[var] = val
        except:
            ValueError('option {} is not recognized.'.format(val))
    
    for var in optional_var_list:
        try:
            val = clean(re.search(kw_pattern.format(var)).group(1))
            var_dict[var] = val
        except:
            print('warning: option {} is not recognized.'.format(val))
    
    return var_dict

if __name__ == '__main__':
    pass