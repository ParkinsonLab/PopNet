'''
Created on Nov 9, 2016

@author: javi

This module converts a multisample .fsa file to the tabular format. Usable mostly for wayne's thing.
'''

import re

def parse(raw_text):
    data_tree = {}
    
    sep = '>'
    sections = re.split(sep, raw_text)[1:]
    length = 0
    
    for section in sections[:10]:
        split = re.split('\n', section)
        name = split[0]
        seq = split[1]
        
        if len(seq) != length:
            print('length is {}'.format(len(seq)))
            length = len(seq)
        data_tree[name] = seq.upper()
    
    return data_tree, length

def toTabular(data_tree, length, outpath):
    
    chr = 'STREP_CHRI'
    pos = 0
    sample_list = sorted(list(data_tree.keys()))
    with open(outpath, 'w') as output:
        
        header = '\t'.join(['#CHROM', 'POS'] + sample_list)
        output.write(header)
        output.write('\n')
        
        for n in range(length):
            line = '\t'.join([chr, str(n)] + [data_tree[x][n] for x in sample_list])
            output.write(line)
            output.write('\n')
        
        
if __name__ == '__main__':
    
    directory = '/data/new/javi/Strep'
    input_name = '236SPN_align.fas.txt'
    output_name = '236SPN_tabular.txt'
    
    with open('/'.join([directory, input_name]), 'r') as input:
        data_tree, length = parse(input.read())
    
    toTabular(data_tree, length, '/'.join([directory, output_name]))
    
    print('FSA Conversion completed')