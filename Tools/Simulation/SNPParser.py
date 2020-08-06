'''
Created on Jul 6, 2015

@author: javi

Intended for the parsing of SNPs (The results.txt file) to get some ancestral genomes
for the simulated data
'''

import MCLCounter as mclc
import re
import ChrNameSorter as cns

def parse(file_path):
    '''
    parses the input file into a data tree-like structure
    for ease of picking out a few samples, also the
    positions
    
    return (data_tree, pos_tree)
    '''
    print('Parsing Data')
    
    with open(file_path, 'r') as input:
        data = input.read()
    

    #splitting the text    
    data = data.split('\n')    
    data = [line.split('\t') for line in data]
    tab = sorted(data[0][2:])
        
    #Reconstructing the SNP tree
    chr_list = set([line[0] for line in data[1:-2]])
    data_tree = {sample:{chr:[] for chr in chr_list} for sample in tab}
    pos_tree = {chr:[] for chr in chr_list}
    
    #filling out the tree
    for line in data[1:-2]:       
        chr = line[0]
        pos = int(line[1])
        pos_tree[chr].append(pos)
        for i, e in enumerate(line[2:]):
            data_tree[tab[i]][chr].append(translate(e))
    
    print('Parsing Complete.')
    return (data_tree, pos_tree)

def select(data_complex, num):
    '''
    selects a couple of samples from the whole thing.
    might add some conditions later.
    
    data complex is what the parse function returns,
    namely (data_tree, pos_tree)
    '''
    print('Randomly Selecting Ancestors...')
    import random
    data_tree = data_complex[0]
    pos_tree = data_complex[1]
    
    sample_count = len(data_tree) - 1
    sample_list = sorted(data_tree.keys())
    
#     ind_list = sorted(random.sample(range(sample_count), num))
    ind_list = [sample_list.index('ME49'), sample_list.index('VEG'), sample_list.index('GT1')]
    
    results = {sample_list[i]: data_tree[sample_list[i]] for i in ind_list}
    
    print('done.')
    return (results, pos_tree)

def restructure(data_complex, mode):
    '''
    changes the data from the simulation output to
    the Grigg format
    '''
    
    data_tree = data_complex[0]
    pos_tree = data_complex[1]
    sample_list = data_tree.keys()
    result = {chr: {position: {} for position in pos_tree[chr]} for chr in pos_tree.keys()}
    
    for chr in pos_tree.keys():
        for ind, position in enumerate(pos_tree[chr]):
            for sample in sample_list:
                if mode == 'simupop':
                    result[chr][position][sample] = translate(data_tree[sample][chr][ind])
                elif mode == 'generated':
                    result[chr][position][sample] = data_tree[sample][chr][ind]
                else: raise TypeError('Unrecognized Mode for restructure()')
    return result

def outputGriggFormat(data_tree, file_path):
    '''
    prints a file that looks like what grigg gave us
    '''
    sample_list = sorted(list(list(list(data_tree.values())[0].values())[0].keys()))
    
    with open(file_path, 'w') as output:
        
        #first line
        output.write('\t'.join(['#CHROM', 'POS'] + sample_list) + '\n')
        
        for chr_name, chr in sorted(data_tree.items()):
            for pos_id, position in sorted(chr.items()):
                output.write('\t'.join([chr_name, str(pos_id)] + [position[sample] for sample in sample_list]) + '\n' )
    
def translate(input):
    '''
    translate the alleles between numbers and the characters (ATCG)
    '''
    
    table = {'A': 0, 'T': 1, 'C' : 2, 'G' : 3, 0 : 'A', 1 : 'T', 2 : 'C', 3 : 'G'}
    
    try:
        return table[input]
    except:
        raise ValueError('Invalid Nucleotide')



if __name__ == '__main__':
    pass
    
    
    
    
    

