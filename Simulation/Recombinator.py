'''
Created on Jul 13, 2015

@author: javi

meant to replace the whole simupop thing and just recombine manually, also to record those events.
'''

import random
import string
import re


def recombine(parent_a, parent_b, pos_tree):
    '''
    produce one offspring from two genomes, and record the event.
    assumes the format is the same for all samples: (samples,{chr:[snps]})
    '''
    parents = [parent_a, parent_b]
    counter = 0
    events = []
    
    #adjustable parameters
    rate = 0.01 / 104 #100 times less than experimental
    conversion_rate = 0.5
    conversion_mean = 100000
    conversion_sigma = 10000 #range: 70 - 130kb
    
    offspring = {chr: [] for chr in pos_tree}
    junction_positions = {chr: [] for chr in pos_tree}
    
    #names
    parenta_code = getCode(parent_a[0])
    parentb_code = getCode(parent_b[0])
    offspring_code = ''.join(random.sample(string.ascii_uppercase, 3))
    offspring_name = "_".join([parenta_code, parentb_code, offspring_code])
    
    for chr_name, chr in pos_tree.items():
        last_position = 0
        last_index = 0
        deadzone = 0
        counter = 0
        for ind, position in enumerate(chr):
            if position > last_position:
                if doesRecombine(last_position, position, rate):
                    offspring[chr_name] += parents[counter%2][1][chr_name][last_index:ind]
                    last_position = position
                    last_index = ind
                    counter += 1
                    if doesConvert(conversion_rate):
                        convert_distance = int(random.gauss(conversion_mean, conversion_sigma))
                        convert_pos = sorted(chr, cmp = lambda x, y: abs(x - last_position - convert_distance) - abs(y - last_position - convert_distance))[0]
                        convert_ind = chr.index(convert_pos)
                        offspring[chr_name] += parents[counter%2][1][chr_name][last_index:convert_ind] 
                        counter += 1
                        events.append((offspring_name, 'conv', parents[counter%2][0], parents[(counter+1)%2][0], chr_name, str(position), str(convert_pos)))

                        last_position = convert_pos
                        last_index = convert_ind                        
                    else:
                        events.append((offspring_name, 'crs', parents[counter%2][0], parents[(counter+1)%2][0], chr_name, str(position), parents[counter%2][0]))   
        offspring[chr_name] += parents[counter%2][1][chr_name][last_index:]
        
        
    return ((offspring_name, offspring), events)
    
    
def cycle(population, pos_tree, size):
    '''
    one cycle. Produces the desired number of offspring
    population structure: {strain_name: {chrs: [snps]}}
    '''
    
    #population will be a mix of clonal offspring and admixed offspring
    all_events = []
    
    while(len(population) < size):
        parent_a, parent_b = select(population.keys())
        offspring, events = recombine((parent_a, population[parent_a]), (parent_b, population[parent_b]), pos_tree)
        population[offspring[0]] = offspring[1]
        all_events += events
        
    return population, all_events


#Helpers
def select(strain_list):
    '''
    picks two strains to recombine, and record the event
    '''
    return tuple(random.sample(strain_list, 2))

def doesRecombine(first, last, rate):
    a = random.random()
    b = rate
    return a < b

def doesConvert(rate):
    return random.random() < rate

def formatName(name):
    return '_'.join(['-','-', name])

def getCode(name):
    return re.match('.+?[_].+?[_](.+?)$', name).group(1)
    
if __name__ == '__main__':
    pass