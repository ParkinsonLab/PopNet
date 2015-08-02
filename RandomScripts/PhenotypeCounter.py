'''
Created on Jul 21, 2015

@author: javi
'''
import re
import numpy as np

def parseGroups(file):
    '''
    parses the groups file
    '''
    
    import GroupComposition as gc
    return gc.reverseGroups(gc.loadGroups(file, None))

def parse(file, groups):
    '''
    parses the input tsv into a matrix
    {test:Group:Sample:value}
    '''
    
    with open(file, 'r') as input:
        data = input.read()
    
    lines = re.split('\n', data)
    
    tests = re.split('\t', lines[3])[3:]
    result = {test:{} for test in tests}
    
    for line in lines[4:-1]:
        linesplit = re.split('\t', line)
       
        sample = linesplit[0]
        group = groups[sample]
        for test, value in zip(tests, linesplit[3:]):
            if group not in result[test]: result[test][group] = {}
            if value == '': value = 0
            result[test][group][sample] = float(value)
    
    return result
    
def variance(values):
    
    values = np.array(values)
    return np.var(values)


def significantGroups(data_tree):
    '''
    finds groups that relate together in a 
    significant number of tests
    '''
    
    results = {group: 0 for group in data_tree.values()[0].keys()}
    
    for test, groups in data_tree.items():
        for group, samples in groups.items():
            results[group] += variance(samples.values())
    
    return sorted(results.keys(), key = lambda x: results[x])
                
    
    
def goodTests(data_tree, group):
    '''
    finds tests that relate well with groups
    '''
    
    tests = data_tree.keys()
    results = {test:0 for test in tests}
    
    for test in results:
        var = variance(data_tree[test][group].values())
        results[test] += var
    
    return [(test, value) for test, value in sorted(results.items(), key = lambda x: x[1], reverse = True)][:10]
    
    

if __name__ == '__main__':
    directory = '/data/new/javi/yeast/phenotypes/'
    input_name = 'phenotypes.csv'
    groups_name = 'groups.txt'
    output_name = 'sig_pheno.txt'
    
    input_path = directory + input_name
    groups_path = directory + groups_name
    output_path = directory + output_name
    
    groups = parseGroups(groups_path)
    data_tree = parse(input_path, groups)
    
    sig_groups = significantGroups(data_tree)
    
    
    randomkey = data_tree.keys()[0]
    
    with open(output_path, 'w') as output:
        for group in sig_groups:
            if len(data_tree[randomkey][group]) > 1:
                good_tests = goodTests(data_tree, group)
                output.write('@{}\n'.format(group))
                output.write('\n'.join(["\t".join([test, str(value)]) for test, value in good_tests]))
                output.write('\n')
            
    print('PhenoCounter finished.')


