'''
Created on Jan 13, 2017

@author: Javi
'''
'''
purpose is to match phenotype with genotype
'''
import re 
def loadTabNetwork(path):
    '''
    the return value is a diction with the form:
    {'groups': list, 'names': list, 'data':list of split lines }
    '''
    
    with open(path, 'r') as input:
        data = input.read()
    
    lines = re.split('\n', data)[:-1]
    
#     groups = re.split('\t', lines[0])[2:]
    names = re.split('\t', lines[0])[2:]
    results = []
    for line in lines[2:]:
        items = re.split('\t', line)
        results.append(items)
        
    return {'names': names, 'data':results}

def loadPheno(path):
    '''
    This file is in CSV format
    return format is {'names':names, 'id':ids, 'data':results}
    '''
    
    with open(path, 'r') as input:
        data = input.read()
        
    lines = re.split('\n', data)[:-1]
    
    names = re.split(',', lines)[2:]
    
    results = []
    ids = []
    for line in lines[1:]:
        items = re.split(',', line)
        pheno_id = items[0:3]
        data = items[3:]
        results.append(data)
        ids.append(pheno_id)
    
    return {'names':names, 'id':ids, 'data':results}

def getPhenoPattern(table, name):
    '''get the pattern, according to some cutoff, 
    of a particular phenotype'''
    
    cutoff = []
    
    index = table['id'].index(name)
    pheno = table['data'][index]
    
    
    
    pass

def matchGenoPattern(table, pattern):
    '''given a particular pattern in the format
    of a nested list, which each member of the second
    level list should be the same'''
    def toIndex(pattern, names):
        results = []
        for group in pattern:
            temp = []
            for item in group:
                temp.append(names.index(item))
            results.append(temp)
        return results
        
    def matchline(line, pattern):
        line = line[2:]
        groups = sorted(list(set(line)))
        results = [[] for x in groups]
        ind_pattern = toIndex(pattern, names)
        for index, e in enumerate(line):
            results[groups.index(e)].append(index)
        
        for p in ind_pattern:
            tmp = None
            for g in results:
                if set(p).issubset(set(g)):
                    tmp = g
                    continue
            if tmp is not None:
                results.remove(tmp)
            else:
                return False
        return True
    
    data = table['data'] 
    names = table['names']
    results = []
    for line in data:
        match = matchline(line, pattern)
        if match:
            results.append(line[:2])
        
    return results

def outputMatches(matches):
    
    for match in matches:
        print('{0}\t{1}\n'.format(match[0], match[1]))
    
    print('{0} matches'.format(len(matches)))
    
    target = ['@SC_CHRXIII', '30000']
    if target in matches:
        print('target in match')
    else:
        print('target NOT in match')


if __name__ == '__main__':
    tab_path = '/data/new/javi/yeast/pipeline/new_heatmaps/cytoscape/tabNetwork.tsv'
    pheno_path = ''
    pheno_name = ''
    
#     pheno_table = loadPheno(pheno_path)
    geno_table = loadTabNetwork(tab_path)
    
#     pattern = getPhenoPattern(pheno_table, pheno_name)
    pattern = [['YPS128_S288C', 'YPS606_S288C', 'UWOPS83_787_3_S288C'], ['UWOPS03_461_4_S288C','UWOPS05_217_3_S288C','UWOPS05_227_2_S288C'],\
               ['K11_S288C', 'Y12_S288C', 'Y9_S288C'], ['UWOPS87_2421_S288C', 'DBVPG6040_S288C', 'YS2_S288C', 'BC187_S288C', 'L_1374_S288C'],\
               ['Y55_S288C', 'NCYC110_S288C', 'DBVPG6044_S288C','SK1_S288C']]
    
#     pattern = [['YPS128_S288C', 'YPS606_S288C', 'UWOPS83_787_3_S288C']]
    
    matches = matchGenoPattern(geno_table, pattern)
    
    outputMatches(matches)
    