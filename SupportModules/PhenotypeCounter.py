'''
Created on Jul 21, 2015

@author: javi
'''
import re
import numpy as np
import scipy as sp
from scipy.stats import ttest_ind
from matplotlib import pyplot as plt
from matplotlib import mlab
import itertools as it
from scipy.cluster import vq
import GroupComposition as gc
import functools as ft
    
def parseGroups(file):
    '''
    parses the groups file
    '''
    groups = gc.loadGroups(file, None)
    del groups['NONE']
    return groups, gc.reverseGroups(groups)

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
    
def filter(data_tree):
    '''
    filters out tests that have missing data
    '''
    
    for test, test_data in list(data_tree.items()):
        for group_data in list(test_data.values()):
            if 0 in list(group_data.values()):
                print('0 detected in {}'.format(test))
                del data_tree[test]
                break
    
    return data_tree
    
    
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
    vars = {test:0 for test in tests}
    means = {test:0 for test in tests}
    
    for test in vars:
        var = variance(data_tree[test][group].values())
        mean = np.mean(data_tree[test][group].values())
        vars[test] += var
        means[test] += mean
        
    return [(test, var, means[test]) for test, var in sorted(vars.items(), key = lambda x: x[1])]
    
def tTest(data_tree, test, group_a, group_b):
    
    threshold = 0.05
   
    samples_a = np.array(data_tree[test][group_a].values())
    samples_b = np.array(data_tree[test][group_b].values())

    t_stat, p_val = ttest_ind(samples_a, samples_b, equal_var = False)
    
    return p_val
    

def anova(data_tree, test, groups):
    pass

def findGroup(sample_info, sample):
    for group in sample_info.keys():
        if sample in sample_info[group].keys(): return group
    
    raise ValueError('sample not in a group')


def normalize(matrix):
    '''given a 2D numpy array, normalize it by column'''
    
    result = np.zeros_like(matrix)
    
    for col_num in range(matrix.shape[1]):
        head = max(matrix[:,col_num])
        tail = min(matrix[:,col_num])
        
        
        if head == 0.: head = 1
        avg = np.average(matrix[:,col_num])
        
        for ind, row in enumerate(matrix):
            result[ind, col_num] = (matrix[ind,col_num] - tail) / (avg - tail)
    
    return result


def formatToMatrix(data_tree):
    
    sample_info = list(data_tree.values())[0]
    
    tests = sorted(data_tree.keys())
    groups = sorted(sample_info.keys())
    samples = sorted(ft.reduce(lambda x, y: x + y, [list(sample_info[group].keys()) for group in groups]))
    
    result = np.zeros((len(samples), len(tests)))
    
    for s_ind, sample in enumerate(samples):
        for t_ind, test in enumerate(tests):
            result[s_ind, t_ind] = data_tree[test][findGroup(sample_info, sample)][sample]

    
    return tests, samples, result
    

def formatToCluster3(tests, samples, result, group_info, outpath):
    '''
    formats to something readable by cluster 3.0
    '''
    
    rev_groups = gc.reverseGroups(group_info)
    
    new_names = []
    for index, sample in enumerate(samples):
            new_name = '_'.join([rev_groups[sample], sample])
            new_names.append((index, new_name))
    new_names = sorted(new_names, key = lambda x: x[1])
    
    with open(outpath, 'w') as output:
        output.write('\t'.join([''] + tests))
        output.write('\n')

        for index, sample in enumerate(new_names):   
            output.write(sample[1] + '\t')
            output.write('\t'.join([str(x) for x in result[sample[0]]]))
            output.write('\n')
    
    print('Cluster3 File Completed')
            
# 
# def formatToMCL(tests, samples, result, group_info, outpath):
#     
#     
    
    
    
def findPC(results):
    
    combinations = it.combinations(range(results.shape[1]), results.shape[0] - 1)
    
    max = 0
    comb = None
    
    for combination in combinations:
        print('start')
        subset = results[:, combination]
        subset = subset[:, np.sum(subset, axis = 0) > 0]
        print(mlab.PCA(subset).fracs)
        subset_var = sum(mlab.PCA(subset).fracs[0:2])
        
        if subset_var > max:
            max = subset_var
            comb = combination
    
    print('findPC: PC has a variance of ' + subset_var)
    print('findPC: PC is composed of ' + ','.join(comb))
    
    return comb
        
def compareClusters(clusters_a, clusters_b):
    
    results = []
    
    for cluster in clusters_a:
        loc_scores = []
        if len(cluster) == 0:
            results.append(0)
        else:
            for clusterb in clusters_b:
                loc_score = 0
                for item in cluster:
                    if item in clusterb: loc_score += 1.
                    else: loc_score -= 0.5
                loc_score = loc_score / float(len(cluster))
                loc_scores.append(loc_score)
            results.append(max(loc_scores))
    
    return sum(results) / float(len(clusters_a))

def compareClustersExact(clusters_a, clusters_b):
    '''
    rewards exact matches
    '''
    
    results = []
    
    for cluster in clusters_a:
        loc_scores = []
        if len(cluster) == 0:
            results.append(0)
        else:
            for clusterb in clusters_b:
                loc_score = 0                
                if set(cluster) == set(clusterb): loc_score += 1.
                loc_scores.append(loc_score)
            results.append(max(loc_scores))
    return sum(results)

def cluster(matrix, preclusters):
    
    results = []
    
    for x in range(matrix.shape[1]):
        score = 0
        k = len(preclusters)
        subset = matrix[:,x]
        kmeans_result = vq.kmeans2(subset, subset[:k], iter = 50, minit = 'matrix', )
        
        clusters = []
        
        for y in range(k):
            clusters.append([])
            
        for i, z in enumerate(kmeans_result[1]):
            clusters[z].append(i)
        
        score = compareClusters(clusters, preclusters)
        results.append((score, x, clusters))
    
    results = sorted(results, key = lambda x: x[0], reverse = True)
    
    return results
                   
def phenoCount():
    directory = '/data/new/javi/yeast/phenotypes/'
    input_name = 'phenotypes.csv'
    groups_name = 'PhenoGroups.txt'
    output_name = 'sig_pheno.txt'
    
    input_path = directory + input_name
    groups_path = directory + groups_name
    output_path = directory + output_name
    
    group_dict, rev_groups = parseGroups(groups_path)
    data_tree = parse(input_path, group_dict)
    
#     sig_groups = significantGroups(data_tree)
    sig_groups = ['UWOP', 'PW5']
    
    randomkey = data_tree.keys()[0]
    
    with open(output_path, 'w') as output:
        for group in sig_groups:
            if len(data_tree[randomkey][group]) > 1:
                good_tests = goodTests(data_tree, group)
                output.write('@{}\n'.format(group))
                output.write('\n'.join(["\t".join([test, str(mean), str(var), str(tTest(data_tree, test, sig_groups[0], sig_groups[1]))]) for test, var, mean in good_tests]))
                output.write('\n')


def translateToIndex(clusters, samples):
    
    results = []
    
    for cluster in clusters:
        translated = []
        for item in cluster:
            translated.append(samples.index(item))
        results.append(translated)
        
    return results


def translateFromIndex(clusters, samples):
    
    results = []
    
    for cluster in clusters:
        translated = []
        for item in cluster:
            translated.append(samples[int(item)])
        results.append(translated)
        
    return results

def performPCA():
    directory = '/data/new/javi/yeast/phenotypes/'
    input_name = 'phenotypes.csv'
    groups_name = 'PhenoGroups.txt'
    output_name = 'pca.txt'
    
    input_path = directory + input_name
    groups_path = directory + groups_name
    output_path = directory + output_name
    
    group_dict, revgroup_dict = parseGroups(groups_path)
    data_tree = parse(input_path, group_dict)
    sample_info = data_tree.values()[0]
    
    tests, groups, samples, result = formatToMatrix(data_tree)
    color_list = ['red', 'blue', 'green', 'yellow', 'brown', 'black', 'purple']
    color_dict = {group: color_list[index] for index, group in enumerate(group_dict.keys())}
    sample_colors = [color_dict[findGroup(sample_info, sample)] for sample in samples]
    
    comb = findPC(result)
    pca_matrix = mlab.PCA(result[:, comb])
    
    plt.scatter(pca_matrix.Y[:,1], pca_matrix.Y[:,2], sample_colors)
    print('weights' + str(pca_matrix.Wt))
    print('Frac' + str(pca_matrix.fracs))
    plt.show()



def performCluster():
    
    directory = '/data/new/javi/yeast/phenotypes/'
    input_name = 'phenotypes.csv'
    groups_name = 'PhenoGroups2.txt'
    output_name = 'clus5.txt'
    
    input_path = directory + input_name
    groups_path = directory + groups_name
    output_path = directory + output_name
    
    group_dict, rev_groups = parseGroups(groups_path)
    data_tree = parse(input_path, rev_groups)
    
    #filtering
    data_tree = filter(data_tree)
    
    sample_info = list(data_tree.values())[0]
    
    tests, samples, result = formatToMatrix(data_tree)  
    
    result = normalize(result)
    
#     #color stuff
#     color_list = ['red', 'blue', 'green', 'yellow', 'brown', 'black', 'purple']
#     color_dict = {group: color_list[index] for index, group in enumerate(group_dict.keys())}
#     sample_colors = [color_dict[findGroup(sample_info, sample)] for sample in samples]
    
    
    
    pre_clusters = translateToIndex(group_dict.values(), samples)
     
#     No longer taking only the top 100!
#     best_tests = cluster(result, pre_clusters)[:100]
    best_tests = cluster(result, pre_clusters)
    best_test_info = [(tests[id[1]], id[0], id[2]) for id in best_tests]
     
#     printResults(best_test_info, samples)
    
    best_test_names = [x[0] for x in best_test_info]
    filtered_result = np.zeros((len(samples), len(best_tests)))
    ordered_best_test_names = []
    #filter out everything except for the best tests, then cluster using cluster3
    
    counter = 0
    for ind, test in enumerate(tests):
        if test in best_test_names:
            filtered_result[:,counter] = result[:,ind]
            counter += 1
            ordered_best_test_names.append(test)
    
    
    formatToCluster3(ordered_best_test_names, samples, filtered_result, group_dict, output_path)
    
def printResults(results, samples):
    for x in results:
        print(x[0])
        print(x[1])
        
        clusters = translateFromIndex(x[2], samples)
        for cluster in clusters:
            print(cluster)
    
        
if __name__ == '__main__':
    
    performCluster()
            
    print('PhenoCounter finished.')


