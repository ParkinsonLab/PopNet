'''
Created on Apr 28, 2014

@author: javi
'''
'''in fulfillment of John's request for a chromosome coloring scheme where each section
is colored according to which group is based on.. MCL. zzz.
Use the createColorTable and calculateColor from ClusterPattern to visualize.'''

'''Strain -> Dict
number represents the group ID, as specified in this function'''
import numpy as np
import multiprocessing as mp
import pandas
import string
import re
import copy

groups = {}

def condenseToGroupMatrix(dfs, group_names, overall_clusters, sample_list):
    '''
    take the output from at.clustersToMatrix and the results of overall clustering
    returns a bunch of matrices that are NxM where N is the n_samples and M is
    n_groups. 
    implements the logic of determining whether a sample belong to a group
    within one block
    '''
    pool = mp.Pool(None, condenseInit, (group_names, overall_clusters, sample_list))
    res = pool.map(condenseWorker, dfs)

    return res

def condenseInit(_group_names, _overall_clusters, _sample_list):
    '''
    passes some variables
    '''
    global group_names
    global overall_clusters
    global sample_list
    global cutoffs

    group_names = _group_names
    overall_clusters = _overall_clusters
    sample_list = _sample_list

    membership_cutoff = 0.3
    cutoffs = [np.ceil(membership_cutoff * len(x)) for x in overall_clusters]

def condenseWorker(input_df):
    '''
    takes a cluster and change the columns to membership to a group.
    that bit of logic is here.
    '''
    def getSelfGroup(sample):
        for i, g in enumerate(overall_clusters):
            if sample in g:
                return i

    def getMatchingGroup(sample):
        '''
        logic for determining which groups a sample can belong to
        we're sticking with at least 30 percent, but that's actually done
        in the parent function. 
        '''

        return np.array([np.sum(input_df.loc[sample, group]) for group in overall_clusters]) > cutoffs

    global group_names 
    global overall_clusters 
    global sample_list 
    global cutoffs
 
    base_df = pandas.DataFrame(np.zeros((len(sample_list), len(group_names))), index = sample_list, columns = group_names)

    for sample in sample_list:
        tmp = getMatchingGroup(sample)
        if sum(tmp) <= 0:
            tmp = np.array([0 for x in range(len(group_names))])
            tmp[getSelfGroup(sample)] = 1
        base_df.loc[sample] = tmp

    return base_df


def getChromosomePaintings(matrices, breaks, overall_clusters, group_names, sample_list):
    '''
    [df/ndarray], [[string]], [string] > ???
    entry methods to getting the chr paintings
    we now have the condensed matrices
    '''
    #TODO
    #Separate into worker processes
    MAX_SCORE = 40
    PENALTY = -8

    pool = mp.Pool(1, paintWorkerInit, (matrices, breaks, overall_clusters, group_names, sample_list, MAX_SCORE, PENALTY))
    res = pool.map(paintWorker, sample_list)

    return res

def paintWorkerInit(_matrices, _breaks, _overall_clusters, _group_names, _sample_list, _max_score, _penalty):
    global matrices
    global breaks
    global overall_clusters
    global group_names 
    global sample_list 
    global max_score 
    global penalty

    matrices = _matrices
    breaks = _breaks 
    overall_clusters = _overall_clusters 
    group_names = _group_names
    sample_list = _sample_list  
    max_score = _max_score  
    penalty = _penalty

def paintWorker(sample):
    '''
    paints one sample
    '''

    def getSelfGroup(s):
        '''
        finds which group this sample came from
        '''
        for i, group in enumerate(overall_clusters):
            if sample in group:
                return i
    
    def backTrack(assignment, chr_mat_slice):
        '''
        int, int, array -> int
        finds the place to go back to after the 'race'.
        basically finds the last 1 in the assignment column
        '''
        try:
            return 1 + np.argwhere(chr_mat_slice[:, assignment]==1)[-1][0]
        except IndexError:
            print('backtrack error at', assignment, chr_mat_slice)
            return 'ERROR'


    #TODO
    global matrices
    global breaks 
    global overall_cluster
    global group_names 
    global sample_list 
    global max_score 
    global penalty 

    
    
    #make a mask to hide scores for your own group
    self_group = getSelfGroup(sample)
    self_mask = np.full(len(group_names), 1)
    self_mask[self_group] = 0

    #iterate through the condensed matrices of this one sample
    data = np.array([matrix.loc[sample] for matrix in matrices])
    data[data==0] = penalty

    #transform the data to accomodate for self group
    min_score = penalty * (len(group_names) - 1)
    masked_data = data * self_mask
    self_hits =  np.sum(masked_data, axis = 1) <= min_score
    self_misses = ~self_hits
    masked_data[self_hits, self_group] = 1
    masked_data[self_misses, self_group] = penalty 
    data = masked_data

    #put a black thing as the beginning
    segments = [(10, -1)]
    c = 0 #cutoff
    i = 0 #iterator
    for b in breaks:
        chr_mat = data[i:b]
        chr_mat[chr_mat == 0] = penalty
        score = np.zeros(len(group_names))
        dropped = np.full_like(score, 1)
        while i < b:
            pos_data = chr_mat[i]
            # score = np.clip((score + pos_data) * dropped, 0, max_score) #updates the score, but puts the dropped ones to zero.
            pos_score = pos_data * dropped
            score = score + pos_score
            
            # #diagnostic
            print(sample, c, i, score, pos_data)

            dropped[np.where(score <= 0)] = 0.
            if np.sum(dropped) == 0:
                assignment = np.argmin(score) #find the index with the last non-zero value
                l = backTrack(assignment, chr_mat[c:i]) #length of this segment
                if l == 'ERROR':
                    print('error', b, c, i, score, pos_data, sample, chr_mat[c:i])
                    raise Exception

                c = c + l #move the cutoff up
                i = c #move the iterator up
                segments.append((l, assignment)) #segment contains only length and assignment
                score = np.zeros_like(score)
                dropped = np.full_like(score, 1)
            else:
                score = np.clip(score, 1, max_score)
                i += 1
        else:
            #chooses the highest tally as the assignment for the last one.
            # print(sample, c, i)
            assignment = np.argmax(score)
            segments.append((i - c, assignment))
            segments.append((1, -1))
            c = i

    
    #flatten
    return segments

'''matrix (condensed but not aggregate), output file -> None (output)
outputs the whole graph to tabular format
input_type should bear the format {chr: strain, [(start, end, group)]}'''
def tab_output(composition, samplelist, colortable, blocksize, outpath):
    import CytoscapeEncoder as ce
    import copy
    samplelist = sorted(samplelist)
    matrix = copy.deepcopy(composition)
    
    #expansion step
    for chr in matrix:
        for sample in samplelist:
            temp_result = []
            for section in matrix[chr][sample]:
                start = section[0]
                end = section[1]
                color = colortable[section[2]]
                for x in range(end - start + 1):
                    temp_result.append(ce.toHexColor(color))
            matrix[chr][sample] = temp_result
    
    #output stepmatrix
    with open(outpath, 'w') as output:
        output.writelines(['\t'.join(['Chromosome', 'Position'] + samplelist)])
        output.write('\n')
        for chr in matrix:
            for index, position in enumerate(zip(*[matrix[chr][sample] for sample in samplelist])):
                num = str(index * blocksize)
                output.write('\t'.join([chr, num] + list(position)))
                output.write('\n')     
    
    return