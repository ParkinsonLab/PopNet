'''
Created on Oct 11, 2013

@author: javi
'''

#system imports
import numpy as np
import pandas
import os

from pathlib import Path
from multiprocessing import Pool
from sklearn.metrics import pairwise_distances
from sklearn.metrics.pairwise import cosine_similarity
from sklearn.preprocessing import normalize

#module imports
from IOTools import mcl, buildMatrix

       
def primaryCluster(df, sample_list, cluster_params, logger):
    '''
    main function for primary clustering
    we're returning a flat list of blocks with a second list of chr break points
    '''
    def groupBySection(positions, section_length):
        '''
        ndarray, float -> [ndarray]
        return a nested list containing positions within sections in sequence
        '''
        def findBI(positions, boundary):
            #ndarray, float -> int
            #so we're finding the index of the position just before the boundary.
            #this should work since we're only cutting off the back
            #edge case of where the last position happens to be valid
            dists = positions - boundary
            try:
                idx = np.argmax(dists[dists <= 0]) + 1
            except ValueError:
                print('Warning: empty block')
                idx = 0

            #if we want to return the last position, return -1 isntead
            if idx == len(positions) - 1:
                return -1
            else:
                return idx
        
        #find break points
        start = 0
        end = positions[-1]
        boundaries = [section_length * (x+1) for x in range(int(np.ceil((end - start) / section_length)))]
        boundary_idc = [findBI(positions, boundary) for boundary in boundaries]
        boundary_idc.insert(0, 0) #put a zero in the front
        res = [positions[boundary_idc[i - 1]: boundary_idc[i]] for i in range(1, len(boundary_idc))]

        return res

    tab_path = Path('./tab.txt')
    recordTab(sample_list, tab_path)

    section_length = cluster_params.getSectionLength()
    ival = cluster_params.getIVal()
    pival = cluster_params.getPiVal()

    res = []
    chr_breaks = []
    chr_names = []
    prev = 0
    # pool = Pool(processes = 6, initializer=mclInit, initargs=(sample_list, tab_path, ival, pival))


    #slicing and handing out jobs
    global df_chr_global
    for chr_name, df_chr_global in df.groupby(level=0, sort=False):
        with Pool(processes = 4, initializer=mclInit, initargs=(sample_list, tab_path, ival, pival)) as pool:
            print(chr_name)
            chr_names.append(chr_name)

            positions = np.array(df_chr_global.index.get_level_values('POS'))
            sections = groupBySection(positions, section_length)

            clusters = pool.starmap(mclWorker, [(i, (slice(None),section)) for i, section in enumerate(sections)], chunksize = 5)
            res += clusters

            chr_breaks.append(prev + len(clusters))
            prev += len(clusters)
    
    #post-processing
    analysis = analyzeClusters(res)
    logger.info(analysis)

    return res, chr_names, chr_breaks

def mclInit(sample_list_arg, tab_path_arg, ival_arg, pival_arg):

    global sample_list
    global tab_path
    global ival
    global pival

    sample_list = sample_list_arg
    tab_path = tab_path_arg
    ival = ival_arg 
    pival = pival_arg

    return


def mclWorker(i, section_pos_list):
    '''
    main function of an mcl worker
    gets an mcl result from a block
    '''

    def matrixFromBlock(section):
        '''
        takes one block (N positions) and returns a similarity matrix
        uses pairwise_distances, metric = cityblock

        need to make sure samples are rows though
        '''
        onehot = pandas.get_dummies(section)
        dists = cosine_similarity(onehot.values)
        #possible reformat the table here
        #rescale values to 0, 1
        dists = dists - np.min(dists)
        dists = dists / np.max(dists)

        return dists

    

    #use initialized variables
    global sample_list
    global tab_path
    global ival
    global pival

    # print('thread {0} is doing section {1}'.format(os.getpid(), i))
    # get dists
    try:

        section = df_chr_global.loc[section_pos_list, :]

        if section.empty:
            dists = np.full((len(sample_list), len(sample_list)), 1)
        else:  
            dists = matrixFromBlock(section.T)

        #reformat into mcl.   
        mcl_matrix = buildMatrix(dists, sample_list)
        
        # call mcl
        groups = mcl(mcl_matrix, tab_path, ival, pival, raw=False)

        return groups
    except:
        print('mcl worker crashed in {0}'.format(os.getpid()))


    
def recordTab(sample_list, tabpath):
    #writes the single tab file
    with open(tabpath, 'w') as persistentTab:
        for index, key in enumerate(sample_list):
            persistentTab.write("{0} {1}\n".format(str(index), key))
                    


def analyzeClusters(clusters):
    '''
    given a nested list of clusters, return some summary stats
    currently we do average clusters and unclustered blocks, unclustered where average size = 1
    '''

    def clusterStats(cluster):
        '''returns the average number of clusters formed'''
        #ag.decomposed assumed to give you a nested list back

        # return np.mean([len(x) for x in cluster])
        return len(cluster)

    avg_sizes = np.array(list(map(clusterStats, clusters)))
    unclustered = np.sum(avg_sizes[avg_sizes == 1])
    total = len(clusters)

    template = "\
Primary clustering results:\n\
Total: {0}\n\
Unclustered: {1}\n\
Average size: {2}\n"

    return template.format(total, unclustered, np.mean(avg_sizes))


if __name__ == '__main__':

    pass
