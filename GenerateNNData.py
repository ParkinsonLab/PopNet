'''
This module outputs the primary cluster results 'as-is', to be used by the neural net. We can
do the rest of popnet as well, I guess
'''
import pandas as pd
import numpy as np


def genNNData(clusters, chr_names, chr_breaks, s1_params, sample_list, out_path):

    res = np.zeros((len(clusters), len(sample_list)))

    for i, cluster in enumerate(clusters):
        for j, line in enumerate(cluster):
            for e in line:
                res[i, sample_list.index(e)] = j
    
    #create index
    sl = s1_params.getSectionLength()
    c = 0
    idx = []
    itr = iter(clusters)
    for i, chr in enumerate(chr_names):
        b = chr_breaks[i]
        for x in range(c, b):
            idx.append((chr, sl * x))
        c = b

    df = pd.DataFrame(res, columns = sample_list, index = idx)
    print(df)
    print(idx)

    df.T.to_csv(out_path, sep='\t', header=True)

if __name__ == "__main__":
    #TESTING ONLY
    import sys

    clusters = [
        [['A', 'B', 'C'], ['D', 'E']],
        [['A', 'B'], ['C', 'D'], ['E']],
        [['B', 'C', 'E'], ['A', 'D']],
        [['C', 'D', 'E'], ['A', 'B']]
    ]

    sample_list = ['A', 'B', 'C', 'D', 'E']

    out_path = '/d/data/plasmo/training_data/testdata.tsv'

    genNNData(clusters, sample_list, out_path)
