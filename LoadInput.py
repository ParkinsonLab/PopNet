'''
Created on Feb 16, 2014

@author: javi
'''
import re
import copy
import os
from collections import Counter as counter
import pandas
import numpy as np
import multiprocessing as mp
import logging
import sklearn
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from ChrNameSorter import sortNames


def loadToPandas(hdf_path, tsv_path, reference, filtering=True):
    
    logger = logging.getLogger()
    
    #load df. tries to get hdf first. Otherwise go to tsv and write a feather file.
    try:
        df = pandas.read_hdf(hdf_path, key='genome')
    except:
        df = pandas.read_csv(tsv_path, sep='\t')
        #position as int
        df['POS'] = df['POS'].apply(int)
      
        if (reference != None and df.columns[2] == 'REF'):
            renamed_columns={'#CHROM': 'CHROM', 
                            'REF': reference}
        else:
            renamed_columns = {'#CHROM': 'CHROM'}
        df = df.rename(index=str, columns=renamed_columns)

        #filter drift
        if filtering:
            logger.info('filtering drift')
            msk = df.apply(filter, axis = 1)
            df_filtered = df.loc[msk]
            m = len(df) - sum(msk)
            n = len(df)
            logger.info('{0} positions filtered out of {1} for drift ({2}%)'.format(m, n, round(m / n, 2)))
        else:
            df_filtered = df

        #multiindex
        chrs = sortNames(list(set(df_filtered['CHROM'])))
        cat = pandas.Categorical(df_filtered['CHROM'], chrs)
        df_filtered['CHROM'] = cat
        df_filtered = df_filtered.set_index(['CHROM', 'POS'])


        df_filtered.sort_index(0, ['CHROM', 'POS'], ascending=True, inplace=True)
        
        
        #writeout
        df = df_filtered
        df.to_hdf(hdf_path, key='genome')#move to the end
    
    return df, sorted(list(df.columns))

def filter(x):
    '''
    x is the row. Return 1 if you want to keep.
    '''
    #filter for drift, i.e. only 1
    #Is the oTHER WAY: FALSE is DRIFT to generate a mask for loc
    def filterDrift(x):
        threshold = 2
        s = set(x)
        if len(s) < 2:
            return False
        else:
            c = 0
            for e in s:
                if sum(x == e) >= threshold:
                    c += 1
                if c > 1:
                    return True
            return False

    x = x[2:]
    res = filterDrift(x) 

    return res     



    
if __name__ == '__main__':    
 pass
            
