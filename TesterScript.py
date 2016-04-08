'''
Created on Sep 18, 2013

@author: javi
'''
import os
from os import listdir
from os.path import isfile, join
import re
import numpy as np
from subprocess import call
import time
import random
import csv
import SimilarityContigs as simcon
import Proteins.HMMParser as hmp
import Proteins.DomainFamily as dmf
import Proteins.SequenceSelect as ss
import Proteins.SRSCytoscape as srsce
import multiprocessing as mp


def list_maker():
    return [1,2,3,4,5], [6,7,8,9,0]

def c(w, x, y, z):
    print((w, x, y, z))
    return x
    
def b(var_tuple):
    
    var_a = var_tuple[0]
    var_b = var_tuple[1]
    var_c = var_tuple[2]

#     var_a = 1
#     var_b = 2
#     var_c = 3

    var_d = 99
    
    return c(var_a, var_b, var_c, var_d)

def a(i):
    
    def list_maker2(var_c, list_1, list_2):
#         return [(var_c, x, list_2) for x in list_1]
        return [(var_c, i, list_2) for i in list_1]
    
    n = 3
    m = 4
    o = 5
    list_a, list_b = list_maker()

    pool = mp.Pool()
    param_list = list_maker2(n, list_a, list_b)
    
    results = pool.map(b, param_list)
    print(results)

if __name__ == '__main__':
    
    directory = '/data/new/javi/'
    os.chdir(directory)
    tempname = 'zzz.txt'
    with open(tempname, 'w') as input_type:
        input_type.write('hi')
    os.remove(tempname)
    





