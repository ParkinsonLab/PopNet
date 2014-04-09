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

#  The indices where each party's data appears in a 4-element list.
NDP_INDEX = 0
GREEN_INDEX = 1
LIBERAL_INDEX = 2
CPC_INDEX = 3

# A list of the indices where each party's data appears in a 4-element list.
PARTY_INDICES = [NDP_INDEX, GREEN_INDEX, LIBERAL_INDEX, CPC_INDEX]
PARTY_STRINGS = ["NDP", "GREEN", "LIBERAL", "CPC"]

# for simulation data (part 2):
POPULATIONS = [5, 50, 500, 5000]
POP_DESCRIPTORS = ['tiny', 'small', 'medium', 'large']

def str_to_party_num(party_string):
    '''(str) -> int
    Return the numerical index associated with a given party.'''

    if party_string == "NDP":
        return NDP_INDEX
    elif party_string == "GREEN":
        return GREEN_INDEX
    elif party_string == "LIBERAL":
        return LIBERAL_INDEX
    elif party_string == "CPC":
        return CPC_INDEX


def num_to_party_str(party_int):
    '''(int) -> str
    Return the string associated with a given party index.'''

    return PARTY_STRINGS[party_int]


def parties_to_strings(parties_by_ints):
    '''(list of int) -> list of str
    Given a list of parties by their ints, returns the list of associated strings.
    For example: [2, 2, 1] -> ['LIBERAL', 'LIBERAL', 'GREEN']'''
    
    strs = []
    for p in parties_by_ints:
        strs.append( PARTY_STRINGS[p] )
    return strs


def break_ties(winners):
    '''(list of integers) -> int

    Randomly chooses a winner from the list of possible winners.'''
    return random.choice(winners)


def make_matrix(n, m):
    '''(int, int) -> list of NoneType

    Creates and returns an n by m list, filled with None.'''

    matrix = [None] * n
    for i in range(len(matrix)):
        matrix[i] = [None] * m

    return matrix


def read_range_borda(filename, num_ridings, voters_per_riding): #works
    votes = []
    votes_in_ridings = make_matrix(num_ridings, voters_per_riding)
    with open(filename, 'r') as f: #works to create a single list of ints
        first_row = next(f) 
        r = csv.reader(f)
        for row in r: 
            for index in range(len(row)): 
                if index>1:
                    votes += [int(row[index])]
        
        x = np.array(votes)
        new_x = np.reshape(x, ((len(votes)/4), 4))
        borda_ballots = new_x.tolist()
        
        for i in range(num_ridings):
            for j in range(voters_per_riding):
                votes_in_ridings[i][j] = borda_ballots[i*voters_per_riding + j]
    print "done"    
    return votes_in_ridings

def diamond(max_itr):
    return rec_diamond(max_itr - 1, 0)

def rec_diamond(max_itr, cur_itr):
    if not cur_itr < max_itr:
        return "*" * (2 * cur_itr + 1) + "\n"
    
    spaces = " " * (max_itr - cur_itr)
    stars =  "*" * (2 * cur_itr + 1)
    line = spaces + stars + spaces + "\n"
    return line + rec_diamond(max_itr, cur_itr + 1) + line

def one_itr(matrix):
    result = np.empty_like(matrix)
    length = len(matrix)
    for x in xrange(length):
        for y in xrange(length):
            ncount = 0
            for a in xrange(x-1, x+2):
                for b in xrange(y-1, y+2):
                    if not out_of_bound(a, b, matrix):
                        ncount += matrix[a, b]
            result[x, y] = transform(matrix[x,y], ncount)
    return result
            
    
def out_of_bound(x, y, matrix):
    return not (0 <= x < len(matrix) and 0 <= y < len(matrix[0]))

def transform(value, neighbors):
    if value == 0:
        if neighbors == 3:
            return 1
    else:
        if 3 <= neighbors <= 4:
            return 1
    return 0
        
def game_of_life(matrix, num):
    print matrix
    for x in xrange(num):
        matrix = one_itr(matrix)
        print matrix
        
        
if __name__ == '__main__':
    
#     a = numpy.ones(36)
#     a[::2] = 0
#     print a.reshape(6, 6)
    input = np.zeros((5, 5))
    input[1:4, 2] = 1
    game_of_life(input, 5)
    
        
    
    
            


    





