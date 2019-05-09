import unittest

import pandas as pd 
import numpy as np 
import ChromosomePainter as gc

class gcTestCase(unittest.TestCase):

    def setUp(self):
        pass
    
    def testCondense(self):
        sample_list = ['A','B','C']
        overall_clusters = [['A', 'C'],['B']]
        matrices = [pd.DataFrame([[1,0,1],[0,1,0],[1,0,1]], columns = ['A', 'B', 'C'], index = ['A', 'B', 'C'])]*2
        group_names = ['X', 'Y']
        condensed = gc.condenseToGroupMatrix(matrices, group_names, overall_clusters, sample_list)

        # print(condensed)

    def testChromPainting(self):
        sample_list = ['A']
        overall_clusters = [['A', 'C'],['B'], ['D', 'E']]
        data = [[0,0,1],
                [0,0,1],
                [0,0,1],
                [0,0,1],
                [0,0,0],
                [0,0,0],
                [0,0,0],
                [0,1,0],
                [0,1,0],
                [0,1,0],
                [0,0,0],
                [0,1,0],
                [0,1,0],
                [0,0,0],
                [0,0,1],
                [0,1,1],
                [0,0,1],
                [1,0,1],
                [1,0,1],
                [1,0,1],
                [1,0,1],
                [1,0,1],
                [1,0,1],
                [0,1,1],
                [0,1,1],
                [0,1,1],
                [0,1,1],
                [0,1,1],
                [0,0,0],
                [0,0,0],
                [0,0,0],
                [0,0,0]]

        condensed_matrices = [pd.DataFrame([row], columns = ['X', 'Y', 'Z'], index = ['A']) for row in data]
        chr_breaks = [len(data)]
        group_names = ['X', 'Y', 'Z']
        # condensed_matrices = gc.condenseToGroupMatrix(matrices, group_names, overall_clusters, sample_list)
        paintings = gc.getChromosomePaintings(condensed_matrices, chr_breaks, overall_clusters, group_names, sample_list)
        expected = [[(10, -1), (4, 2), (3, 0), (3, 1), (1, 0), (2, 1), (1, 0), (14, 2), (4, 0), (1, -1)]]
        assert paintings == expected

if __name__ == '__main__':
    unittest.main()