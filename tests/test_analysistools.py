import unittest
import AnalysisTools as at 
import IOTools as io
import pandas as pd
import numpy as np

class ATTestCase(unittest.TestCase):

    def setUp(self):
        pass
    
    def testClusterToMatrix(self):
        sample_list = ['A', 'B', 'C']
        clusters = [[['A', 'C'], ['B']], [['A', 'C'], ['B']]]
        res = at.clustersToMatrix(clusters, sample_list)
        expected = [pd.DataFrame([[1,0,1],[0,1,0],[1,0,1]], columns = ['A', 'B', 'C'], index = ['A', 'B', 'C'])]*2

        for x, y in zip(res, expected):
            assert np.all(np.equal(x.index, y.index))
            assert np.all(np.equal(x.values, y.values))
    
    def testOverallMatrix(self):
        sample_list = ['A', 'B', 'C']
        clusters = [[['A', 'C'], ['B']], [['A', 'C'], ['B']]]
        res = at.clustersToMatrix(clusters, sample_list)
        overall = at.overallMatrix(res)
        expected = pd.DataFrame([[2,0,2],[0,2,0],[2,0,2]], columns = ['A', 'B', 'C'], index = ['A', 'B', 'C'])
        assert np.all(np.equal(overall.index, expected.index))
        assert np.all(np.equal(overall.values, expected.values))

    def testColorTable(self):
        group_names = ['A', 'B']
        overall_clusters = [['W', 'X'], ['Y', 'Z']]
        sample_list = ['W', 'X', 'Y', 'Z']
        color_table = at.createColorTable(group_names, overall_clusters, sample_list)

        print(color_table)
    
    def testOutPut(self):
        overall = pd.DataFrame([[2,0,2],[0,2,0],[2,0,2]], columns = ['A', 'B', 'C'], index = ['A', 'B', 'C'])
        io.writeOverallMatrix(overall, 'zzz.txt')
        
if __name__ == '__main__':
    unittest.main()