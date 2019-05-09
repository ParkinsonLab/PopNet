import unittest
import pandas as pd

import Exporter as exporter 
import AnalysisTools as at 

class _TestCase(unittest.TestCase):

    def setUp(self):
        pass
    
    def testExport(self):
        overall_matrix = pd.DataFrame([[15,15,2,2],[15,15,2,2],[2,2,15,15],[2,2,15,15]], columns = ['W', 'X', 'Y', 'Z'], index = ['W', 'X', 'Y', 'Z'])
        overall_clusters = [['W', 'X'], ['Y', 'Z']]
        group_names = ['A', 'B']
        sample_list = ['W', 'X', 'Y', 'Z']
        color_table = at.createColorTable(group_names, overall_clusters, sample_list)

        composition = [[(10, -1), (2, 0), (1, -1)], [(10, -1), (2, 0), (1, -1)], [(10, -1), (2, 1), (1, -1)], [(10, -1), (2, 1), (1, -1)]]
        prefix = 'genome'
        output = exporter.parse(overall_matrix, color_table, composition, group_names, overall_clusters, sample_list, prefix)

if __name__ == '__main__':
    unittest.main()