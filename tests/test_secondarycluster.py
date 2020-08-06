import unittest
import pandas as pd
import SecondaryCluster as sc
import ParamWrapper as pw 
import numpy as np

class SecondaryClysterTestCase(unittest.TestCase):

    def setUp(self):
        pass

    def testGroup(self):
        overall_matrix = pd.DataFrame([[15,15,2,2],[15,15,2,2],[2,2,15,15],[2,2,15,15]], columns = ['W', 'X', 'Y', 'Z'], index = ['W', 'X', 'Y', 'Z'])
        tab_path = 'tab2.txt'
        group_path = 'groups.txt'
        s2_params = pw.ParamWrapper()
        s2_params.setIMax(10)
        s2_params.setIMin(2)
        s2_params.setIStep(0.5)
        s2_params.setPiMax(10)
        s2_params.setPiMin(1)
        s2_params.setPiStep(0.5)
        group_names, overall_clusters = sc.group(overall_matrix, tab_path, group_path, s2_params, autogroup=True)

        expected_gn = ['A', 'B']
        expected_oc = [['W', 'X'], ['Y', 'Z']]
        
        assert expected_gn == group_names
        assert np.all(np.equal(expected_oc, overall_clusters))
    
if __name__ == '__main__':
    unittest.main()