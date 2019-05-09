import sys
import unittest
import pandas as pd
import numpy as np
from pathlib import Path

import LoadInput


class testLoadInput(unittest.TestCase):

    def setUp(self):
        pass
    
    def testLoadToPandas(self):
        path = Path('testdata.tsv')
        res, sample_list = LoadInput.loadToPandas(Path('testdata.hdf'), Path('testdata.tsv'), 'ME49')

        idx_vals = [
            ('TGME49_chrXIV', 1287), 
            ('TGME49_chrXIV', 1404),
            ('TGME49_chrXII', 1524),
            ('TGME49_chrXII', 1958),
            ('TGME49_chrIII', 3030),
            # ('TGME49_chrIII', 390077),
            ('TGME49_chrI', 1857791),
            ('TGME49_chrI', 1857922)
            ]

        
        vals = [
            ['G', 'A', 'A', 'A', 'G'],
            ['T', 'C', 'C', 'C', 'T'],
            ['T', 'C', 'C', 'C', 'T'],
            ['G', 'A', 'A', 'A', 'G'],
            ['G', 'T', 'T', 'T', 'G'],
            # ['G', 'C', 'C', 'C', 'C'],
            ['A', 'G', 'G', 'G', 'A'],
            ['C', 'C', 'G', 'G', 'C']
            ]

        idx = pd.MultiIndex.from_tuples(idx_vals, names=['CHROM', 'POS'])
        expected = pd.DataFrame(np.array(vals), index = idx, columns = ['ME49', 'S1', 'S2', 'S3', 'S4'])


        expected, res = expected.align(res)
        print(res)
        print(expected)
        assert np.all(np.equal(res.index, expected.index))

        assert np.all(np.equal(res.values, expected.values))
    
    # def testFilter(self):
    #     data = np.array([[1,2,'C', 'C', 'G', 'G', 'C'], [1,2,'G', 'A', 'A', 'A', 'G'], [1,2,'G', 'G', 'G', 'G', 'G'], [1,2,'C', 'A', 'C', 'C', 'C']])
    #     res = list(map(LoadInput.filter, data))
    #     expected = [True, True, False, False]

    #     assert res == expected



    
if __name__ == '__main__':
    unittest.main()