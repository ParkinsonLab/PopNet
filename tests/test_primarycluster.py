import unittest
import PrimaryCluster as pc
import LoadInput
import IOTools
import logging
import ParamWrapper as pw
from multiprocessing import Pool
from pathlib import Path 

class testPrimaryCluster(unittest.TestCase):

    def setUp(self):
        pass

    def testPrimaryCluster(self):
        df, sample_list = LoadInput.loadToPandas('testdata2.hdf', 'testdata2.tsv', 'A', filter=False)
        logger = logging.getLogger()
        cluster_params = pw.ParamWrapper()
        cluster_params.setIVal(4)
        cluster_params.setPiVal(1.5)
        cluster_params.setSectionLength(5)
        res, chr_breaks = pc.primaryCluster(df, sample_list, cluster_params, logger)

        assert res == [[['A', 'C'], ['B']], [['A', 'C'], ['B']]]
        assert chr_breaks == [2]

    def testMCL(self):
        tab_path = Path('./tab.txt')
        ival = 30
        pival = 18

        df_chr, sample_list = LoadInput.loadToPandas('testdata2.hdf', 'testdata2.tsv', 'A', filter=False)
        sections = [[1,2,3,4,5,6,7,8,9]]
        pool = Pool(initializer=pc.mclInit, initargs=(sample_list, tab_path, ival, pival))
         
        IOTools.writeTab(sample_list, 'tab.txt')
        clusters = pool.map(pc.mclWorker, [df_chr.loc[(slice(None),section),:] for section in sections], chunksize=10)
        expected = [[['A','B','C']]]
        
        assert clusters == expected

if __name__ == '__main__':
    unittest.main()