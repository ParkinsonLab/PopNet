'''
Created on Feb 22, 2015

@author: javi
'''

from matplotlib import pyplot as plt
from matplotlib import mlab
import numpy as np

if __name__ == '__main__':
    
    a = [[1,1,1,1,2,6],
         [7,9,6,7,5,2],
         [3,2,1,3,2,1],
         [3,8,7,4,2,1],
         [3,7,5,6,2,10],
         [8,3,1,2,2,1],
         [1,2,1,2,3,1]]
    arr = np.array(a)
    
    b = mlab.PCA(arr)
    
    plt.plot(b.Y[:,0], b.Y[:,1], 'o', color = 'red')
    print(b.Wt)
    print('fracs\n' + str(b.fracs))
    plt.show()
    
