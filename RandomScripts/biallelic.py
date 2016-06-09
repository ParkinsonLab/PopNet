'''
Created on May 5, 2016

@author: javi
'''

if __name__ == '__main__':
    import numpy as np
    
    file_path = '/data/new/javi/plasmo/pipeline/matrix/results2.txt'
    
    data = np.genfromtxt(file_path, skip_header=1, delimiter='\t', dtype='string')[:,2:]
    count = 0
    
    for i, row in enumerate(data):
        if i%10 == 0:
            print('processed {} lines.'.format(i))
        tmp = set(row)
        if len(tmp) > 2:
            print(tmp)
            count += 1
    
    print('{0} out of {1} positions have more than two alleles.'.format(count, len(data)))