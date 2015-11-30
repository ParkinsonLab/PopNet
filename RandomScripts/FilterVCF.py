'''
Created on Mar 10, 2015

@author: javi
'''
import os
import re

def filter(filename, kw, val, outputname):
    count=0
    with open(filename) as input, open(outputname, 'w') as output:
        for line in input:
            if line.startswith('#') or int(re.search('{}=([0-9]+?);'.format(kw), line).group(1)) > val:
                output.write(line)
            else:
                count+=1
    print('filtered {:d} SNPs'.format(count))
if __name__ == '__main__':
    
    '''
    Invoke: python FilterVCF.py [WorkingDirectory] [KW] [INT] [Filename] [Outputname]
    
    filter by a particular keyword in the vcf INFO (e.g. DP) followed by the minimum value
    '''
    import sys
    args = sys.argv
    print(args)
    workingDir = args[1]
    filename = args[4]
    outputname = args[5]
    #filter variables
    kw = args[2]
    val = int(args[3])
    
    os.chdir(workingDir)
    filter(filename, kw, val, outputname)
    
    
        