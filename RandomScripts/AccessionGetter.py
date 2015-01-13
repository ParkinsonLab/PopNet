'''
Created on Nov 12, 2014

@author: javi

Gets the run accession from secondary sample accession. Or modify for other uses. 
'''

if __name__ == '__main__':
    prefix = 'ERS'
    pattern = '(ERS[0-9]+)\D'
    masterpattern = '(?m)^(\w+)\s+(\w+)$'
    import re
    
    file = '/home/javi/RawAccessions.txt'
    outfile = '/home/javi/Accessions.txt'
    masterfile = '/home/javi/MasterAccList.txt'
    with open(file, 'r') as input:
        data = input.read()
    
    with open(masterfile, 'r') as master:
        masterdata = master.read()
        
    hits = re.findall(pattern, data)
    masterdict = {}
    for item in re.findall(masterpattern, masterdata):
        masterdict[item[0]] = item[1] 
    
    with open(outfile, 'w') as output:
        for i in hits:
            if i in masterdict:
                output.write('{0}\t{1}\n'.format(i, masterdict[i]))
            else:
                print(i + " is not in the master file!")