'''
Created on Nov 12, 2014

@author: javi
'''

if __name__ == '__main__':
    logfile = '/data/new/javi/plasmo/ena/log.txt'
    accessionFile = '/home/javi/Accessions.txt'
    accessionPattern = '(?m)^(\w+)\s+(\w+)$'
    pathfile = '/data/new/javi/plasmo/ena/paths.txt'
    filePattern = 'ftp://ftp.sra.ebi.ac.uk/vol1/fastq/{0}/{1}/{1}_{2}.fastq.gz\n'
    dest = '/data/new/javi/plasmo/ena'
    
    with open(accessionFile, 'r') as acc:
        import re 
        accessions = re.findall(accessionPattern, acc.read())
    
    with open(pathfile, 'w') as paths:
        for item in accessions:
            for x in range(1,3):
                paths.write(filePattern.format(item[1][:6], item[1], x))

    
    from subprocess import call
    from os import chdir
    chdir(dest)
    call(['wget', '-i', pathfile, "-p", logfile, "-nd"])