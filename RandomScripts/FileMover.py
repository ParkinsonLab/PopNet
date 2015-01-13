'''
Created on Oct 28, 2014

@author: javi
'''

if __name__ == '__main__':
    path = "/data/new/javi/yeast/SECONDNewGenomes/cere/strains"
    dest = "/data/new/javi/yeast/SECONDNewGenomes/fasta"
    import os
    import shutil
    os.chdir(path)
    if not os.path.isdir(dest):
        os.mkdir(dest)
    
    for dir in os.listdir(path):
        if dir != dest and os.path.isdir(dir):
            for file in os.listdir("{0}/{1}/assembly/".format(path, dir)):
                if file.endswith("fasta"):
                    filepath = "{0}/{1}/assembly/{2}".format(path, dir, file)
                    shutil.move(filepath, dest)
                    print("{0} moved".format(file))