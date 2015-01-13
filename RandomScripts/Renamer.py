'''
Created on Oct 28, 2014

@author: javi
'''

if __name__ == '__main__':
    
    path = "/data/new/javi/yeast/SECONDNewGenomes/cere/strains"
    oldname = "genome.fa"
    
    import os
    os.chdir(path)
    for dir in os.listdir(path):
        print("processing {}".format(dir))
        oldfile = "{0}/assembly/{1}".format(dir, oldname)
        newfile = "{0}/assembly/{0}.fasta".format(dir)
        if os.path.isfile(oldfile):   
            os.rename(oldfile, newfile)
        elif not os.path.isfile(newfile):
            print("bad folder?? @ {}".format(dir))
        
    print("end of script")