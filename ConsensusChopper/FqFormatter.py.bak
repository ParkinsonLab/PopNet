'''
Created on Sep 23, 2013

@author: javi
'''
import os
from os import listdir
from os.path import isfile, join
import re
from FqBlock import FqBlock

def format(directory, f):
    source = open("f", "r").read()
    temp = re.match("{(@.*?\n\+\n")

if __name__ == '__main__':
     #modify these as needed
    directory = raw_input("Please specify working directory, in full \n")
    
    #Do not modify
    os.chdir(directory)  
     
    if not (os.path.isdir(directory + "/Formatted")):
        os.mkdir("Formatted")
    
    chopdirectory = directory + "/Formatted"
        
    log = open(chopdirectory + "/log.txt", "w")
    log.write("Run inputs are: \n%s" % (directory))
    
    onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
    for f in onlyfiles:
        print "\nProcessing %s ..." % f
        log.write("\nProcessing %s ... " % f)
        
        try: 
            if f.split(".")[1] == "fq":
                format(directory, f)
            else:
                log.write("\nerror, %s is not a fq file" % f)
        except:
            log.write("\nerror, %s is not a fq file" % f)
 
              
    print "\nend of script"