'''
Created on Sep 16, 2013

@author: javi
'''

# Breaks a given FQ sequence into 150 bp fragments, with 50bp overlab on either side.
# Input: name of fastq file
# Output: chopped up fastq file, separate copy, in a subfolder called "Chopped"
# Notes:  Only works for fastQ!

import os
from os import listdir
from os.path import isfile, join
import re
from FqBlock import FqBlock
   
 

    
def chop(newdirectory, newchopdirectory, newsourcename, newidentifier, newlength, newoverlap, newlog):
# initialize
    directory = newdirectory
    sourcename = newsourcename
    chopdirectory = newchopdirectory
    source = None
    product = None
    identifier = newidentifier
    length = newlength
    overlap = newoverlap
    nextblockid = ""
    
# this is now done in main
#     directory = raw_input("Please specify working directory, in full \n")
#     sourcename = raw_input("please specify source file name, in full \n")
#     identifier = raw_input("Please input sequence identifier, not including the @ symbol \n")
#     length = input("Please specify the length of moving window \n")
#     overlap = input("Please specify the length of overlap \n")
#     chopdirectory = directory + "/Chopped"
    
     
    source = open(sourcename, "r")
    print "\nreading file %s ..." % sourcename
    data = source.read()
    print "done."
    
    product = open("%s/%s-chopped.fq" % (chopdirectory, sourcename.split(".")[0]), "w")
    
    # identifying blocks. Each block processes and writes itself into the file.
    eof = False
    currentblock = FqBlock(length, overlap)
    
    while(not eof):
        g = re.match("(?is)(@%s.*?)\n([a-zA-Z\n]*?)\n\+\n(.*?)(@%s.*)" % (identifier, identifier), data)
        if not g:
            g = re.match("(?is)(@%s.*?)([a-zA-Z\n]*?)\n\+\n(.*?)$" % (identifier), data)
            eof = True
        currentblock.setSeqid(g.group(1))
        currentblock.setSequence(g.group(2).replace('\n', ''))
        currentblock.setQscore(g.group(3).replace('\n', ''))
        currentblock.write(product)
        try:
            data = g.group(4)
        except:
            log.write("\n%s reached eof" % sourcename)
        
            
#         currentblock.reset()
#         currentblock.setSeqid(nextblockid)
#         inblock = True
#         
#         while(inblock):
#             temp = source.readline()
# #             print temp
#             if (re.match("@" + identifier + ".*?\n", temp) or temp == ""):
#                 inblock = False
#                 if not currentblock.seqid == "":
#                     currentblock.write(product)
#                 nextblockid = temp.rstrip('\n')
#                 
#             else:
#                 currentblock.addRaw(temp.rstrip('\n'))
#                 
#             if temp == "":
#                 eof = True   
            


if __name__ == '__main__':
    #modify these as needed
    directory = raw_input("Please specify working directory, in full \n")
    identifier = raw_input("Please input sequence identifier, not including the @ symbol \n")
    length = input("Please specify the length of moving window \n")
    overlap = input("Please specify the length of overlap \n")
    
    #Do not modify
    os.chdir(directory)  
     
    if not (os.path.isdir(directory + "/Chopped")):
        os.mkdir("Chopped")
    
    chopdirectory = directory + "/Chopped"
        
    log = open(chopdirectory + "/log.txt", "w")
    log.write("Run inputs are: \n%s\n%s\n%s\n%s\n" % (directory, identifier, length, overlap))
    
    onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
    for f in onlyfiles:
        print "\nProcessing %s ..." % f
        log.write("\nProcessing %s ... " % f)
        
        try: 
            namesplit = f.split(".")
            if namesplit[len(namesplit)-1] == "fq":
                chop(directory, chopdirectory, f, identifier, length, overlap, log)
            else:
                log.write("\nerror, %s is not a fq file" % f)
        except:
            log.write("\nerror, %s is not a fq file" % f)
 
              
    print "\nend of script"