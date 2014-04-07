'''
Created on Sep 17, 2013

@author: javi
'''
#Temporarily holds a fastq block prior to writing into the file.
#Starts out blank
#Sequence and qscores are kept separately, and are both extensible. You can't really delete though.

import re

class FqBlock:
    '''
    classdocs
    '''


    def __init__(self, blocklength, overlap):
        self.sequence = ""
        self.qscore = ""
        self.seqid = ""
        self.blocklength = blocklength
        self.overlap = overlap
        
    def build(self):
    #This is the main part of the class. Takes the original block and 
    #outputs the block divided into X basepair blocks with Y overlap on either side.    
          
        tempseqlist = []
        tempqscorelist = []
        
        while(len(self.sequence) > self.blocklength):
            tempseqlist.append(self.sequence[:self.blocklength])
            self.sequence = self.sequence[self.blocklength-self.overlap:]
            tempqscorelist.append(self.qscore[:self.blocklength])
            self.qscore = self.qscore[self.blocklength-self.overlap:]
        
        tempseqlist.append(self.sequence)
        tempqscorelist.append(self.qscore)
        
        output = ""
        for x in range (0, len(tempseqlist)):
            output = output + "%s %d\n%s\n+\n%s\n" % (self.seqid, x,  tempseqlist[x], tempqscorelist[x])
        
#         print output
        return output
    
    def write(self, output):
        output.write(self.build())
        
    def setSequence(self, newsequence):
        self.sequence = newsequence  
    
    def setQscore(self, newqscore):
        self.qscore = newqscore
            
    def addRaw(self, newraw):
        self.raw = self.raw + newraw
        
    def setSeqid(self, newseqid):
        self.seqid = newseqid
        
    def __print(self):
        print self.build()
        
    def reset(self):
        self.__init__(self.blocklength, self.overlap)
        