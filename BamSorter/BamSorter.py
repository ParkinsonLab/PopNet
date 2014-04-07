'''
Created on Oct 3, 2013

@author: javi
'''

import os
from os import listdir
from os.path import isfile, join
import re




if __name__ == '__main__':
    
#init variables
    folder=""

#load all bam files. probably use bamtools to parse them into txt files. 

#set one of them as the "lead". The first one will do. 

#define a regex to represent a "block". in bam files this is just one tab delimited line. 

#read the reference location and start site in the lead block, then search for that block in the other files.

#grab the results and write them into one file.  

#do this until the end of file is reached. 


#new function? call indels for each group of blocks. 

