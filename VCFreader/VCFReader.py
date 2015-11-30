'''
Created on Oct 4, 2013
@author: javi
'''
import os
import sys
from os import listdir
from os.path import isfile, join
import re

'''    
Regarding the Regexs: 
1. Flags: Multiline, Dotall
#(?:\w+?\t){9}    skips all the ##INFO stuff, and finds the line with the column headers.
                  then grabs all the headers except the first 9 (which are predefined)
                  
((?:\w+?\t)+\w+?)    after skipping the first 9, this part grabs all the remaining headers (sample names) as group(1)

\n(.*)    Skips the newline character, and grabs the remainder of the document (the reads) as group(2)

2. Flags: None
Separates out each sample name (from the whole line) and puts them into an array.

3. Separates out each read and put into array

This serves to help with mapping each number to the sample. But I wonder if we even need to.. 

'''

def analyzeVCF(data):
    separatedata = re.search("(?ms)#(?:\w+?\t){9}((?:\w+?\t)+\w+?)\n(.*)", data)
    if separatedata:
        samplelist = re.findall("\w+?(?=(?:\t|$))", separatedata.group(1))
        readlist = re.split("\n", separatedata.group(2))

        
    pass 


if __name__ == '__main__':

#Load the VCF file.
    path = ""
    directory=""
    filename=""
    pathfrag = re.match("(.*)/(.*?)\.(.*)$")
    try:
        if pathfrag == True and (pathfrag.group(3) == "vcf" or pathfrag.group(3) == "bcf"):
            directory = pathfrag.group(1)
            filename = pathfrag.group(2)
            os.chdir(directory)
            os.mkdir("%s/analyzed" % directory)
#Loads the entire file to memory, then passes it to the analyzeVCF function to do the work. 
#analyzeVCF should return a single string to be written to the output file at directory/analyzed/filename.txt
            data = open("filename","r").read()
            open("analyzed/%s.txt","w").write(analyzeVCF(data))     
        else
            print "%s is not a valid file" % path
            sys.exit()
    except IndexError:
        print "%s is not a valid file path" % path
        sys.exit()
    


#Set Chromosome Track (scaffolds too I guess). You'll also need to set the map for the samples. 

#Set read head, increment by 500 each time.

#Because the VCF should be sorted by chromosome number, just read by line. Each line will either count for this group or the next one. 

#For each group, you'll need to store all the type information from the genotype column of the VCF file. This will probably be in the form of an array. 

#You'll need a map to link the information to filename, as they would not be in the same line. The order is the same though.. 


separatedata = re.search("(?m)#(?:\w+?\t){9}((?:\w+?\t)+\w+?)\n((?:.*\n)+)", data)
