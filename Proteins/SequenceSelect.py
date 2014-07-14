'''
Created on Jun 18, 2014

@author: javi
'''
'''picks out the sequences from the fasta file that were hits in the HMM search. 

Overview:

Requires: filename, identifying the strain. Folder of all the fastas. folder with all hmmr files. 

1. read the hmmr files, get name + id of all proteins
2. for each protein id, look in the corresponding fasta file and grab that chunk. 
3. append the stain name to the front of the protein ID so it'd be like ARI_TGME49_27900
4. write one file for each strain.
'''
import re
import os 
'''String -> (name, list)
reads the HMMR file'''
def parseHMMR(filepath):
    with open(filepath, 'r') as input:
        data = input.read()
        proteins = re.findall("(?m)^>>((.+?)\s+.+)$", data)
        print("{0} proteins read from {1}".format(len(proteins), re.match("^(.+?)[.].+$", filepath).group(1)))
    return proteins

'''String, list -> list
the returned list already have the strain name attached to it'''        
def findSequences(fastaFilepath, seqList, strain):
    with open(fastaFilepath, 'r') as input:
        data = input.read()
        results = []
        for protein in seqList:
            ID = ">{0}".format(protein[0])
            seq = re.search("(?s)>{0}\s.+?\n(.+?)(?=\n>|$)".format(protein[1]), data).group(1)
            results.append("\n".join([ID, seq]))
        return results
    
'''string, list -> none (output)
writes the fasta files'''
def writeFasta(filepath, list):
    with open(filepath, 'w') as output:
        output.write("\n".join(list))
        print("{0} proteins wrote to {1}".format(len(list), re.match("^(.+?)[.].+$", filepath).group(1)))

'''string, string, string, string -> none (output)
primary running method for this script'''
def select(hmmrpath, fastapath, outputpath, strain):
    seqList = parseHMMR(hmmrpath)
    sequences = findSequences(fastapath, seqList, strain)
    writeFasta(outputpath, sequences)



def filter(filepath, outpath):
    with open(filepath) as input, open(outpath, 'w') as output:
        data = input.read()
        lines = re.split("\n", data)
        results = []
        for index, line in enumerate(lines):
            good = ["Group_{0}:".format(index+1)]
            try:
                removeHeader = re.match("^.+?: (.+)$", line).group(1)
            except:
                continue
            elements = re.split("\s", removeHeader)
            for e in elements:
                if not re.match("^ham", e):
                    good.append(e)
            results.append(good)
        
        print("filtered {0} lines".format(index))
        
        count = 0
        for line in results:
            output.write(" ".join(line) + "\n")
            
def toName(filepath, hmmrDirectory, outpath):
    hmmrFiles = [ x for x in os.listdir(hmmrDirectory) if os.path.isfile("/".join([hmmrDirectory,x]))]
    proteins = []
    proteinDict = {}
    for f in hmmrFiles:
        proteins += parseHMMR("/".join([hmmrDirectory,f]))
    for protein in proteins:
        proteinDict[protein[1]] = re.search("\s+?(\w.+)$", protein[0]).group(1)
    
    with open(filepath) as input, open(outpath, 'w') as output:
        data = input.read()
        lines = re.split("\n", data)
        results = []
        for index, line in enumerate(lines):
            try:
                removeHeader = re.match("^.+?: (.+)$", line).group(1)
            except:
                continue
            elements = re.split("\s", removeHeader)
            results.append([re.split("[|]", x) for x in elements])
    
        for line in results:
            towrite = []
            lineNames = []
            for element in line:
                strain = element[0]
                name = proteinDict[element[1]]
                
                srs = re.search("SRS.*$", name)
                if srs:
                    name = srs.group(0)
                
                if name in lineNames:
                    name = "-"
                else: lineNames.append(name)
                towrite.append("{1}_{0}".format(strain, name))
                
            output.write(",".join(towrite) + "\n")

        
    