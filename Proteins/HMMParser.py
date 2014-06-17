'''
Created on Jun 13, 2014

@author: javi
'''
'''Container script for HMM parsing routines'''
import re
import numpy

'''filepath -> dictionary of nested dictionary
reads a hmm search file and get the protein domain scores.

return structure

{proteinID : {domainID: (family, score)}}
'''
def read(filepath, family):
    results = {}
    pattern = re.compile('>> (.+?)\n\n(?=>>|\n)', re.DOTALL)
    
    with open(filepath) as f:
        fStr = f.read()
        matches = re.findall(pattern, fStr)
        for block in matches:
            pid, domainScores = parseBlock(block, family)
            results[pid] = domainScores
    
    print(str(len(matches)) + " proteins from " + filepath + " from " + family + " read.")
    return results
            
'''string -> nested dict
given a block representing one protein, parse info into {name: {domainnum :  (family, score, length, cys, degenerate)}}'''
def parseBlock(block, family):
    results = {}
    name = re.match("(.+)\s", block).group(1)
    
    
    domains = re.search("(?s)^.*?----\n   (.+)(?=\n\n  Alignments)", block).group(1).split("\n   ")
    domSeqs = re.split("==", block)[1:]
    
    for domain, domSeq in zip(domains, domSeqs):
        fields = re.split("[\s]+", domain)
        score = float(fields[4])
        length =  int(fields[10]) - int(fields[9])
        sequence = "".join(re.findall("{0}\s+\d+\s+(.*?)\s.*?".format(re.split(" ", name)[0]), domSeq))
        cys = sequence.upper().count("C")
        
        if score > 0.00001 or length < 90 or cys < 4:
            degenerate = True
        else:
            degenerate = False
        
        results[int(fields[0])] = (family, score, length, cys, degenerate)
    
    return name, results
        
        
    
    