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
            if len(domainScores) > 0:
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
    
    count = 1
    for domain, domSeq in zip(domains, domSeqs):
        fields = re.split("[\s]+", domain)
        score = float(fields[5])
        coords =  (int(fields[9]), int(fields[10]))
        sequence = "".join(re.findall("{0}\s+\d+\s+(.*?)\s.*?".format(re.split(" ", name)[0]), domSeq))
        cys = sequence.upper().count("C")
        
        if score > 0.00001 or coords[1] - coords[0] < 90:
            continue
        
        if cys < 4:
            degenerate = True
        else:
            degenerate = False
        
        results[count] = (family, score, coords, cys, degenerate)
        count += 1
    
    return name, results
        
        
    
    