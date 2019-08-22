'''
Created on Apr 28, 2014

@author: javi
'''
'''in fulfillment of John's request for a chromosome coloring scheme where each section
is colored according to which group is based on.. MCL. zzz.
Use the createColorTable and calculateColor from ClusterPattern to visualize.'''

'''Strain -> Dict
number represents the group ID, as specified in this function'''
import numpy as np
import string
import re
import copy



def loadGroups(file, strain):
    with open(file, 'r') as input_type:
        data = input_type.read()
    global groups
    groups = {}
    for pos, section in enumerate(re.split('@', data)[1:]):
        lines = re.split('\n', section)
        lineSplit = re.split("\t", lines[1])
        groups[getGroupID(pos)] = lineSplit
#     groups['NONE'] = ['NONE']
    
    for name, group in groups.items():
        if strain in group or len(group) == 0:
            del groups[name]
            if 'NONE' in groups:
                del groups['NONE']
    
    return groups

def getGroups(strain = None):

    global groups
    if len(groups) > 0:
        if strain:
            return {k:v for k, v in groups.items() if strain not in v}
        else:
            return groups
    else:
        print('Groups not loaded.')

def getGroupID(number):
    alph = string.ascii_uppercase
    result = alph[number%26]
    
    multi = number // 26
    while multi > 0:
        result = alph[multi%26 - 1] + result
        multi = multi // 26
    
    return result 
        
        

#     #this function has changed. Now used solely for getting groups internally.
#     #must call the loadGroups function first. They are separate because one is also used externally.
#     groups = {}
#      
#     if len(groups) > 0:
#         return groups
#      
#     groups = {}
#     groups['S288C'] = ['BY4742', 'BY4741', 'CLIB324', 'W303', 'FL100', 'CEN.PK113-7D', 'SIGMA1278B', 'S288C']
#     groups['UC5'] = ['PW5', 'KYOKAI7', 'YJM269', 'YPS163', 'T7', 'EC9-8', 'Y10', 'UC5', 'ZTW1']
#     groups['AWR'] = ['VIN13', 'T73', 'CBS7960', 'RM11-1A', 'EC1118', 'JAY291', 'AWRI1631', 'M22', 'CLIB215', 'VL3', 'LALVINQA23', 'AWRI796', 'FOSTERSO', 'FOSTERSB', 'YJM789', 'CLIB382']
#     groups['NONE'] = ['NONE']
#     for name, group in groups.items():
#         if strain in group:
#             del groups[name]
#             if 'NONE' in groups:
#                 del groups['NONE']
#      
#     return groups
    
    '''8 genomes groups
    groups = {}
    groups['ME49'] = ['ME49']
    groups['VEG'] = ['VEG.S']
    groups['TYPEX'] = ['3045.', '3142.', 'ARI.V', 'TGSKN','P89.S' ]
    groups['GT1'] = ['GT1.S']
    groups['NONE'] = ['NONE']
    for name, group in groups.items():
        if strain in group:
            del groups[name]
            if 'NONE' in groups:
                del groups['NONE']
    
    return groups
    '''

    
    '''64 genome final groups
    groups = {}
    groups['GUYS'] = ["BRC_TgH_18003_GUY-MAT", "BRC_TgH_18001_GUY-DOS", "BRC_TgH_18021", \
                  "BRC_TgH_18009", "BRC_TgH_18002_GUY-KOE", "GUY-2003-MEL", "GUY-2004-ABE", \
                  "RUB", "VAND"]
    groups['TgH'] = ["BRC_TgH_20005", "CASTELLS", "TgH26044", "BRC_TgH_21016"]
    groups['ME49'] = ["COUG", "GUY-2004-JAG1", "TgCat_PRC2", "ARI", "RAY", "PRU", \
                   "B41", "ME49", "B73"]
    groups['VEG'] = ["TgShUS28", "TgCkGy2", "ROD", "M7741", "VEG", "SOU", "G662M"]
    groups['p89'] = ["TgCatBr64", "p89", "TgCatBr3", "TgCatBr15"]
    groups['TgCats'] = ["TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
                   "TgCatBr44", "MAS", "TgCatBr1", "TgCatBr34"]
    groups['GAL-DOM1/2'] = ["CAST", "FOU", "GAB5-2007-GAL-DOM1", \
                   "GAB3-2007-GAL-DOM2", "BOF"]
    groups['GAL-DOM10'] = ["GAB2-2007-GAL-DOM2", "GAB5-2007-GAL-DOM6", "GAB1-2007-GAL-DOM10", "GAB3-2007-GAL-DOM9"]
    groups['GT1'] = ["TgDogCo17", "GT1", "RH-JSR", "RH-88", "TgCkCr1", "TgCkBr141", "CAST"]
    groups['MISC']  = ["TgRsCr1", "TgCtCo5", "TgCkCr10", "TgCatBr72", "TgCATBr9", "TgCatBr26"]
    groups['NONE'] = ['NONE']
'''    
    
    
    '''ORIGINAL GROUPS
    groups['GUYS'] = ["BRC_TgH_18003_GUY-MAT", "BRC_TgH_18001_GUY-DOS", "BRC_TgH_18021", \
                  "BRC_TgH_18009", "BRC_TgH_18002_GUY-KOE", "GUY-2003-MEL", "GUY-2004-ABE", \
                  "RUB", "VAND"]
    groups['TgH'] = ["BRC_TgH_20005", "CASTELLS", "TgH26044", "BRC_TgH_21016"]
    groups['ME49'] = ["COUG", "GUY-2004-JAG1", "TgCat_PRC2", "ARI", "RAY", "PRU", \
                   "B41", "ME49", "B73", "SOU"]
    groups['VEG'] = ["TgShUS28", "TgCkGy2", "ROD", "M7741", "VEG"]
    groups['p89'] = ["G662M", "p89", "TgCatBr3", "TgCatBr15"]
    groups['TgCats'] = ["TgCatBr64", "TgCATBr5", "TgCatBr10", "TgCatBr18", "TgCatBr25", \
                   "TgCatBr44", "MAS", "TgCatBr1"]
    groups['GAL-DOM'] = ["TgCatBr72", "TgCATBr9", "TgCatBr26", "FOU", "GAB5-2007-GAL-DOM1", \
                   "GAB3-2007-GAL-DOM2", "GAB3-2007-GAL-DOM9", "GAB1-2007-GAL-DOM10", \
                   "GAB5-2007-GAL-DOM6", "GAB2-2007-GAL-DOM2", "BOF"]
    groups['GT1'] = ["CAST", "TgDogCo17", "GT1", "RH-JSR", "RH-88", "TgCkCr1", "TgCkBr141", \
                   "TgRsCr1", "TgCtCo5", "TgCkCr10"]
                   '''

'''(dict of groups) -> list of tuples
representing the group name of each strain for ease of iteration'''
def expandGroups(groups):
    results = []
    for group, strains in groups.items():
        for strain in strains:
            results.append((strain, group))
    return results

'''dict of groups -> dict of strains
reverses the dictionary so we can search by strain'''
def reverseGroups(groups):
    results = {}
    for group, strains in groups.items():
        for strain in strains:
            results[strain] = group
    return results

'''(string, 2-Nested list, string) -> int
returns the index, in the first dimension, of the group that this strain is most
closely related to, in this particular block'''    
def matchToGroup(strain, groups, block):
    counts = []
    for name, group in groups.items():
        counts.append((countMatches(strain, group, block), name))
    result = sorted(counts, reverse=True)
    #requires that the result return has a minimum of two members, else return the flag "SELF"
    if result[0][0] < 0.3:
        return "SELF"
    else: 
        return result[0][1]



'''(string, 2-Nested list, string) -> int
USED FOR CONTIG BUILDING ONLY!!! RETURNS EVERYTHING
returns the index, in the first dimension, of the group that this strain is most
closely related to, in this particular block'''    
def matchToGroupFull(strain, groups, block):
    counts = []
    cluster = ''
    for line in block:
        if strain in line:
            cluster = line #the strain must be in one of the lines. 
            break
    
    for name, group in groups.items():
        counts.append((countMatches(strain, group, cluster), name))
        
    result = sorted(counts, reverse=True)
    return result

'''(string, list, block) -> int
counts how many of the strains in the group clusters together with
the given strain within a line'''
def countMatches(strain, group, line):
    if strain in line:
        return set(group).intersection(line)
#             return len(set(group).intersection(line)) / float(len(group))
    return set([])

def findContigsComposition(strain, dataTree):
    
    minGroupPercent = 0.3 #min fraction of members present
    minDensity = 50 #min number of SNPs in this region for it to be relevant
    PENALTY_CONST = 1000.
#     penalty = reduce(lambda x, y: x+y, [len(e) for e in dataTree.values()]) / PENALTY_CONST #Penalty for a mismatch to the current match
#     #set to 0.1% of total.
#original is penalty = 8, max = 5 * penalty
    penalty = 8
    maxScore = 5 * penalty
#     print('Gap Penalty used is {}'.format(penalty))
    all_groups = getGroups() 
    other_groups = getGroups(strain)
    revGroups = reverseGroups(getGroups(""))
    self_group = revGroups[strain]
    
    def update(scores, filtered_matches):
        for racer in scores:
            if racer in filtered_matches:
                if scores[racer] <= maxScore:
                    scores[racer] += 1
            else:
                scores[racer] -= penalty
        
        for racer, score in list(scores.items()):
            if score <= 0:
                if max(scores.values()) > 0:
                    del scores[racer]
        
        return scores
    
    def update_tally(tally, scores):
        for racer in tally:
            if racer in scores:
                tally[racer].append(max(0,scores[racer]))
            else:
                tally[racer].append(0)
                
        return tally
        
    def backtrack(scores, tally, previous_matches, steps):
        for filtered_matches in previous_matches[-1 * steps :]:
            scores = update(scores, filtered_matches)
            tally = update_tally(tally, scores)
        #clears the previous matches list, because we are done backtracking.
        previous_matches = previous_matches[-1 * steps:]
        
        return scores, tally, previous_matches
    
    def filter_match(match):
        #matchlist = (number of matches, name of group)
        match_count = len(match[0])
        return float(match_count) / len(groups[match[1]]) >= minGroupPercent and match_count >= 2

    def find_window_boundary(tally):
        deltas = [tally[0]] + [tally[x] - tally[x-1] for x in range(1, len(tally))]
        for index, n in enumerate(deltas[::-1]):
            if n > 0:
                return index
  
    results = {}
    base_scores = {x:0 for x in all_groups}
    base_tally = {x:[] for x in all_groups}
    for name, chr in sorted(dataTree.items()):
        results[name] = chrResults = []
        scores = copy.deepcopy(base_scores)
        tally = copy.deepcopy(base_tally)
        length = 0
        previous_matches = []
        blocks = iter(chr)
        try:
            while True:
                block = next(blocks)
                length += 1           
                matches = matchToGroupFull(strain, other_groups, block)
                #Attempts a more sophisticated way to build the longest possible.
                #by simulating a race
                filtered_matches = [x[1] for x in matches if filter_match(x)]
                if len(filtered_matches) < 1:
                    filtered_matches.append(self_group) #appends self
                        
                scores = update(scores, filtered_matches) #update score. Eliminate dropped racers unless no more positives
                tally = update_tally(tally, scores)
                previous_matches.append(filtered_matches)
                if max(scores.values()) <= 0: #race finish when no more positives remain or end of chr.
                    winner = sorted(scores.items(), key=lambda x:x[1], reverse=True)[0][0] #The highest score at the end is the winner.
                    backtrack_steps = find_window_boundary(tally[winner]) #calculates steps needed to backtrack
                    for x in range(length - backtrack_steps): #write to results
                        chrResults.append(frozenset([winner]))
                    scores = copy.deepcopy(base_scores)
                    tally = copy.deepcopy(base_tally)                   
                    length = backtrack_steps
                    scores, tally, previous_matches = backtrack(scores, tally, previous_matches, backtrack_steps)
        except StopIteration:
            winner = sorted(scores.items(), key=lambda x:x[1], reverse=True)[0][0] #The highest score at the end is the winner.
            for x in range(length): #write to results
                chrResults.append(frozenset([winner]))
                    
    return results

   

'''([strains], dataTree) -> [[(begin, end, color]]
a variation on multicomposition that is divided by chromosomes. Used for generating the cytoscape view.'''
def cytoscapeComposition(strains, dataTree, groupfile, grey_dict):
    resultsList = {}
    loadGroups(groupfile, "")
    for strain in strains:
        resultsList[strain] = findContigsComposition(strain, dataTree)
#         resultsList[strain] = findComposition(strain, dataTree)
    
    for strain, item in resultsList.items():
        for key, value in item.items():
            resultsList[strain][key] = [next(iter(x)) for x in value]
        
        for pos_tup in grey_dict[strain]:
            chr = pos_tup[0]
            pos = pos_tup[1]
            resultsList[strain]['@' + chr][pos] = 'GREY' #looks like I always add a @ before chr names after reloading
            print('{0} at {1} {2} was Grey\'d'.format(strain, chr, pos))
    
    compiledResults = {}
    for chr in sorted(list(resultsList.values())[0].keys()):
        chrList = {}
        for strain in resultsList.keys():
            chrList[strain] = condense(resultsList[strain][chr])
        compiledResults[chr] = chrList
    
    return compiledResults
        
'''list -> list of tuples
similar to what we had in the original snps sorter thing. condenses the list by grouping
together identical elements and storing their indices.'''
def condense(list):
    resultList = []
    currentIndex=0
    for item in list:
        if len(resultList) > 0 and resultList[-1][2] == item:
            resultList[-1] = (resultList[-1][0], resultList[-1][1]+1, resultList[-1][2])
        else:
            resultList.append((currentIndex, currentIndex, item))
            currentIndex = resultList[-1][1] + 1
    return resultList

'''list of matrices -> matrix
condenses all the chrs into one, for the whole genome picture'''
def aggregate(matrixList, organism):
    import ChrNameSorter as cns
    strains = list(matrixList.values())[0].keys()
    results = {}
    for strain in strains:
#          results[strain] = [(-30, 0, 'SPACER')]
        results[strain] = []
         
    for chrName, chr in sorted(matrixList.items(), key=lambda x: cns.getValue(x[0])):
        for strain, data in chr.items():
            results[strain] += [(-10, 0, 'SPACER')]
            results[strain] += data
    return results


'''matrix (condensed but not aggregate), output file -> None (output)
outputs the whole graph to tabular format
input_type should bear the format {chr: strain, [(start, end, group)]}'''
def tab_output(composition, samplelist, colortable, blocksize, outpath):
    import CytoscapeEncoder as ce
    import copy
    samplelist = sorted(samplelist)
    matrix = copy.deepcopy(composition)
    
    #expansion step
    for chr in matrix:
        for sample in samplelist:
            temp_result = []
            for section in matrix[chr][sample]:
                start = section[0]
                end = section[1]
                color = colortable[section[2]]
                for x in range(end - start + 1):
                    temp_result.append(ce.toHexColor(color))
            matrix[chr][sample] = temp_result
    
    #output stepmatrix
    with open(outpath, 'w') as output:
        output.writelines(['\t'.join(['Chromosome', 'Position'] + samplelist)])
        output.write('\n')
        for chr in matrix:
            for index, position in enumerate(zip(*[matrix[chr][sample] for sample in samplelist])):
                num = str(index * blocksize)
                output.write('\t'.join([chr, num] + list(position)))
                output.write('\n')     
    
    return


'''file -> {chr:[density of blocks]]}
used to load the densities file'''
def loadDensity(infile):
    results = {}
    
    with open(infile, 'r') as input_type:
        data = input_type.read()
        chrs = re.split("@", data)[1:] #first one will be empty for the first symbol.
        for chr in chrs:
            chrSplit = re.split("\n", chr)[:-1] #the last one will be empty for the ending newline
            chrName = chrSplit.pop(0)
            results[chrName] = []
            for line in chrSplit:
                lineSplit = re.split("\t", line)
                results[chrName].append(int(lineSplit[1]))
    return results

#used to load density files with individual informations (from alignment files, for example)
#into the tree format for comparision
#currently hacked to make everything 100. Fix when appropriate!
def loadMultiDensity(infile):
    
    results = {}
    
    with open(infile, 'r') as input_type:
        data = input_type.read()
        lines = re.split("\n", data)
        columns = {}
        for index, element in enumerate(re.split("\t", lines[0])[2:]):
            columns[index] = element
        currChr = ""
        
        for line in lines[1:-1]:
            lineSplit = re.split("\t", line)
            chr = lineSplit[0]
            if chr != currChr:
                results[chr] = chrList = []
                currChr = chr
            temp = {}
            for index, element in enumerate(lineSplit[2:]):
#                 temp[columns[index]] = int(element)
#             chrList.append(temp)

                temp[columns[index]] = 500 #this is the hack
            for x in range(10):   #also this is because the current density file looks every 100K as opposed to 10K like I wanted
                chrList.append(temp)
            
    return results
            
        
            
        

if __name__ == '__main__':
    all_strains = 'ABCDEFGHIJKLMN'
    a,b,c,d,e,f,g,h,i,j,k,l,m,n = tuple(all_strains)
    global groups
    groups = {'G1':[a,b,c,d], 'G2':[e,f,g,h,i], 'G3':[j,k], 'G4':[l,m,n]}
        
    block1 = [[a,b,c,d], [e,f,g,h,i], [j,k], [l,m,n]]
    block2 = [[a,b,c,d], [e,f,g,h,i,j], [l,m,n,k]]
    block3 = [[a,b], [e,f,g,h,i,c,d], [j,k], [l,m,n]]
    block4 = [[a,b,c,d], [e,f,g,h,i], [j,k,l,m,n]]
    
    typeI = [block1]*4 + [block2]*4 + [block3]*4 + [block4]*4
    typeII = [block1]*2 + [block2]*10 + [block3]*2 + [block4]*7
    
    data_tree = {'CHRI': typeI, 'CHRII': typeII}
    result = {x:findContigsComposition(x, data_tree) for x in all_strains}
    print('\n'.join(map(lambda x: '{0}:{1}'.format(x[0], x[1]), sorted(result.items()))))
    
