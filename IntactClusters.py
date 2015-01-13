'''
Created on Apr 2, 2014

@author: javi

This script is used to produce a particular view on the Java based Chr view, to determine whether a group of strains cluster together
over the entire genome. This is both a discovery (when you have some leads) and a validation tool (when you think something things
should be grouped together). The results produced here was not used for anything, and this lead may be discontinued for all I know. 
'''



#depends on MCLCounter!
'''(list of names, list of lines) -> int
checks whether the designated cluster is intact. 
for each line in clusters, check whether all element of target is in that line.'''
def isIntact(target, clusters):
    #criteria is 75% of the members must cluster together
    minHits = int(0.6 * len(target))
    for line in clusters:
        hits = 0
        for member in target:
            if member in line:
                hits += 1
        if hits > minHits: return 1
    return 0

'''(dict of name:list, list of lines) -> dict of name:int
calculates the scores for the current block. Looks through each target in the targets, 
and then call isIntact for each. Puts the result, 1 or 0 for y/n, in the results dict.
WARNING: FUNCTION ONLY WORKS IF CRITERIA IS > 50%'''     
def calculateBlock(targets, clusters):
    results = {}
    for targetName, target in list(targets.items()):
        results[targetName] = isIntact(target, clusters)
        if results[targetName] > 1: print(('ERROR: TARGET INTACT IN MORE THAN 1 CLUSTER. \n{0}\n{1}'.format(targetName, clusters)))
    return results

'''(dict of list of dicts, filename) -> None, writes to file
writes the results of the 'process' function to the newly created outfile in the same
format as similarity.py. This is to be read by the java SNPViewer to generate the
chr view.'''   
def write(resultsTree, outfile):
    #something to put on the first part of each line
    selfName = "GENOME"
    
    with open(outfile, "wb") as output:
        #write the hit list as the keys within the first block 
        output.write('$\n' + '\n'.join(list(list(resultsTree.values())[0][0].keys())) + '\n$\n')
        
        #write the rest
        for chrName, chr in list(resultsTree.items()):
            output.write('%s\n' % chrName)
            count = 0
            for block in chr:
                output.write('#%d\n'%count)
                output.write('\n'.join(['%s - %s: %f'%(selfName, itemName, value) for itemName, value in list(block.items())]))
                output.write('\n')
                count+=1

'''(None) -> List
placeholder function for generating the target list. Right now I'm just hardcoding it in'''
def generateTargets():
    targets = {}
    targets['ME'] = ['B73', 'ME49', 'PRU', 'B41', 'ARI', 'RAY']
    targets['RH'] = ['CAST', 'TgCkCr1', 'RH-88', 'RH-JSR', 'GT1', 'TgCkCr10', 'TgCkBr141']
    targets['GUY'] = ['BRC_TgH_18002_GUY-KOE', 'GUY-2004-ABE', 'BRC_TgH_18003_GUY-MAT', 'BRC_TgH_18009', 'BRC_TgH_18021', 'GUY-2003-MEL', 'BRC_TgH_18001_GUY-DOS']
    return targets

'''(list of lists, dict of lists, file name) -> None (writes to file)
takes a clusterTree as processed by a reader, calls calculate block
on each cluster, and writes the result to the outfile by write()
The result tree format will be:
{chrName1:['target1':value1, 'target2', value2...], chrname2...}
'''    
def process(clusterTree, outfile):
    resultsTree = {}
    targets = generateTargets()
    for name, chr in list(clusterTree.items()):
        currChr = []
        resultsTree[name] = currChr
        for block in chr:
            blockResults = calculateBlock(targets, block)
            currChr.append(blockResults)
    write(resultsTree, outfile)

if __name__ == '__main__':
    import MCLCounter
    infile = '/data/javi/Toxo/64Genomes/Counting/persistentResult.txt'
    tabfile = '/data/javi/Toxo/64Genomes/Counting/persistentMatrix.tab'
    outfile = '/data/javi/Toxo/64Genomes/Counting/intact.sim'
    clusterTree = MCLCounter.toMatrix(MCLCounter.loadClusters(infile, tabfile)[0])
    process(clusterTree, outfile)
    print("End of Intact Clusters Script")
    
    
    
    