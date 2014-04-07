'''
Created on Feb 17, 2014

@author: javi
'''
import re
import os
import copy
import Chr_NexusEncoder
import NexusEncoder
import numpy
import matplotlib

#path = results.analyzed
#tabpath = the tab file

def loadClusters(path, tabpath):
    print "loading clusters..."
    with open(path, "r") as source, open(tabpath, "r") as tab:
        #create the sampleList for making the matrix
        sampleList = []
        rawSampleList = filter(None, re.split("\n",tab.read()))
        for element in rawSampleList:
            ts = re.split(" ", element)
            sampleList.append(ts[1])
            
        rawData = source.read()
        results = {}
        chrs = re.findall("(?s)(@.+?)\n(.+?)\n(?=$|@)", rawData)
        for chr in chrs:
            chrBranch = []
            results[chr[0]] = chrBranch
            segments = filter(None, re.split("\n\n", chr[1])) #each "matrix" is separated by a blank line.
            for segment in segments:
                chrBranch.append(segment)
    print "done"    
    return (results, sampleList)



def distribution(treeTuple, outputpath):
    print "calculating distribution..."
    results = {}
    data = treeTuple[0]
    sampleList = treeTuple[1]
    
    for chrName, chr in data.items():
        chrBranch = []
        results[chrName] = chrBranch
        for segment in chr:
            split = filter(None, re.split("\n", segment))
            position = int(split[0][1:])
            size = len(split) - 1
            chrBranch.append((position, size))
    
#     #Write it all out
#     with open(outputpath, "w") as outputFile:
#         for name, chr in results.items():
#             outputFile.write("%s\n"%name)
#             for line in chr:
#                 outputFile.write("%i\t%i\n" % (line[0], line[1]))
#     print "done"
    return results

def graph_distribution(results):
    pdf = PdfPages('Distribution_plots.pdf')
#    sizedPdf = PdfPages('%s_300Plots.pdf')

    count = 0
    for chrName, chr in results.items():
        
        positions = [x[0] for x in chr]
        values = [x[1] for x in chr]
        plt.figure(count)
        plt.title("Cluster Size Distribution for %s"%chrName)
        plt.plot(positions, values)
        plt.ylim(ymax = 10, ymin = 0)
        pdf.savefig() 
#        plt.ylim(ymax = 300, ymin = 0)
#        sizedPdf.savefig()
        count += 1
     
    pdf.close()
#    sizedPdf.close()
    plt.close("all")
    print "done"
                
def findDifferences(treeTuple, outpath, blackList):
    '''(results, sampleList), string -> (results, sampleList)
    finds positions where the clustering pattern changes by a lot. Possibly need to be calibrated.
    '''
    print "finding major differences"
    with open(outpath, "w") as output:
        data = treeTuple[0]
        
        resultMatrix = {}
        sampleList = treeTuple[1]
        for x in sampleList:
            resultMatrix[x] = 0
        
        results = {}
        
        for name, chr in data.items():
            chrResultBranch = {}
            results[name] = chrResultBranch
            prevCluster = filter(None, re.split("\n", chr[0])[1:])
            for cluster in chr[1:]:
                clusterSplit = filter(None, re.split("\n", cluster))
                score = diffScore(clusterSplit[1:], prevCluster)
                chrResultBranch[int(clusterSplit[0][1:])] = len(score)
                
                #score matrix
                for element in score:
                    if element == '#1':
                        print "bad"
                    resultMatrix[element] += 1
                
#toggle for which type you want - comparing against pre or compare against static                        
                prevCluster = clusterSplit[1:]
        
        return (results, resultMatrix)
        
                

# def diffScore(lines, prevLines): #clusterSplit of the current cluster, clusterSplit of the previous cluster
#     '''list, list -> int
#     calculates a measure of how different the two clustering patterns are.
#     It counts how many different pair-associations exist between the
#     two patterns. May be normalized later.
#     As things stand, this is biased towards migration between large clusters. 
#     '''
#     pairs = []
#     results = []
# 
#     
#     for line in lines:
#         pairs.extend(returnAllPairs(filter(None, re.split("\t", line))))
#     
#     for line in prevLines:
# #        print line
#         pairs.extend(returnAllPairs(filter(None, re.split("\t", line))))
#     
#     for pair in pairs:
#         if pairs.count(pair) == 1:
#             results.append(pair)
#             print pair
#     print "--"
#     return results
            
def diffScore(lines, prevLines):
    a = [re.split("\t", line) for line in lines]
    b = [re.split("\t", line) for line in prevLines]
    results = set([])
    matched = matchClusters(a, b)
    for pair in matched:
        for element in pair[1]:
            if element not in pair[0]:
                results.add(element)
                
    
    return results
    
    
def matchClusters(lines, prevLines):
    #the lines here are split. Returns a list of tuples with matching groups.
    results = []
    matched = []
    for line in lines:
        max = 0
        match = []
        for prevLine in prevLines:
            curLength = len(set(line).intersection(set(prevLine)))
            if curLength > max and prevLine not in matched:
                match = prevLine
                max = curLength
        results.append((line, match))
        matched.append(match)
    
    for line in prevLines:
        if line not in matched:
            results.append(([], line))
    
    return results  
        
def graphResults(name, results):
    print "graphing..."
    pdf = PdfPages('%s.pdf'%name)
#    sizedPdf = PdfPages('%s_300Plots.pdf'%results.name)
    count = 1 
    for chrName, chr in results.items():
        positions = [x[0] for x in sorted(chr.items())]
        values = [x[1] for x in sorted(chr.items())]
        plt.figure(count)
        plt.title("Difference Scores for %s"%chrName)
        plt.plot(positions, values)
        pdf.savefig() 
#        plt.ylim(ymax = 300, ymin = 0)
#        sizedPdf.savefig()
        count += 1
     
    pdf.close()
#    sizedPdf.close()
    plt.close("all")
                
def count(treeTuple, outputPath):
    print "counting..."
    with open(outputPath, "w") as output:
        #create the sampleList for making the matrix
        sampleList = treeTuple[1]
        data = treeTuple[0]
        results = {}
        baseResultMatrix = {}
        for a in sampleList:
            baseResultMatrix[a] = {}
            for b in sampleList:
                baseResultMatrix[a][b] = 0
        #Organizing the data and input into matrix.         
        for name, chr in data.items():
            #reset
            resultMatrix = copy.deepcopy(baseResultMatrix)
            for segment in chr:
                split = filter(None, re.split("\n", segment))
                for line in split[1:]:
                    pairs = returnAllPairs(filter(None, re.split("\t", line)))
                    for x, y in pairs:
                        resultMatrix[x][y] += 1
                        if x != y:
                            resultMatrix[y][x] += 1
            results[name] = resultMatrix
        #return/outputs the result
            output.write(chr[0])
            output.write("\n")
            output.write(repr(sampleList))
            output.write("\n")
            output.write(repr(resultMatrix))
            output.write("\n")
#             for x in resultMatrix:
#                 output.write(str(resultMatrix[x]))
#                 output.write("\n")
        print("done")
        return (results, sampleList)
    
def returnAllPairs(list):
    if len(list) == 1:
        return [(list[0], list[0])]
    result = []
    for e in list[1:]:
        result.append((list[0], e))
    result += [(list[0], list[0])]
    return result + returnAllPairs(list[1:])

def aggregate(tree):
    '''combines all the information from all the chromosomes'''
    print "aggregating results..."
    chrs = tree.items()
    for chr in chrs[1:]:
        for sampa in chr[1]:
            for sampb in chr[1][sampa]:
                chrs[0][1][sampa][sampb] += chr[1][sampa][sampb] 
                print chrs[0][1][sampa][sampb] > chr[1][sampa][sampb]
    result = {"Genome": chrs[0][1]}
    print result
    print "done"
    return result
                        
def summary_distribution(results):
    summary = {}
    for name, chr in results.items():
        for dp in chr:
            num_of_clusters = dp[1]
            if num_of_clusters not in summary:
                summary[num_of_clusters] = 0
            summary[num_of_clusters] += 1
    print sorted(summary.items())
            
def printMatrix(matrix, outfile):
    '''[[e]] -> None
    [x[1] for x in matrix.items()]
    writes a representation of the matrix into the output file'''
    with open(outfile, "w") as output:
        header = [x for x in matrix]
#        row_format ="{:>30}" * (len(header) + 1)
        row_format = "{:>30}" * 2 + "\n"
#        output.write(row_format.format("", *header) + "\n")
        for name, info in matrix.items():
#            output.write(row_format.format(name, *[x[1] for x in info.items()]) + "\n")
            output.write(row_format.format(name, str(info)))
            
def printDiff(filename, results):
    with open(filename, "w") as out:
        for key, item in results.items():
            out.write(key + "\n")
            for key, item in item.items():
                out.write("{0} : {1}\n".format(key, item))
                                
'''(chrName:list of strings) -> chrName: List of Matrices
Breaks the each "clustering pattern" into a 2D list of elements.'''
def toMatrix(clusters):
    results = {}
    for chrName, chr in clusters.items():
        results[chrName] = [[y.split('\t') for y in x.split('\n')[1:]] for x in chr]
    return results
        
                    
if __name__ == '__main__':    
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
#     #Test Data
#     path = "/home/javi/testzone/Griggs Stuff/persistentResults.txt"
#     tabpath = "/home/javi/testzone/Griggs Stuff/persistentMatrix.tab"
#     outputpath = "/home/javi/testzone/Griggs Stuff/Diff.txt"
#     os.chdir("/home/javi/testzone/Griggs Stuff")
#     treeTuple = loadClusters(path, tabpath)
            
    #To Load the whole genome. 
    path = "/data/javi/Toxo/64Genomes/Counting/persistentResult.txt"
    tabpath = "/data/javi/Toxo/64Genomes/Counting/persistentMatrix.tab"
    outputpath = "/data/javi/Toxo/64Genomes/Counting/Count.txt"
    os.chdir("/data/javi/Toxo/64Genomes/Counting")

    diffPath = "/data/javi/Toxo/64Genomes/Counting/diffs/diff.txt"
#    treeTuple = count(loadClusters(path, tabpath), outputpath)
    treeTuple = loadClusters(path, tabpath)
    
#     #To Encode to Nexus
#     wholegenome = aggregate(treeTuple[0])
#     Chr_NexusEncoder.clusterMatrixOutput(wholegenome, treeTuple[1])
#     Chr_NexusEncoder.nexusMatrixOutput(wholegenome, treeTuple[1])
#     NexusEncoder.nexusOutput(wholegenome)
    
#     #To calculate distribution
#     path = "/data/javi/Toxo/64Genomes/Counting/persistentResult.txt"
#     tabpath = "/data/javi/Toxo/64Genomes/Counting/persistentMatrix.tab"
#     outputpath = "/data/javi/Toxo/64Genomes/Counting/distribution.txt"
#     os.chdir("/data/javi/Toxo/64Genomes/Counting")
#     summary_distribution(distribution(loadClusters(path, tabpath), outputpath))
    
    # To calculate diff score
    results = findDifferences(treeTuple, outputpath, [])
#    print results
    graphResults("Differences", results[0])
    printDiff(diffPath, results[0])
    printMatrix(results[1], "diffmatrix.txt")
    print results
    print "end of MCLCounter script"