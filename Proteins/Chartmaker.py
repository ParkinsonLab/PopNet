'''
Created on Jul 20, 2014

@author: Javi
'''
import re
'''makes some charts based on the clusters. borrows the loadclusters method from
SRSCytoscape. really should extract that.'''

'''(list, 2D list) -> matrix
makes a chart that describes which SRS is present in which strain'''
def distributionChart(translatedMatrix, strainList):
    resultMatrix = []
    clusters = translatedMatrix
    iterated = []
    
    while clusters:      
        for cluster in clusters:
            #presort into both separate names
            separated = {}
            for element in cluster:
                if element[1] in separated:
                    separated[element[1]].append(element[0])
                else:
                    separated[element[1]] = [element[0]]
                    
            for name, strains in separated:            
                results = {}

                for strain in strainList:
                    results[strain] = 0

                for strain in strains.values():
                    results[strain] += 1
                    
                resultMatrix.append((name, results))
    #add the .1 and stuff
    sortedResults = sorted(resultMatrix, key=lambda x: x[0])
    return (sortedResults, strainList)
                        
''' returns the matrix (with IDs turned into tuples) and the strain list
make a 2D matrix of tuples, stating the strain and the SRS name'''
def translateMatrix(data, composition):
    strainList = []
    translatedMatrix = []
 
    for cluster in data:
        current = []
        translatedMatrix.append(current)
        for element in cluster:
            strain = re.split("_", element)[0]
            if strain not in strainList:
                strainList.append(strain)
            name = composition[element][0]
            translatedMatrix.append((strain, name))
    
    return (translatedMatrix, strainList)

def printChart(matrix, strainList, outfile):
    strainList.insert("", 0)
    with open(outfile, "w") as output:
        previous = ""
        count = 0
        output.write(",".join(strainList) + "\n")
        for item in strainList:
            header = item[0]
            if header == previous:
                header += ".{0}".format(str(count+1))
                count += 1
            else:
                previous = item[0]
                count = 0
            
            body = item[1].insert(0,header)
            output.write(",".join([str(x) for x in body]) + "\n")
            
            
        