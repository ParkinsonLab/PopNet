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
    
      
    for cluster in clusters:
        #presort into both separate names
        separated = {}
        for element in cluster:
            if element[1] in separated:
                separated[element[1]].append(element[0])
            else:
                separated[element[1]] = [element[0]]
                
        for name, strains in separated.items():            
            results = {}

            for strain in strainList:
                results[strain] = 0

            for strain in strains:
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
        
        for element in cluster:
            
            nameSplit = re.split("_", composition[element][0])
            strain = nameSplit[0]
            name = nameSplit[1]
            if strain not in strainList:
                strainList.append(strain)
            current.append((strain, name))
        translatedMatrix.append(current)
        
    return (translatedMatrix, strainList)

def printDistributionChart(matrix, strainList, outfile):
    strainList.insert(0,"")
    with open(outfile, "w") as output:
        previous = ""
        count = 0
        output.write(",".join(sorted(strainList)) + "\n")
        for item in matrix:
            header = item[0]
            if header == previous:
                header += ".{0}".format(str(count+1))
                count += 1
            else:
                previous = item[0]
                count = 0
            
            body = header
            for entry in sorted(item[1].items(), key=lambda x: x[0]):
                body+=",{0}".format(str(entry[1]))
            output.write(body + "\n")
            
def domainStats(composition):
    results = {}
    for ID, item in composition.items():
        strain = re.split("_", ID)[0]
        domains = item[1]
        if strain not in results:
            results[strain] = {}
        for domain in domains:
            if domain not in results[strain]:
                results[strain][domain] = 0
            results[strain][domain] += 1
    return results

def printDomainStats(data, outfile):
    with open(outfile, "w") as output:
        famList=["fam_{0}".format(str(x)) for x in range(1, 9)]
        header="," + ",".join(famList) + ",Total\n"
        output.write(header + "\n")
        
        for strain, domainInfo in sorted(data.items()):
            output.write(strain + ",")
            info = []
            for fam in sorted(famList):
                if fam in domainInfo:
                    info.append(str(domainInfo[fam]))
                else:
                    info.append(str(0))
            total = sum([int(x) for x in info])
            info.append(str(total))
            output.write(",".join(info) + "\n")
        
        
        