'''
Created on Jun 13, 2014

@author: javi
'''
import Proteins.HMMParser as hmp
import Proteins.DomainFamily as dmf
import Proteins.SequenceSelect as ss
import Proteins.SRSCytoscape as srsce
import Proteins.Chartmaker as cm
import os
import re


if __name__ == '__main__':
#     #Job: merging multiple hmm families. This script specifically deals with the file/folder structure.
#     directory = "/data/javi/Proteins/results"
#     resultPath = "/data/javi/Proteins/HMMResolved"
#     try:
#         os.mkdir(resultPath)
#     except:
#         pass
#          
#     print("reading data")
#     data = {}
#     for dir in os.walk(directory).next()[1]:
#         family = dir
#         path = "/".join([directory, dir])
#         files = [ x for x in os.listdir(path) if os.path.isfile("/".join([path,x]))]
#         for file in files:
#             strain = re.split("\.",file)[0]
#             if strain not in data:
#                 data[strain] = []
#             data[strain].append(hmp.read("/".join([path, file]), family))
#          
#     print("analyzing")
#     results = {}
#     for strain, info in data.items():
#         print("resolving strain: " + strain)
#         results[strain] = dmf.choose(info)
#         dmf.printTree(results[strain], "{0}/{1}.hmmr".format(resultPath, strain))
#             
#         
#     print("HMM reorganiation completed.")
#     
#     # Job: selecting all the found sequences from fasta file.
#     hmmrDirectory = "/data/javi/Proteins/HMMResolved"
#     fastaDirectory = "/data/javi/Proteins"
#     outputDirectory = "/data/javi/Proteins/SelectedFastas"
#         
#     hmmrFiles = [ x for x in os.listdir(hmmrDirectory) if os.path.isfile("/".join([hmmrDirectory, x]))]
#     for hmmr in hmmrFiles:
#         print("selecting proteins for {0}".format(hmmr))
#         strain = re.match("^(.+?).hmmr", hmmr).group(1)
#         hmmrpath = "/".join([hmmrDirectory, hmmr])
#         fastapath = "/".join([fastaDirectory, strain + ".aa.fsa"])
#         outpath = "/".join([outputDirectory, strain + ".hits.fasta"])
#         ss.select(hmmrpath, fastapath, outpath, strain)
#     print("selection completed")


#     #Job: filtering out hammodia
#     filepath="/data/javi/Proteins/orthoMCL/final.txt"
#     outpath="/data/javi/Proteins/orthoMCL/noHam.txt"
#     ss.filter(filepath, outpath)
#     print("done filtering")

#     #Job: changing the names.
#     filepath="/data/javi/Proteins/orthoMCL/groups.txt"
#     outpath="/data/javi/Proteins/orthoMCL/named.csv"
#     hmmrDirectory="/data/javi/Proteins/HMMResolved"
#     ss.toName(filepath, hmmrDirectory, outpath)

#     #Job: cytoscape encoding
#     filepath="/data/javi/Proteins/orthoMCL/groups.txt"
#     outpath="/data/javi/Proteins/SRSCytoscape.xgmml"
#     hmmDirectory="/data/javi/Proteins/HMMResolved"
#     clusters = srsce.loadClusters(filepath)
#     cytotext = srsce.encode(clusters, hmmDirectory, 'SRS Proteins')
#     with open(outpath, 'w') as output:
#         output.write(cytotext)
#     print("cytoscape file completed.")

#     #job: Make distrubution chart
#     filepath="/data/javi/Proteins/orthoMCL/groups.txt"
#     outpath="/data/javi/Proteins/chart.csv"
#     hmmrDirectory="/data/javi/Proteins/HMMResolved"
#     clusters = srsce.loadClusters(filepath)
#     domains = srsce.getNodeDomains(hmmrDirectory)
#     translated = cm.translateMatrix(clusters, domains)
#     data = cm.distributionChart(translated[0], translated[1])
#     cm.printDistributionChart(data[0], data[1], outpath)
    
    #job: Make domainsstats chart
    outpath="/data/javi/Proteins/domainStat.csv"
    hmmrDirectory="/data/javi/Proteins/HMMResolved"
    composition = srsce.getNodeDomains(hmmrDirectory)
    data = cm.domainStats(composition)
    cm.printDomainStats(data, outpath)
    
    