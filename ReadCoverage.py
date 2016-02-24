'''
Created on Oct 23, 2013

@author: javi
'''
import re
import os
from os import listdir
from os.path import isfile, join
import numpy
import matplotlib
import traceback
import copy
import ChrTranslator as ct

def calcReadCoverage(file, output, minOutput, mode):
    #dataparsing
    outputName = output.name
    noHeader = re.findall("(?m)^([^@].*)\n", file.read()) #removed header, divided into lines.
    file.close()
    
    dataTree = {}
    avgTree = {}
    #Basically the idea is that, each time a read covers a position, you add 1 to that position's read coverage. 
    for x in range(0, len(noHeader)):
        currLine = noHeader[x]
        lineSegment = re.split("\t", currLine)
        #the chr name is translated using vct
        info = [ct.translate(lineSegment[2], mode=mode), int(lineSegment[3]), len(lineSegment[9])] #Info is: chromosome, leftmost position, length of the match.
        
        if not info[0] in dataTree:
            dataTree[info[0]] = {}
            avgTree[info[0]] = {}
        currentBranch = dataTree[info[0]]
        avgBranch = avgTree[info[0]]
        
        for x in range(info[1], info[1] + info[2]):
            if not x in currentBranch:
                currentBranch[x] = 0
                avgBranch[x] = 0
            currentBranch[x] += 1
        
    for chrName in dataTree:
        if re.search("(?i)chr", chrName):
            
            #smooth it out by moving average
            #record
            
            output.write("@%s\n"%chrName)
            minOutput.write("@%s\n"%chrName)
            
            
            currBranch = dataTree[chrName]
            currAvgBranch = avgTree[chrName]
            highest = sorted(currBranch.keys())[-1] + 1 #correction of 1 due to the behavior of range()
            for x in range(0, highest):
                upperBound = min([x+100, highest])
                lowerBound = max([x-100, 0])
                tempValue = 0.0
                for z in range(lowerBound, upperBound):                
                    if z in currBranch:
                        tempValue += currBranch[z]
                value = tempValue / (upperBound - lowerBound)
                currAvgBranch[x] = value
                 
            for x in range(0, highest):
                try:
                    val = currBranch[x]
                except KeyError:
                    val = 0
                
                try:
                    avg = currAvgBranch[x]
                except KeyError:
                    print("no average for %d computed!"%x)
                    avg = 0
                
                output.write("%s\t%d\t%f\n"%(chrName, x, avg))
                if val < 20 or val > (avg*3):    
                    minOutput.write("%s\t%d\t%f\n"%(chrName, x, val))
           
    output.close()
    minOutput.close()
    
    
def graphResults(results):
    
    pdf = PdfPages('%s_Plots.pdf'%results.name)
    sizedPdf = PdfPages('%s_300Plots.pdf'%results.name)
    allData = results.read()
    chrList = re.findall("(?m)^[@](.*?)\n", allData)
    count = 1 
    for chrName in chrList:
        chrData = re.findall("(?m) ^%s\t([0-9.]*?)\t([0-9.]*?)\n"%chrName, allData)
        positions = [x[0] for x in chrData]
        values = [x[1] for x in chrData]
        plt.figure(count)
        plt.title("Read Coverage for %s"%chrName)
        plt.plot(positions, values)
        pdf.savefig() 
        plt.ylim(ymax = 300, ymin = 0)
        sizedPdf.savefig()
        count += 1
     
    pdf.close()
    sizedPdf.close()
    plt.close("all")
    results.close()
        

if __name__ == '__main__':
    
    #modify these as needed
#    directory = raw_input("Please specify working directory, in full \n")
    directory = "/data/new/javi/yeast/scinetDownload/sams"
    mode = 'yeast'
    
    #Do not modify
    os.chdir(directory)  
    
    onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
    matplotlib.use("pdf")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    
    for f in onlyfiles:
        print("\nProcessing %s ..." % f)        
        try: 
            namesplit = f.split(".")
            if namesplit[len(namesplit)-1] == "sam":
                outputName = "%s_coverage"%namesplit[0]
                calcReadCoverage(open(f, "r"), open("%s.txt"%outputName, "w+"), open("%s.min"%outputName, "w"), mode)
#                 graphResults(open("%s.txt"%outputName, "r"))
        except Exception as e:
            print(traceback.print_exc())
    
    print("end of script")
    

    
   