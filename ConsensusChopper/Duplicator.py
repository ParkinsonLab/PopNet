'''
Created on Oct 9, 2013

@author: javi
'''
import os

if __name__ == '__main__':
       #modify these as needed
    directory = input("Please specify working directory, in full \n")
    identifier = input("Please input sequence identifier, not including the @ symbol \n")
    length = eval(input("Please specify the length of moving window \n"))
    overlap = eval(input("Please specify the length of overlap \n"))
    
    #Do not modify
    os.chdir(directory)  
     
    if not (os.path.isdir(directory + "/Chopped")):
        os.mkdir("Chopped")
    
    chopdirectory = directory + "/Chopped"
        
    log = open(chopdirectory + "/log.txt", "w")
    log.write("Run inputs are: \n%s\n%s\n%s\n%s\n" % (directory, identifier, length, overlap))
    
    onlyfiles = [ f for f in listdir(directory) if isfile(join(directory,f)) ]
    for f in onlyfiles:
        print("\nProcessing %s ..." % f)
        log.write("\nProcessing %s ... " % f)
        
        try: 
            namesplit = f.split(".")
            if namesplit[len(namesplit)-1] == "fq":
                chop(directory, chopdirectory, f, identifier, length, overlap, log)
            else:
                log.write("\nerror, %s is not a fq file" % f)
        except:
            log.write("\nerror, %s is not a fq file" % f)