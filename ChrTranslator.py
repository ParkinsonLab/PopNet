'''
Created on Oct 6, 2014

@author: javi

This is a translator for the chromosome names used in vcf files. Basically, different reference genomes have different ways to name their
chromosomes. Hence this serves as a way to give them unified names in the format of 

Strain_ChrXX

as things stand we might have to just rewrite this for every new dataset. 
'''

import re

#for the plasmodium one
def translate(input):
    if re.search("M76611", input):
        return "Pf3D7_ChrM"
#     inpattern = "(.+?)[_](.+?)[_].*" #for plasmodium
#     inpattern = 'REF[\|](.+?)_(.+?)[\|]'
    inpattern = 'TGME49_(.+)$'
    mobject = re.search(inpattern, input.upper())

    
    try:
#         #plasmodium way
#         strain = mobject.group(1)
#         number = int(mobject.group(2))
#         roman = "CHR" + toRoman(number)

# #         yeast way
#         strain = mobject.group(1)
#         roman = mobject.group(2)

        #Toxo
        strain = 'TGME49'
        roman = mobject.group(1)
    except:
        print("Illegal Chr Name @ {}".format(input))
        return "unknown"
        
    returnstring = "{0}_{1}".format(strain, roman)
    return returnstring.upper()

def toRoman(number):
    returnstring=''
    table=[['X',10],['IX',9],['V',5],['IV',4],['I',1]]

    for pair in table:

        while number-pair[1]>=0:

            number-=pair[1]
            returnstring+=pair[0]

    return returnstring

#testing purposes only
if __name__ == '__main__':
    chr = "TGME49_chrXII"
    print(translate(chr))