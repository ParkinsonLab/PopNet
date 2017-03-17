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
def translate(input, **kwargs):
    if re.search("M76611", input):
        return "Pf3D7_ChrM"
    
    mode = kwargs['mode']
    
    try:
        if mode == 'plasmodium':
        #plasmodium way
            inpattern = "(.+?)[_](.+?)[_].*"
            mobject = re.search(inpattern, input.upper()) 
            strain = mobject.group(1)
            number = int(mobject.group(2))
            roman = "CHR" + toRoman(number)
        elif mode == 'yeast':
        #yeast way
            inpattern = 'REF[\|](.+?)_(.+?)[\|]'
            mobject = re.search(inpattern, input.upper())
            strain = mobject.group(1)
            roman = mobject.group(2)
        elif mode == 'toxoplasma':
        #Toxo
            inpattern = 'TGME49_(.+)$'
            mobject = re.search(inpattern, input.upper())
            strain = 'TGME49'
            roman = mobject.group(1)
        elif mode == 'strep':
            strain = 'STREP'
            roman = 'CHRI'
        else:
            raise RuntimeError('unknown mode specified')
            import sys
            sys.exit()
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

def reverseTranslate(numStr):
    '''(str) -> str
    gives a partial, but correctly formatted
    chromosome name for matching purposes'''
    
    returnstring = 'CHR'
    try:
        lastchar = ''
        number = int(numStr)
    except TypeError:
        lastchar = numStr[-1]
        number = int(numStr[:-1])
    
    return returnstring + toRoman(number) + lastchar
    
def isEuk(organism):
    euks = ['toxoplasma', 'plasmodium', 'yeast']
    if organism in euks: return True
    else: return False
          
#testing purposes only
if __name__ == '__main__':
    chr = "TGME49_chrVIIA"
    print(translate(chr, mode='toxoplasma'))
