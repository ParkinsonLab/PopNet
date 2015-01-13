'''
Created on Aug 7, 2014

@author: javi
'''
'''sorts chr names

you need to change the offsets for each dataSet!!!!!
'''
import re

romanNumeralMap = (('X',  10),
                   ('IX', 9),
                   ('V',  5),
                   ('IV', 4),
                   ('I',  1),
                   ('M',  0))

def getValue(name):
    
    #current: tg snps
    pattern = '.+_CHR(.+)'
    endPattern = ".*[XIVM]$"
    
    tempValue = re.match(pattern, name.upper()).group(1)
    
    if re.match(endPattern, tempValue):
        lastchar = None
    else:
        lastchar = tempValue[-1]
        tempValue = tempValue[:-1]
    
    value = romanToNum(tempValue)
    offset = getOffSet(value, lastchar)
    
    
    return value + offset

def romanToNum(roman):
    result = 0
    index = 0
    for numeral, integer in romanNumeralMap:
        while roman[index:index+len(numeral)] == numeral:
            result += integer
            index += len(numeral)
    return result

    
def getOffSet(num, letter):
#     #for plasmo and yeast
#     offsetList = {}
#     letterValue = {}
    
    #for Toxo
    offsetList = {1: 1, 7: 1}
    letterValue = {"A": 0, "B": 1, "C": 2}
    
    offset = 0
    for chr, value in offsetList.items():
        if num > chr:
            offset += value
    
    if letter:
        offset += letterValue[letter]
        
    return offset

if __name__ == '__main__':
    a = "@PF3D7_CHRIV"
#    b = "@TGME49_CHRVIIA"
    print(getValue(a))
#    print(getValue(b))