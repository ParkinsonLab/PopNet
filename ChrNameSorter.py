'''
Created on Aug 7, 2014

@author: javi
'''
'''sorts chr names

you need to change the offsets for each dataSet!!!!!
'''
import re
import string

romanNumeralMap = (('X',  10),
                   ('IX', 9),
                   ('V',  5),
                   ('IV', 4),
                   ('I',  1),
                   ('M',  0))

def getValue(name):

    name = name.upper()
    pattern = '.+_(CHR[XIMV]+[ABCDE]?|[XIMV]+|[0-9]+[A-Z]?)(_.+|$)'
    offset = 0
     
    section = re.search(pattern, name)

    if section is None:
        raise ValueError
    else:
        section = section.group(1)

    #Removes CHR if applicable
    if section.startswith('CHR'):
        section = section[3:]
    
    #Removes suffix if applicable
    if section[-1] in string.ascii_uppercase and section[-1] not in 'IXVM':
        offset = string.ascii_uppercase.index(section[-1]) / 10.
        section = section[:-1]
    
    #Try it as numeric, go for roman if it doesn't work.
    try:
        val = int(section) + offset
    except ValueError:
        val = romanToNum(section) + offset

    return val

def romanToNum(roman):
    result = 0
    index = 0
    for numeral, integer in romanNumeralMap:
        while roman[index:index+len(numeral)] == numeral:
            result += integer
            index += len(numeral)
    return result

    
def getOffSet(num, letter, mode):
#     #for plasmo and yeast
#     offsetList = {}
#     letterValue = {}

    if mode == 'toxoplasma':
        offsetList = {1: 1, 7: 1}
        letterValue = {"A": 0, "B": 1, "C": 2}
    elif mode == 'yeast' or mode == 'plasmodium' or mode == 'strep':
        offsetList = {}
        letterValue = {}
    else:
        raise TypeError('Unknown type for chr translation')
    
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