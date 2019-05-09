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
offset = 0

def sortNames(chr_names):
    '''
    main sort function
    '''

    try:
        return sorted(chr_names, key=lambda x: getValue(x))
    except ValueError:
        return chr_names

def getValue(name):
    global offset
    tmp_offset = 0
    name = name.upper()
    pattern = '.+_(CHR[XIMV]+[ABCDE]?|[XIMV]+|[0-9]+[ABCDE]?)(_.+|$)'
    suffix = 'ABCDE'
    
    section = re.search(pattern, name)

    if section is None:
        raise ValueError
    else:
        section = section.group(1)

    #Removes CHR if applicable
    if section.startswith('CHR'):
        section = section[3:]
    
    #Removes suffix if applicable
    if section[-1] in suffix:
        tmp_offset = suffix.index(section[-1])
        section = section[:-1]
    
    #Try it as numeric, go for roman if it doesn't work.
    try:
        val = int(section) + offset + tmp_offset
    except ValueError:
        val = romanToNum(section) + offset + tmp_offset

    offset += min(1, tmp_offset)
    return val

def romanToNum(roman):
    result = 0
    index = 0
    for numeral, integer in romanNumeralMap:
        while roman[index:index+len(numeral)] == numeral:
            result += integer
            index += len(numeral)
    return result