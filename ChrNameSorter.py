'''
Created on Aug 7, 2014

@author: javi
'''
'''sorts chr names
'''
import re
import string

romanNumeralMap = (('X',  10),
                   ('IX', 9),
                   ('V',  5),
                   ('IV', 4),
                   ('I',  1),
                   ('M',  0))

def sortNames(chr_names):
    '''
    main sort function
    '''

    try:
        return sorted(chr_names, key=lambda x: getValue(x))
    except ValueError:
        return chr_names

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


if __name__ == '__main__':
    #testing
    print(getValue('Pf3D7_06_v3'))