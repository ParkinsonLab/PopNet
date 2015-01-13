'''
Created on Oct 20, 2014

@author: javi
'''

import re
import copy
import os
from collections import Counter as counter


def scan(dataTree):
    print("scanning for drift")
    count = 0
    driftcount = 0
    for name, chr in dataTree.items():
        for posNum, position in chr.items():
            count+= 1
            if len(position)<=2:
                del chr[posNum]
                driftcount+= 1
    percent = float(driftcount) / float(count) * 100
    print("{0} drift positions deleted out of {1}, for a rate of {2} percent".format(driftcount, count, percent))
    
    return dataTree


# '''(list of chars) -> boolean
# see if there are any SNPs which aren't shared by at least 2 other strains.
# Those SNPs are considered drifts and shouldn't be included'''
# def isDrift(list):
#     minStrains = 2
#     uniques = set(list)
#     if len(uniques) > 2:
#         return False
#     else:
#         for element in uniques:
#             if list.count(element) < minStrains:
#                 return True
#         return False

