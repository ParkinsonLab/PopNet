'''
Created on Apr 14, 2014

@author: javi

utilities file containing common functions used to read/write.
'''

'''dumps a list to a readable format'''
def dumpList(list, outpath):
    with open(outpath, "wb") as output:
        output.write(bytes(recursiveListJoin(list)).encode('utf-8'))
    

def recursiveListJoin(list):
    import types
    if not isinstance(list, types.ListType):
        return repr(list)
    return "\n[{0}]".format(", ".join([recursiveListJoin(x) for x in list]))

def dumpDict(dict, outpath):
    with open(outpath, "wb") as output:
        output.write(bytes(recursiveDictJoin(dict)).encode('utf-8'))

def recursiveDictJoin(dict):
    import types
    if not isinstance(dict, types.DictType):
        return recursiveListJoin(dict)
    return "\n".join(["\n{0}:{{\n{1}}}".format(x, recursiveDictJoin(y)) for x, y in dict.items()])
        