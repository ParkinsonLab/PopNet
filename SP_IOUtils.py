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

def listJoin(list):
    if len(list) < 2:
        return list
    return "\n".join(list)

def dumpDict(input_dict, outpath):
    with open(outpath, "wb") as output:
        output.write(bytes(recursiveDictJoin(input_dict)).encode('utf-8'))


def recursiveDictJoin(input_dict):
    import types
    if not isinstance(input_dict, types.DictType):
        return listJoin(input_dict)
    return "\n".join(["\n{0}:{{\n{1}}}".format(x, recursiveDictJoin(y)) for x, y in sorted(input_dict.items())])

        
def dumpDataTree(dataTree, outpath):
    with open(outpath, "wb") as output:
        for name, chr in sorted(dataTree.items()):
            output.write(bytes(name + ":\n" + recursiveDictJoin(chr) + "\n").encode('utf-8'))