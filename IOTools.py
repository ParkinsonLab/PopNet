'''
An aggregate script containing functions for outputting stuff
'''
import os
import numpy as np
from subprocess import check_output
           
def mcl(currString, tabpath, iValue, piValue, raw=False):
    '''
    (dict, list, num, num) -> [[string]]
    helper function for repeatedly running mcl over an array of I and PI values.
    '''
    def mcl_pi():
        return check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_nopi():
        return check_output(["mcl", tempname, "-use-tab", tabpath, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_notab_pi():
        return check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-pi", str(piValue), "-q", "x", "-V", "all", '-te', '1'])
    
    def mcl_notab_nopi():
        return check_output(["mcl", tempname, "-I", str(iValue), "-o", "-", "-q", "x", "-V", "all", '-te', '1'])
    
    tempname = 'temp{}.mci'.format(os.getpid())
    
    with open(tempname, 'w') as temp:
        temp.write(currString)
    
    if tabpath is None and piValue > 0:
        result = mcl_notab_pi()
    elif tabpath is None:
        result = mcl_notab_nopi()
    else:
        if piValue > 0:
            result = mcl_pi()
        else:
            result = mcl_nopi()
            
    result = bytes.decode(result)
    os.remove(tempname)

    if raw:
        return result
    
    else: 
        results = [line.split('\t') for line in result.rstrip('\n').split("\n")]
        return results


def readTab(tab_file):
    '''(str)->list
    loads a tab file for mcl'''
    with open(tab_file, 'r') as input:
        sample_list = [line.split(' ')[1].rstrip('\n') for line in input if line is not '']
    return sample_list

def writeTab(sample_list, tab_file):
    '''
    writes the tab file
    '''
    out_string = '\n'.join(['{0} {1}'.format(i, sample) for i, sample in enumerate(sample_list)])

    with open(tab_file, 'w') as output:
        output.write(out_string)


def writeGroup(names, groups, out_path):
    '''
    list, [list], path -> None
    we're expected synced group with the name of that group
    '''
    if len(names) != len(groups):
        raise ValueError('writeGroup encounted mismatched names and groups')
    
    out_string = '\n'.join(['@{0}\n{1}'.format(name, "\n".join(group)) for name, group in zip(names, groups)])
    
    with open(out_path, 'w') as output:
        output.write(out_string)


def writeClusters(clusters, chrs, clusters_path):
    '''
    [[string]], list -> None
    '''
    
    def oneChr(clusters_in_chr):
        return '\n'.join(['#{0}\n{1}'.format(i, cluster) for i, cluster in enumerate(clusters)])

    if len(clusters) != len(chrs):
        raise ValueError('writeclusters encountered mismatch between clusters and chr list')
    
    output_string = '\n'.join(['{0}\n{1}'.format(chrs[i], oneChr(chr)) for i, chr in clusters])
    
    with open(clusters_path) as output:
        output.write(output_string)

# def clmInfo(currString, tabpath, params):
#     '''runs clminfo in with a range of values, and finds the pattern with the least amount of change over time.'''
    
#     #helper functions
#     ###
    
#     matrixName = 'in.mci'
#     with open(matrixName, 'w') as tmp:
#         tmp.write(currString)
        
#     iVals = stepList(params.getIMin(),params.getIMax(),params.getIStep(), reverse=True)
#     piVals = stepList(params.getPiMin(),params.getPiMax(),params.getPiStep())
    
#     files = []
#     for i in iVals:
#         files += setupClm(currString, i, piVals)
# #     files = reduce(lambda x, y: x + y, files) #make 2d list flat
    
#     result = check_output(['clm', 'info', matrixName] + files, close_fds=True)
#     result = bytes.decode(result)
#     #returns the clm dist (variance of information) but not using that atm.
# #    result = subp.check_output(['clm', 'dist', '-mode', 'sj', '--chain'] + files, close_fds=True)
    
#     for file in files:
#         os.remove(file)
#     os.remove(matrixName)
        
#     return result

'''(dict, list) -> string
builds the matrix (.mci) string to be used in mcl
'''
def buildMatrix(matrix, sample_list):
    
    def buildHeaderRow(ls):
        '''takes the list of sample idx and returns the mclrows bit'''
        ls.insert(0, '')
        ls.append('$')
        print(ls)
        return ' '.join(ls)


    def buildRow(n, ls):
        '''takes a list and its index and returns the formatted row'''
        ls = ['{0}:{1}'.format(i, x) for i, x in enumerate(ls) if i != n]
        ls.append('$')
        return '{0}\t{1}'.format(str(n), ' '.join(ls))

    template = "\
(mclheader\n\
mcltype matrix\n\
dimensions {0}x{0}\n\
)\n\
(mclmatrix\n\
begin\n\
{1}\n\
)"

    if matrix.shape[0] != matrix.shape[1] or matrix.shape[0] != len(sample_list):
        raise ValueError('buildMatrix encountered malformed matrix!')

    # dom_text = buildHeaderRow([str(x) for x in range(len(sample_list))])
    row_text = [buildRow(i, row) for i, row in enumerate(matrix)]

    return template.format(len(sample_list), '\n'.join(row_text))

def writeOverallMatrix(matrix, outfile):
    '''[[e]] -> None
    [x[1] for x in matrix.items()]
    writes a representation of the matrix into the output file'''
    matrix.to_csv(outfile, sep='\t')

def writeTabularPainting(composition, chrs, section_length, sample_list, path):

    def makeLine(data, *args):
        #put the data first
        args = [str(e) for e in args]
        data = [str(e) for e in data]
        return '\t'.join(args + data) + '\n'
    
    def expand(data):
        #now we gotta expand the painting instead so it's a bit of a shame
        result = []
        for e in data:
            if e[1] >= 0:
                result += [e[1]] * e[0]
            else:
                result.append(-1) 
        return result

    chr_itr = iter(chrs)
    comp_array = np.array([expand(sample) for sample in composition])
    chr = None
    c = 0

    with open(path, 'w') as output:
        output.write(makeLine(sample_list, 'CHR', 'POS'))
        for row in comp_array.T:

            if np.sum(row) <= 0:
                try:
                    # print('Terminating at {0}'.format(row))
                    chr = next(chr_itr)
                    # print('Moving to {0} after {1} positions'.format(chr, c))
                    c = 0
                    
                except:
                    pass
            else:
                # print(row)
                pos = c * section_length
                output.write(makeLine(row, chr, pos))
                c += 1

def writePrimaryClusters(chr_names, chr_breaks, clusters, path):

    prev = 0
    with open(path, 'w') as output:
        for name, break_point in zip(chr_names, chr_breaks):
            output.write('#{0}\n'.format(name))
            for i, cluster in enumerate(clusters[prev:break_point]):
                output.write('\n'.join([str(i)] + [' '.join(c) for c in cluster] + ['\n']))
            prev = break_point
