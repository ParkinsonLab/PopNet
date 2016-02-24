'''
This function outputs the color table to a file, to assist with other scripts that
would like this kind of information.
'''


def outputColors(groups, colorTable, outpath):
    import CytoscapeEncoder as ce
    with open(outpath, 'w') as output:
        for group in groups:
            output.write('{name}\t{value}\n'.format(name=group, value=ce.toHexColor(colorTable[group])))
    print('Colors output completed.')