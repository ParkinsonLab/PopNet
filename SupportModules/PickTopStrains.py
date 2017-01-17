'''
Created on Dec 18, 2014

@author: javi
'''

if __name__ == '__main__':
    
    import re
    import os
    from shutil import move
    
    masterpath = '/data/new/javi/plasmo'
    accpath = masterpath + '/MasterAccList.txt'
    llpath = masterpath + '/loclist.txt'
    fdirectory = '/data/new/javi/plasmo/ena'
    outdir = fdirectory + '/selected'
    outpath = masterpath + '/finalList.txt'
    
    #load the masterlist
    mdict = {}
    with open(accpath, 'r') as accfile:
        lines = re.split('\n', accfile.read())[1:-1]
    for line in lines:
        pair = re.split('\t', line)
        mdict[pair[1]] = pair[0]
    
    #load the loclist
    locdict = {}
    loclist= []
    with open(llpath, 'r') as llfile:
        sections = re.split('@', llfile.read())[1:]
    for section in sections:
        lines = re.split('\n', section)[:-1]
        loclist.append(lines[0])
        cleanedlines = [re.split('\s', x) for x in lines[1:]]
        for line in cleanedlines:
            locdict[line[0]] = " ".join(line[1:])
    
    #build approved list
    files = []
    os.chdir(fdirectory)
    for file in os.listdir(fdirectory):
        try:
            size = os.path.getsize(file)
            acc = re.match("([\S].+?)(?:_)", file).group(1)
            loc = locdict[mdict[acc]]
            files.append((file, size, acc, loc))
        except:
            continue
    files = sorted(files, key=lambda x:x[1], reverse=True)
    
    approvecount = {}
    approvedlist = []
    for loc in loclist:
        approvecount[loc] = 0
        
    for file in files:
        loc = file[3]
        acc = file[2]
        if acc not in approvedlist and approvecount[loc] < 20:
            approvedlist.append(file[2])
            approvecount[loc] += 1
    
    #move files
    if not os.path.isdir(outdir):
        os.mkdir(outdir)
        
    for file in files:
        filename = file[0]
        acc = file[2]
        if acc in approvedlist:
            move(filename, outdir)
            
    #create the log
    files = sorted(files, key=lambda x:locdict[mdict[x[2]]])
    with open(outpath, 'w') as outfile:
        for file in files:
            acc = file[2]
            if acc in approvedlist:
                outfile.write("{0} {1} {2} {3}\n".format(acc, mdict[acc], locdict[mdict[acc]], file[1]))
            
    print('script end')
    