import os,sys
import argparse
import re
import readBed

if __name__=='__main__':
    wd = '../outputs/cufflink.400000/rawbin.4.noann/'
    infile = os.path.join(wd,'junctions.bed')
    outdir = os.path.join(wd,'cat/')

    rsup = []
    csup = []
    reader = readBed.BEDReader(infile)
    count = 0
    for r in reader:
        name = r['name']
        #print name
        m=re.match(r'^JUNC\[(?P<read>\d+)\]_\d+_(?P<signal>[ACGT]{2}-[ACGT]{2}).+',name)
        count+=1
        if m:
            if int(m.group('read')) == 0:
                csup.append(r)
            else:
                rsup.append(r)
    print len(csup)
    print len(rsup)
    print len(csup)+len(rsup)
    print count
    writer = readBed.BEDWriter(os.path.join(outdir,'rsup.bed'))
    writer.writerows(rsup)
    writer = readBed.BEDWriter(os.path.join(outdir,'csup.bed'))
    writer.writerows(csup)

