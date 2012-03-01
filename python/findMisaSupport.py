import os,sys
sys.path.append('.')
import re
import readBed
import argparse
import math
def overlap(start1,end1,start2,end2,f):
    olength = 0
    if start1 >= start2 and start1 < end2:
        if end1 > end2:
            olength = end2 - start1
        else:
            olength = end1 - start1
    elif start2 >= start1 and start2 < end1:
        if end1 > end2:
            olength = end2 - start2
        else:
            olength = end1 - start2
    minLen = max(math.ceil((end1-start1)*f),math.ceil((end2-start2)*f))
    if olength > minLen:
        return True
    else:
        return False

if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Find the supporting junctions")
    parser.add_argument('-r','--reference',dest='ref',required=True)
    parser.add_argument('-i','--input',dest='infile',required=True)
    try:
        options = parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))
    refReader = readBed.BEDReader(options.ref)
    inReader = readBed.BEDReader(options.infile)
    refJuncs = []
    inJuncs = []
    output = []
    minF = 0.9
    count = 0
    for r in refReader:
        refJuncs.append(r)
    for i in inReader:
        output.append([i,])
        for j in refJuncs:
            if overlap(int(i['chromStart']),int(i['chromEnd']),int(j['chromStart']),int(j['chromEnd']),minF):
                output[count].append(j)
        count+=1
    for i in output:
        for j in i:
            print j['chrom'],'\t',j['chromStart'],'\t',j['chromEnd'],'\t',j['name']
        print '\n'

