#!/usr/bin/python
'''
This program is used to classify the junctions that are found by SpliceMap but not by Rawbin.
'''
import os,sys
sys.path.append('.')
import re
import readBed
import argparse
import pysam

parser = argparse.ArgumentParser(description="Options")
parser.add_argument('-i','--input',dest='inputfile',required=True,help="The path to the file")
parser.add_argument('-o','--output',dest='outputdir',required=True,help='The path for output dir')
parser.add_argument('-r','--reads',dest='readsfile',help='The file that contains the reads')

# classify the junctions based on whether they are supported by uniquely mapped reads or multiply mapped reads,
# or if they are supported by both kinds of junctions.
def uniquevsmulti(junctions):
    unique = []
    multi = []
    both = []
    multiout = []
    uniqueout = []
    bothout = []
    global options

    for j in junctions:
        if j[3] == 0:
            multi.append(j)
            multiout.append(j[-1])
        elif j[4] == 0:
            unique.append(j)
            uniqueout.append(j[-1])
        else:
            both.append(j)
            bothout.append(j[-1])

    f = readBed.BEDWriter(os.path.join(options.outputdir,'multi.bed'))
    f.writerows(multiout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,'unique.bed'))
    f.writerows(uniqueout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,'both.bed'))
    f.writerows(bothout)

    return unique,multi,both

# Differentiate which junctions are supported by multiple reads and which are supported by single read.
def singlevsmany(junctions,filename):
    single = []
    many = []
    singleout = []
    manyout = []
    if filename == 'multi':
        for j in junctions:
            if j[4] == 1:
                single.append(j)
                singleout.append(j[-1])
            else:
                many.append(j)
                manyout.append(j[-1])
    elif filename == 'unique':
        for j in junctions:
            if j[3] == 1:
                single.append(j)
                singleout.append(j[-1])
            else:
                many.append(j)
                manyout.append(j[-1])
    elif filename == 'both':
        for j in junctions:
            if j[3] == 1 and j[4] == 1:
                single.append(j)
                singleout.append(j[-1])
            else:
                many.append(j)
                manyout.append(j[-1])
    global options
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.single.bed'))
    f.writerows(singleout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.many.bed'))
    f.writerows(manyout)

    return single,many


#Check whether both anchors have coverage.
def anchorCov(junctions,filename):
    global options
    #look at the coverage of 10bp from the anchor
    range = 10
    bamfile = pysam.Samfile(options.readsfile,'rb')
    both = []
    single = []
    none = []
    bothout = []
    singleout = []
    noneout = []

    for j in junctions:
        blocksizes = [int(j[-1]['blockSizes'].split(',')[0]),int(j[-1]['blockSizes'].split(',')[1]),]
        blockstarts = [int(j[-1]['blockStarts'].split(',')[0]),int(j[-1]['blockStarts'].split(',')[1]),]
        start1 = int(j[-1]['chromStart']) + blocksizes[0] - range
        end1 = start1 + range
        start2 = int(j[-1]['chromStart']) + blockstarts[1]
        end2 = start2 + range
        chrom = j[-1]['chrom']

        reads1 = []
        reads2 = []
        tempreads = bamfile.fetch(chrom,start1,end1)
        for r in tempreads:
            reads1.append(r)
        tempreads = bamfile.fetch(chrom,start2,end2)
        for r in tempreads:
            reads2.append(r)

        countleft = 0
        countright = 0
        for r in reads1:
            if len(r.cigar) == 1 and r.pos + r.cigar[0][1] >= start1 and r.pos+r.cigar[0][1]<end1:
                countleft+=1
        for r in reads2:
            if len(r.cigar) == 1 and r.pos >= start2 and r.pos < end2:
                countright += 1
        if countleft == 0 and countright == 0:
            none.append(j)
            noneout.append(j[-1])
        elif countleft == 0 or countright == 0:
            single.append(j)
            singleout.append(j[-1])
        else:
            both.append(j)
            bothout.append(j[-1])


    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.1.bed'))
    f.writerows(singleout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.2.bed'))
    f.writerows(bothout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.0.bed'))
    f.writerows(noneout)

    return none,single,both




#Check if there are any junctions nearby.
def nearby(current,junctions,diff,filename):
    unique = []
    uniqueout = []
    nun = []
    nunout = []
    for j in current:
        start = int(j[-1]['chromStart'])
        end = int(j[-1]['chromEnd'])
        count = 0
        for row in junctions:
            tempstart = int(row[-1]['chromStart'])
            tempend = int(row[-1]['chromEnd'])
            if abs(start-tempstart)<diff and abs(end-tempend)<diff:
                count += 1
                if count > 1:
                    break
        if count > 1:
            nun.append(j)
            nunout.append(j[-1])
        else:
            unique.append(j)
            uniqueout.append(j[-1])

    global options
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.uni.bed'))
    f.writerows(uniqueout)
    f = readBed.BEDWriter(os.path.join(options.outputdir,filename+'.nun.bed'))
    f.writerows(nunout)

    return unique,nun


if __name__=="__main__":
    try:
        global options
        options = parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))
    reader = readBed.BEDReader(options.inputfile)
    if options.readsfile:
        bamfile = pysam.Samfile(options.readsfile,'rb')
    line = reader.next()
    while not line['name']:
        line = reader.next()
    junctions = []
    while line:
        match = re.match(r'\((?P<nR>\d+)\)\[(?P<width>\d+)_(?P<nNR>\d+)\]\((?P<nUR>\d+)\/(?P<nMR>\d+)\)',line['name'])
        junctions.append((int(match.group('nR')),int(match.group('width')),int(match.group('nNR')),int(match.group('nUR')),int(match.group('nMR')),line))
        try:
            line = reader.next()
        except StopIteration:
            break
    junctions.sort(key=lambda k:(k[3],k[4],k[0],k[1]))

    unique,multi,both = uniquevsmulti(junctions)

    single,many = singlevsmany(unique,'unique')
    zero,one,two = anchorCov(single,'unique.single')
    uni,nun = nearby(zero,junctions,10,'unique.single.0')
    uni,nun = nearby(one,junctions,10,'unique.single.1')
    uni,nun = nearby(two,junctions,10,'unique.single.2')

    zero,one,two = anchorCov(many,'unique.many')
    uni,nun = nearby(zero,junctions,10,'unique.many.0')
    uni,nun = nearby(one,junctions,10,'unique.many.1')
    uni,nun = nearby(two,junctions,10,'unique.many.2')


    single,many = singlevsmany(multi,'multi')
    zero,one,two = anchorCov(single,'multi.single')
    uni,nun = nearby(zero,junctions,10,'multi.single.0')
    uni,nun = nearby(one,junctions,10,'multi.single.1')
    uni,nun = nearby(two,junctions,10,'multi.single.2')

    zero,one,two = anchorCov(many,'multi.many')
    uni,nun = nearby(zero,junctions,10,'multi.many.0')
    uni,nun = nearby(one,junctions,10,'multi.many.1')
    uni,nun = nearby(two,junctions,10,'multi.many.2')


    single,many = singlevsmany(both,'both')
    zero,one,two = anchorCov(single,'both.single')

    uni,nun = nearby(zero,junctions,10,'both.single.0')
    uni,nun = nearby(one,junctions,10,'both.single.1')
    uni,nun = nearby(two,junctions,10,'both.single.2')

    zero,one,two = anchorCov(many,'both.many')
    uni,nun = nearby(zero,junctions,10,'both.many.0')
    uni,nun = nearby(one,junctions,10,'both.many.1')
    uni,nun = nearby(two,junctions,10,'both.many.2')



