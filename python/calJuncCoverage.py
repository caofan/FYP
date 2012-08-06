#!/usr/bin/python

'''
This program takes a set of junctions and a set of reads. And it will
calculate the coverage of the anchor of the junctions (user specified region) using
the reads provided. And will output a wig format coverage file.
'''

import argparse
import os,sys
sys.path.append(".")
import readBed
import pysam
import csv
from Bio import SeqIO
'''
Read the set of junctions into the memory
'''

SIGNALS = set(["GTAG","GCAG","ATAC"])
def getJunctions(inputfile,printSignal=False):
    reader = readBed.BEDReader(inputfile)
    junctions = []
    leftRegions = {}
    rightRegions = {}
    for r in reader:
        tempJunc = r
        tempJunc['blockCount']= int(tempJunc['blockCount'])
        tempJunc['blockSizes']= tempJunc['blockSizes'].replace(',','\t').strip().split('\t')
        for index in range(len(tempJunc['blockSizes'])):
            tempJunc['blockSizes'][index] = int(tempJunc['blockSizes'][index])
        tempJunc['blockStarts'] = tempJunc['blockStarts'].replace(',','\t').strip().split('\t')
        for index in range(len(tempJunc['blockStarts'])):
            tempJunc['blockStarts'][index] = int(tempJunc['blockStarts'][index])
        tempJunc['chromStart'] = int(tempJunc['chromStart'])
        tempJunc['chromEnd'] = int(tempJunc['chromEnd'])
        tempJunc['score'] = int(tempJunc['score'])
        tempJunc['thickStart'] = int(tempJunc['thickStart'])
        tempJunc['thickEnd'] = int(tempJunc['thickEnd'])
        #get the signal
        junctions.append(tempJunc)
    junctions.sort(key=lambda k:(k['chrom'],k['chromStart'],k['chromEnd']))
    # get the signals of the junctions.
    fasta = SeqIO.read("/home/Ken/mm9/%s.fa"%(junctions[0]['chrom'],),"fasta")
    for junc in junctions:
        juncStart = junc['chromStart'] + junc['blockSizes'][0]
        juncEnd = junc['chromStart'] + junc['blockStarts'][1]
        leftSignal = fasta[juncStart:(juncStart+2)]
        rightSignal = fasta[(juncEnd-2):juncEnd]
        signal = leftSignal + rightSignal
        if signal.seq.reverse_complement().data.upper() in SIGNALS:
            junc['signal'] = signal.seq.reverse_complement().data.upper()
        else:
            junc['signal'] = signal.seq.data.upper()

    if printSignal:
        signals = {}
        for junc in junctions:
            juncStart = junc['chromStart'] + junc['blockSizes'][0]
            juncEnd = junc['chromStart'] + junc['blockStarts'][1]
            leftSignal = fasta[juncStart:(juncStart+2)]
            rightSignal = fasta[(juncEnd-2):juncEnd]
            signal = leftSignal + rightSignal
            if signal.seq.data.upper() in signals:
                signals[signal.seq.data.upper()] += 1
            elif signal.seq.reverse_complement().data.upper() in signals:
                signals[signal.seq.reverse_complement().data.upper()] += 1
            else:
                signals[signal.seq.data.upper()] = 1
        #for key in signals:
        #    print key," ", signals[key]
        return junctions,signals
    return junctions


def getReads(inputfile):
    reader = pysam.Samfile(inputfile,'rb')
    return reader


'''
Get the number of spanning reads for each junction.
The offset defines how far the read could be away from the junctions.
'''
def getSpanRead(junctions,readhandle,offset,readsOutFile):
    outSam = pysam.Samfile(readsOutFile,"wb",template=readhandle)
    for junc in junctions:
        # pysam uses 0-based coordinates. And the coordinates are inclusive.
        reads = readhandle.fetch(junc['chrom'],junc['chromStart']-offset,junc['chromEnd']+offset)
        juncstart = junc['chromStart'] + junc['blockSizes'][0]
        juncend = junc['chromStart'] + junc['blockStarts'][1]
        count = 0
        anchors = [0,0,]
        diff_large = 0
        diff_small = 1000
        for r in reads:
            # in order to take the uncertainty in Rawbin reads into consideration.
            currOffset = offset
            if r.tags:
                for tag in r.tags:
                    if "XO" in tag:
                        currOffset += int(tag[1]/2)
            if len(r.cigar) <= 1:
                continue
            i = 0
            tempstart = r.pos
            tempend = r.pos
            while i+2 < len(r.cigar):
                tempstart = tempend + r.cigar[i][1]
                tempend = tempstart + r.cigar[i+1][1]
                templeftAnchor = r.cigar[i][1]
                temprightAnchor = r.cigar[i+2][1]
                tempdiff = abs(templeftAnchor - temprightAnchor)
                if tempstart >= juncstart - currOffset and tempstart <= juncstart + currOffset and tempend >= juncend-currOffset and tempend <= juncend+currOffset:
                    count += 1
                    outSam.write(r)
                    if tempdiff > diff_large:
                        anchors[0] = min(templeftAnchor, temprightAnchor)
                        diff_large = tempdiff
                    if tempdiff < diff_small:
                        anchors[1] = min(templeftAnchor, temprightAnchor)
                        diff_small = tempdiff
                    break
                i += 2

        junc['spanReads'] = count
        junc['anchorSize'] = str(anchors[0]) + "," + str(anchors[1])
    outSam.close()
    os.system("samtools sort %s %s"%(readsOutFile,readsOutFile.split(".txt.bam")[0]+"_sorted"))
    os.system("samtools index %s"%(readsOutFile.split(".txt.bam")[0]+"_sorted.bam",))

'''
Obtain the coverage information for the junctions.
covrange defines the (juncstart - covrange, juncstart) and (juncend - 1, juncend + covrange - 1)
Defind in this way to include the boundary of the junctions.
'''
def getCovInfo(junctions,readhandle,covrange):
    for junc in junctions:
        juncstart = junc['chromStart'] + junc['blockSizes'][0]
        juncend = junc['chromStart'] + junc['blockStarts'][1]
        # The region includes the juncStart position.
        leftiter = readhandle.fetch(junc['chrom'],juncstart - covrange, juncstart)
        rightiter = readhandle.fetch(junc['chrom'],juncend-1,juncend+covrange-1)
        leftreads = []
        rightreads = []
        for r in leftiter:
            if len(r.cigar) <= 1:
                leftreads.append(r)
        for r in rightiter:
            if len(r.cigar) <= 1:
                rightreads.append(r)
        leftreads.sort(key=lambda k:(k.pos+k.qlen,))
        rightreads.sort(key=lambda k:(k.pos+k.qlen,))
        junc['anchorCov'] = str(len(leftreads))+','+str(len(rightreads))
        #0-base half-open indices

        leftcovered = 0
        rightcovered = 0
        for pos in range(juncstart-covrange,juncstart+1):
            for r in leftreads:
                if pos >= r.pos and pos < r.pos + r.qlen:
                    leftcovered += 1
                    break
        for pos in range(juncend-1,juncend+covrange):
            for r in rightreads:
                if pos >= r.pos and pos <= r.pos + r.qlen:
                    rightcovered += 1
                    break
        junc['baseCovered'] = str(leftcovered)+','+str(rightcovered)
        if len(leftreads) == 0:
            leftcovered = 0
            junc['disLeft'] = -100000
        else:
            junc['disLeft'] = leftreads[-1].pos + leftreads[-1].qlen - 1 - juncstart

        if len(rightreads) == 0:
            rightcovered = 0
            junc['disRight'] = -100000
        else:
            junc['disRight'] = rightreads[0].pos - juncend



def writeToFile(junctions, outputfile):
    FIELDS= ('chrom', 'chromStart', 'chromEnd',
              'name', 'score', 'strand',
              'thickStart', 'thickEnd',
              'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts','spanReads', 'anchorSize','anchorCov', 'baseCovered', 'disLeft', 'disRight')
    f = csv.DictWriter(open(outputfile,'wb'),fieldnames=FIELDS,delimiter='\t',extrasaction='ignore',quoting=csv.QUOTE_NONE)
    for j in junctions:
        blockSizes = ''
        for token in j['blockSizes']:
            blockSizes = blockSizes + str(token) + ','
        blockStarts = ''
        for token in j['blockStarts']:
            blockStarts = blockStarts + str(token) + ','
        blockSizes = blockSizes.replace(',','\t').strip().replace('\t',',')
        blockStarts = blockStarts.replace(',','\t').strip().replace('\t',',')
        j['blockSizes'] = blockSizes
        j['blockStarts'] = blockStarts
        f.writerow(j)

def printout(junctions):
    for j in junctions:
        blockSizes = ''
        for token in j['blockSizes']:
            blockSizes = blockSizes + str(token) + ','
        blockStarts = ''
        for token in j['blockStarts']:
            blockStarts = blockStarts + str(token) + ','
        j['blockSizes'] = blockSizes
        j['blockStarts'] = blockStarts
        print j['chrom'],'\t',j['chromStart'],'\t',j['chromEnd'],'\t',j['name'],'\t',j['score'],'\t',j['strand'],'\t',j['thickStart'],'\t',j['thickEnd'],'\t',j['itemRgb'],'\t',j['blockCount'],'\t',j['blockSizes'],'\t',j['blockStarts'],'\t',j['spanReads'],'\t',j['anchorSize'],'\t',j['anchorCov'],'\t',j['baseCovered'],'\t',j['disLeft'],'\t',j['disRight']


def getStats(junctions):
    noSupport = 0 #a
    lowCov = 0 #b
    nonCanonical = 0 #c, signal not canonical
    ab = 0
    bc = 0
    ca = 0
    abc = 0
    for junc in junctions:
        nslc = True
        nsnc = True
        lcnc = True
        tempAll = True
        if junc['spanReads'] == 0 or int(junc['anchorSize'].split(',')[1]) < 8:
            noSupport += 1
        else:
            nslc = False
            nsnc = False
            tempAll = False
        if junc['disLeft'] == -100000 or junc['disRight'] == -100000:
            lowCov += 1
        else:
            nslc = False
            lcnc = False
            tempAll = False
        if junc['signal'] not in SIGNALS:
            nonCanonical += 1
        else:
            nsnc = False
            lcnc = False
            tempAll = False
        if nslc:
            ab += 1
        if nsnc:
            ca += 1
        if lcnc:
            bc += 1
        if tempAll:
            abc += 1
    print 'NoSupport\tLowCoverage\tNotCanonical\tAB\tBC\tAC\tABC\tALL'
    print '%d\t%d\t%d\t%d\t%d\t%d\t%d\t%d'%(noSupport,lowCov,nonCanonical,ab,bc,ca,abc,len(junctions))

if __name__=="__main__":
    parser = argparse.ArgumentParser(description='''This program takes a set of junctions and a set of reads. And it will
calculate the coverage of the anchor of the junctions (user specified region) using
the reads provided. And will output a wig format coverage file.''')
    parser.add_argument('-j','--junctions',dest="junctions",required=True,help="The junctions to be used")
    parser.add_argument('-r','--reads',dest='reads',required=True,help="The reads file. It can be bam format and should be indexed.")
    parser.add_argument('-a','--offset',dest='offset',type=int,default=10,help="The size of the range allowed for a read to be a spanning read. Default is 10.")
    parser.add_argument('-c','--covrange',dest='cov',type=int,default=10,help="The range of anchors where the coverage should be computed.")
    parser.add_argument('-s','--stats',dest='stats',action="store_true",help="If set, stats will be generated")
    #parser.add_argument('--rbed',dest='rbed',action="store_true",help="If the read is in bed format, need to specify this.")
    parser.add_argument('-o','--output',dest='output',default="stdout",help="The output file")
    try:
        options= parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))
    if not os.path.exists(options.junctions):
        print "The junction file specified not exits."
        sys.exit(1)
    if not os.path.exists(options.reads):
        print "The reads file specified not exits."
        sys.exit(1)
    readsOutfile = "stdout"
    if options.output != "stdout":
        readsOutfile = options.output + ".bam"

    junctions = getJunctions(options.junctions)
    readhandle = getReads(options.reads)
    getSpanRead(junctions,readhandle,options.offset,readsOutfile)
    getCovInfo(junctions,readhandle,options.cov)
    if options.stats:
        getStats(junctions)

    if options.output != "stdout":
        writeToFile(junctions,options.output)
    else:
        printout(junctions)
