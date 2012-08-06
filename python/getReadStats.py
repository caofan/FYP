import sys
import re
import argparse
from collections import defaultdict

def writeRead(out,tokens):
    out.write(">%s\n%s\n"%(':'.join(tokens[0:-1]),tokens[-1]))

def processCig(filename,juncs=None):
    infile = open(filename)
    juncsInReads = defaultdict(int)
    minExon = defaultdict(int)
    maxIntron = 0
    totalReads = 0
    if juncs!=None:
        juncs = list(set(juncs))
        juncs.sort()
        outs = []
        for i in juncs:
            outs.append(open('%s.%d.fa'%(filename,i),'w'))
    for r in infile:
        totalReads += 1
        tokens = r.strip().replace(' ','').split('\t')
        cigar = tokens[3]
        #Ignore reads with insertions or deletions first.
        if 'I' in cigar or 'D' in cigar:
            continue
        juncCount = cigar.count('N')
        useRead = True
        
        #This part is for test purpose, assuming perfect match.
        iteratorN = re.finditer("(?P<n>\d+)N",cigar,flags=re.IGNORECASE)
        for n in iteratorN:
            if n.group('n'):
                intron = int(n.group('n'))
                if intron > maxIntron:
                    maxIntron = intron
                if intron < 25:
                    useRead = False
                    break
        
        if useRead and juncCount > 0:
            minimum=1000
            #iterator = re.finditer("((?P<m>\d+)M(\d+I)*\d+N)|((?P<n>\d+)M(\d+I)*$)",cigar,flags=re.IGNORECASE)
            iterator = re.finditer("(?P<m>\d+)M",cigar,flags=re.IGNORECASE)
            for m in iterator:
                if m.group('m') and int(m.group('m')) < minimum:
                    minimum = int(m.group('m'))
                #if m.group('n') and int(m.group('n')) < min:
                #    minimum = int(m.group('n'))
            juncsInReads[juncCount] += 1
            #if minimum >= 40:
            #    print r
            minExon[minimum] += 1
            if minimum < 19:
                useRead = False
        if juncs != None and useRead and juncCount in juncs:
            writeRead(outs[juncs.index(juncCount)],tokens)
    if juncs != None:
        for o in outs:
            o.close()
    
    print "Number of Juncs in Reads:"
    for key in juncsInReads:
        print "%d\t%d\t%f"%(key,juncsInReads[key],juncsInReads[key]*1.0/totalReads)

    print "\n"

    print "Minimum partition size:"
    for key in minExon:
        print "%d\t%d\t%f"%(key,minExon[key],minExon[key]*1.0/totalReads)

    print "Max Intron Size:"
    print maxIntron


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Calculate the distribution of size of the minimum exonic region in the read and the number of junctions within a read.")
    parser.add_argument('cigfiles',nargs='+',help="The cig files to be processed.")
    parser.add_argument('-j','--junctions',dest="junctions",nargs='*',type=int,help="Output reads with the number of junctions specified in separate files.")
    args = parser.parse_args()
    
    for cigfile in args.cigfiles:
        print cigfile
        processCig(cigfile,args.junctions)
