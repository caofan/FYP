#!/usr/bin/python
'''
Extract reads that spans the splicing junctions.
'''
import os,sys
import readBed
import pysam
import argparse
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqIO.FastaIO import FastaWriter

READLENGTH = 75

def overlap(s1,e1,s2,e2):
    #print s1,' ',e1,' ',s2,' ',e2
    if (s1 >= s2 and s1 <= e2) or (e1 >= s2 and e1 <= e2) or (s2 >= s1 and s2 <= e1) or (e2>=s1 and e2 <= e1):
        return True
    return False

def getSplitReads(juncs,bamhandle):
    splitreads = []
    juncsNoneHit = []
    multisplitreads = []
    for j in juncs:
        blockSizes = []
        for s in j['blockSizes'].replace(',','\t').strip().split('\t'):
            blockSizes.append(int(s))
        blockStarts = []
        for s in j['blockStarts'].replace(',','\t').strip().split('\t'):
            blockStarts.append(int(s))


        #all coordinates are 0-based.
        juncStart = int(j['chromStart']) + blockSizes[0] - 1
        juncEnd = int(j['chromStart']) + blockStarts[1]


        count = 0
        tempreads = bamhandle.fetch(j['chrom'],juncStart-10,juncEnd+10)
        for r in tempreads:
            if len(r.cigar) <= 1:
                continue
            readStart = r.pos+r.cigar[0][1] - 1
            readEnd = r.pos + r.cigar[0][1] + r.cigar[1][1]
            if readStart == juncStart and readEnd == juncEnd:
                count += 1
                splitreads.append(r)
            elif len(r.cigar) > 3:
                i = 2
                #print r
                while i+1 < len(r.cigar):
                    #print readStart, ' ',readEnd
                    readStart = readEnd + r.cigar[i][1]-1
                    readEnd = readEnd + r.cigar[i][1] + r.cigar[i+1][1]
                    if readStart == juncStart and readEnd == juncEnd:
                        count += 1
                        multisplitreads.append(r)
                        break
                    i+=2
                #print readStart,' ',readEnd




        # 0-based, open end [,)
        #tempreads1 = bamhandle.fetch(j['chrom'],juncStart-2,juncStart+2)
        #tempreads2 = bamhandle.fetch(j['chrom'],juncEnd-2,juncEnd+2)
        #reads1 = []
        #reads2 = []
        #for i in tempreads1:
        #    reads1.append(i)
        #for i in tempreads2:
        #    reads2.append(i)
        #spanreads = []
        #print j
        #for r1 in reads1:
            #print r1
        #    for r2 in reads2:
        #        if r1.compare(r2) == 0:
        #            spanreads.append(r1)
        #            #print r1
        #            break
        #count = 0
        # Sam files use 1-based coordinates and needs to be transformed to be 0-based.
        #for r in spanreads:
        #    if len(r.cigar) <= 1:
        #        continue
        #    readStart1 = r.pos-1
        #    readEnd1 = readStart1 + r.cigar[0][1] - 1
        #    readStart2 = r.pos-1 + r.cigar[0][1]+r.cigar[1][1]
        #    readEnd2 = readStart2 + r.cigar[2][1]-1
            # to get te actual position of the end instead of obtaining the next position of the end
            # Require the start and end positions of the read be located within the two anchor regions.
        #    if overlap(readStart1,readEnd1,juncStart-2,juncStart+2) and overlap(readStart2,readEnd2,juncEnd-2,juncEnd+2):
                #print r
        #        splitreads.append(r)
        #        count+=1
        if count == 0:
            #print j
            #print '\n'
            juncsNoneHit.append(j)

    #print len(splitreads)
    return splitreads,multisplitreads,juncsNoneHit


if __name__=='__main__':
    parser = argparse.ArgumentParser(description="Extract reads that spans the splicing junction in bed format.\nIf the -m, -b, -s options are not speficified, the corresponding files will not be generated.")
    parser.add_argument('juncFile',nargs=1,help='The file contains the junctions')
    parser.add_argument('hitsFile',nargs=1,help='The file in bam format that contains the hits information')
    parser.add_argument('-o','--outdir',dest='outdir',required=True,help='The output directory')
    parser.add_argument('-m','--multidir',dest='multireads', help = 'The output directory for multi split reads')
    parser.add_argument('-b','--bamdir',dest='bamdir',help='The output directory for the bam files')
    #parser.add_argument('-s','--samdir',dest='samdir',help='The output direcotry for the sam files')
    parser.add_argument('-n','--nonehit',dest='nonehit',help='The output direcotry for the bed files containing the juncs that has no hit')
    try:
        options = parser.parse_args()
    except IOError,msg:
        parser.error(str(msg))
    juncs = []
    reader = readBed.BEDReader(options.juncFile[0])
    for r in reader:
        juncs.append(r)
    #print options.hitsFile
    bamhandle = pysam.Samfile(options.hitsFile[0],'rb')
    splitreads,multisplitreads,juncsNoneHit = getSplitReads(juncs,bamhandle)


    #Output all supporting reads including the multi-split reads
    seqs = []
    for r in splitreads:
        tempSeq = SeqRecord(Seq(r.seq),id=r.qname,description='')
        seqs.append(tempSeq)
        #print r
    for r in multisplitreads:
        tempSeq = SeqRecord(Seq(r.seq),id=r.qname,description='')
        seqs.append(tempSeq)
    #print len(seqs)
    handle = open(os.path.join(options.outdir,os.path.basename(options.juncFile[0])+'.fa'),'w')
    FastaWriter(handle,wrap=80).write_file(seqs)
    #SeqIO.write(seqs,os.path.join(options.outdir,os.path.basename(options.juncFile[0])+'.fa'),'fasta')
    handle.close()



    if options.bamdir:
    #output the split reads in the junctions.
        output = pysam.Samfile(os.path.join(options.bamdir,os.path.basename(options.juncFile[0])+'.splitreads.bam'),'wb',template=bamhandle)
        for r in splitreads:
            output.write(r)
        output.close()

    #multisplit sequences
    if options.multireads:
        multiSeqs = []
        for r in multisplitreads:
            tempSeq = SeqRecord(Seq(r.seq),id=r.qname,description='')
            multiSeqs.append(tempSeq)
        handle = open(os.path.join(options.multireads,os.path.basename(options.juncFile[0])+'.multisplitreads.fa'),'w')
        FastaWriter(handle,wrap=80).write_file(multiSeqs)
        handle.close()
        output = pysam.Samfile(os.path.join(options.bamdir,os.path.basename(options.juncFile[0])+'.multisplitreads.bam'),'wb',template=bamhandle)
        for r in multisplitreads:
            output.write(r)
        output.close()


    # The nonehit files
    if options.nonehit:
        output = readBed.BEDWriter(os.path.join(options.nonehit,os.path.basename(options.juncFile[0]+'.nonehit.bed')))
        output.writerows(juncsNoneHit)
