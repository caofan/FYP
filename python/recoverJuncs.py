import os,sys
import pysam
from Bio import SeqIO
import re
import readBed
import calJuncCoverage
readsDir = "../supposedVSfound/calJuncCovOut"

outDir = "../supposedVSfound/calJuncCovOut/recovered"

fastaDir = "../mm9-genome/fasta"

juncDir = '../supposedVSfound/'

if not os.path.exists(outDir):
    os.makedirs(outDir)

SIGNALS = set(["GTAG","GCAG","ATAC"])
count = 0
unique = 0

validated = 0

def validate(signals,chrom):
    juncfile = os.path.join(juncDir,'supposed.rawbin.splicemap.%s.new.bed'%(chrom,))
    out = readBed.BEDWriter(os.path.join(outDir,chrom+'_recovered.bed'))
    junctions = calJuncCoverage.getJunctions(juncfile)
    global validated
    for junc in junctions:
        leftindex = junc['chromStart'] + junc['blockSizes'][0] - 1
        rightindex = junc['chromStart'] + junc['blockStarts'][1]
        key = '_'.join([junc['chrom'],str(leftindex),str(rightindex)])
        if key in signals:
            out.writerow(junc)
            validated += 1




if __name__=="__main__":
    for readfile in os.walk(readsDir).next()[2]:
        match = re.match(r"^rawbin_(?P<chrom>chr(\d+|X|Y))_sorted.bam$",readfile)
        if match:# and readfile == "rawbin_chr2_sorted.bam":
            print '\n\n',readfile
            uniqueSignals = set()
            readhandle = pysam.Samfile(os.path.join(readsDir,readfile),"rb")
            chrom = match.group("chrom")
            fasta = SeqIO.read(os.path.join(fastaDir,chrom+".fa"),"fasta")
            try:
                read = readhandle.next()
            except StopIteration:
                continue
            offset = 0
            if read and read.tags:
                for tag in read.tags:
                    if tag[0] == "XO":
                        offset = tag[1]

            leftstart = read.pos
            #rightstart = read.pos + read.cigar[0][1] + read.cigar[1][1]
            #leftend = read.pos + read.cigar[0][1]
            # The region is 0-based, half open.
            juncStartRegion = [read.pos + read.cigar[0][1] - offset, read.pos + read.cigar[0][1]]
            juncEndRegion = [read.pos + read.cigar[0][1] + read.cigar[1][1], read.pos+read.cigar[0][1]+read.cigar[1][1] + offset]
            rightend = read.pos + read.cigar[0][1] + read.cigar[1][1] + read.cigar[2][1] + offset
            name = read.qname
            #skipregion = [read.pos + read.cigar[0][1], read.pos + read.cigar[0][1] + read.cigar[1][1],]
            existSignals = {}
            for read in readhandle:
                if read.tags:
                    for tag in read.tags:
                        if tag[0] == "XO":
                            offset = tag[1]
                #print read
                #tempSkipregion = [read.pos + read.cigar[0][1], read.pos + read.cigar[0][1] + read.cigar[1][1],]
                #skipSize = min(skipregion[1],tempSkipregion[1]) - max(skipregion[0],tempSkipregion[0])
                #if skipregion[1] - skipregion[0] == 0 or tempSkipregion[1] - tempSkipregion[0] == 0:
                #    print read
                weight = 1
                tempname = read.qname
                if tempname == name:
                    weight = 0
                else:
                    name = tempname

                templeftstart = read.pos
                temprightend = read.pos + read.cigar[0][1] + read.cigar[1][1] + read.cigar[2][1] + offset
                #if skipSize*1.0/(skipregion[1]-skipregion[0]) >= 0.9 and skipSize*1.0/(tempSkipregion[1]-tempSkipregion[0]) >=0.9:
                overlapSize = min(temprightend,rightend) - max(templeftstart,leftstart)
                if overlapSize*1.0 / (temprightend - templeftstart) >= 0.9 and overlapSize*1.0/(rightend-leftstart) >= 0.9:
                    #temprightstart = read.pos + read.cigar[0][1] + read.cigar[1][1]
                    #templeftend = templeftstart + read.cigar[0][1]
                    tempjuncStartRegion =  [read.pos + read.cigar[0][1] - offset, read.pos + read.cigar[0][1]]
                    tempjuncEndRegion = [read.pos + read.cigar[0][1] + read.cigar[1][1], read.pos+read.cigar[0][1]+read.cigar[1][1] + offset]
                    if templeftstart < leftstart:
                        leftstart = templeftstart
                    #if temprightstart < rightstart:
                    #    rightstart = temprightstart
                    #if templeftend > leftend:
                    #    leftend = templeftend
                    if temprightend > rightend:
                        rightend = temprightend

                    templeftindex = tempjuncStartRegion[1] - 1
                    while templeftindex >= tempjuncStartRegion[0]:
                        tempLeftSignal = fasta[templeftindex+1:templeftindex+3]
                        tempLeftSize = templeftindex - templeftstart + 1
                        #temprightindex = tempjuncEndRegion[0]
                        temprightindex = temprightend + tempLeftSize - 76
                        while temprightend - temprightindex >= 76-tempLeftSize and temprightindex < tempjuncEndRegion[1]:
                            tempRightSignal = fasta[temprightindex-2:temprightindex]
                            tempsignal = tempLeftSignal + tempRightSignal
                            if tempsignal.seq.data in SIGNALS:
                                tempkey = tempsignal.seq.data+'_'+chrom+"_"+str(templeftindex)+'_'+str(temprightindex)
                                existSignals[tempkey] = existSignals.get(tempkey,0) + weight
                            elif tempsignal.seq.reverse_complement().data in SIGNALS:
                                tempkey = tempsignal.seq.reverse_complement().data+"_"+chrom+'_'+str(templeftindex)+'_'+str(temprightindex)
                                existSignals[tempkey] = existSignals.get(tempkey,0) + weight
                            temprightindex += 1
                        templeftindex -= 1
                    #if tempjuncStartRegion[0] < juncStartRegion[0]:
                    #    juncStartRegion[0] = tempjuncStartRegion[0]
                    #if tempjuncStartRegion[1] > juncStartRegion[1]:
                    #    juncStartRegion[1] = tempjuncStartRegion[1]
                    #if tempjuncEndRegion[0] < juncEndRegion[0]:
                    #    juncEndRegion[0] = tempjuncEndRegion[0]
                    #if tempjuncEndRegion[1] > juncEndRegion[1]:
                    #    juncEndRegion[1] = tempjuncEndRegion[1]
                    #skipregion = [min(skipregion[0],tempSkipregion[0]),max(skipregion[1],tempSkipregion[1]),]
                else:
                    #print name
                    #existSignals = {}
                    print '------------------------------------'
                    print leftstart,' ',rightend
                    for key in existSignals:
                        print key,' ', existSignals[key]

                    if len(existSignals.keys()) > 0:
                        count += 1
                    if len(existSignals.keys()) == 1:
                        unique += 1
                        for key in existSignals:
                            tempkey = '_'.join(key.split('_')[1:])
                            uniqueSignals.add(tempkey)
                    existSignals = {}
                    leftstart = read.pos
                    #rightstart = read.pos + read.cigar[0][1] + read.cigar[1][1]
                    #leftend = read.pos + read.cigar[0][1]
                    juncStartRegion = [read.pos + read.cigar[0][1] - offset, read.pos + read.cigar[0][1]]
                    juncEndRegion = [read.pos + read.cigar[0][1] + read.cigar[1][1], read.pos+read.cigar[0][1]+read.cigar[1][1] + offset]
                    rightend = read.pos + read.cigar[0][1] + read.cigar[1][1] + read.cigar[2][1] + offset
                    #skipregion = [read.pos + read.cigar[0][1], read.pos + read.cigar[0][1] + read.cigar[1][1],]
                    #name = tempname
                    leftindex = juncStartRegion[1]-1
                    #print '--------------------------------------'
                    #print leftstart,' ', ' ',' ',rightend
                    #print juncStartRegion
                    #print juncEndRegion
                    while leftindex >= juncStartRegion[0]:
                        leftSignal = fasta[leftindex+1:leftindex+3]
                        leftSize = leftindex - leftstart + 1
                        #rightindex = juncEndRegion[0]
                        rightindex = rightend + leftSize - 76
                        while rightend - rightindex >= 76-leftSize and rightindex < juncEndRegion[1]:
                            rightSignal = fasta[rightindex-2:rightindex]
                            signal = leftSignal + rightSignal
                            if signal.seq.data in SIGNALS:
                                tempkey = signal.seq.data + "_"+chrom+'_'+str(leftindex)+'_'+str(rightindex)
                                existSignals.setdefault(tempkey,weight)
                            elif signal.seq.reverse_complement().data in SIGNALS:
                                tempkey = signal.seq.reverse_complement().data + '_'+chrom+'_' + str(leftindex) + '_' + str(rightindex)
                                existSignals.setdefault(tempkey, weight)
                            rightindex += 1
                        leftindex -= 1
            print '------------------------------------'
            print leftstart,' ',rightend
            for key in existSignals:
                print key,' ', existSignals[key]

            if len(existSignals.keys()) > 0:
                count += 1
            if len(existSignals.keys()) == 1:
                unique += 1
                for key in existSignals:
                    tempkey = '_'.join(key.split('_')[1:])
                    uniqueSignals.add(tempkey)

            validate(uniqueSignals,chrom)
    print "---------------------Summary-----------------------------"
    print "total ", count
    print "unique ", unique
    print "validated ", validated
