#!/usr/bin/python
# To run this program, Biopython and python 2.7 is required.
#This program is designed to work with the refGene files for annotation and fasta files for chromosome sequence
import io,sys
import argparse
import Bio
from Bio import SeqIO
import random
parser = argparse.ArgumentParser(description="Simulate the reads by uniform distribution.")
parser.add_argument('-r','--refGene',dest='refGeneSource', help='The path to refGene file.',action='store',required=True)
parser.add_argument('-g','--genome',dest='chromosomes', help='The directory where the chromosome fasta files are stored.',required=True)
parser.add_argument('-m', '--mode',dest='mode',help='a: adjacent. u: uniform. c: need not be adjacent. n: nonuniform.',choices=['au','an','cu','cn'], default='au', action='store')
parser.add_argument('-l','--read',dest='readLength', help='the length of the reads to generate, default is 76', type=int, default=76)
parser.add_argument('-c','--coverage',dest='coverage',help='The coverage of the reads, default is 50. The coverage is calculated as readLength*numOfReads/lengthOfSequence. This option is silenced if the mode is not uniform.', type=int, default=50)
parser.add_argument('-o','--output',dest='output', help='the location of the output file',required=True)
parser.add_argument('-p','--pair',dest='pair',action='store_true',help='This option indicates to generate pair-end reads.')
parser.add_argument('--minspan',dest='minSpan',help='This option indicates the minimum length that the two tags of a paired-end read should be separated.',type=int, default=10)
parser.add_argument('--maxspan',dest='maxSpan',help='This option indicates the maximum length that the two tags of a paired-end read could be separated.', type=int,default=300)

def deleteRandom(refGene):
    i=0
    while i != len(refGene)-1:
        if 'random' in refGene[i][2]:
            refGene.remove(refGene[i])
        else:
            i+=1

def getRefGene(path):
    f=open(path)
    tempread = f.readlines()
    f.close()
    refGene=[]
    for line in tempread:
        refGene.append(line.split('\t'))
    deleteRandom(refGene)
    refGene.sort(key=lambda gene:(gene[2]))
    return refGene

def generateReads(refGene):
    currentChr = refGene[0][2]
    sequence=SeqIO.read('%s/%s.fa'%(options.chromosomes,currentChr),'fasta')
    print "Generating reads for chromosome " + sequence.description
    f=open(options.output,'w')
    countChr = 0
    for gene in refGene:
        if gene[2] != currentChr:
            print '%d reads generated for chromosome %s'%(countChr,currentChr)
            currentChr = gene[2]
            sequence=SeqIO.read('%s/%s.fa'%(options.chromosomes,currentChr),'fasta')
            countChr=0
            print "Generating reads for chromosome " + sequence.description
        strand = gene[3]
        numExons = int(gene[8])
        exonStarts = gene[9].split(',')
        exonEnds = gene[10].split(',')
        exons = []
        for i in range(numExons):
            exons.append([int(exonStarts[i]),int(exonEnds[i])])
        exons.sort(key=lambda e:(e[0],e[1]))
        if options.mode=='au':
            countChr += adjUni(exons,strand,sequence,f)
        elif options.mode=='an':
            pass
        elif options.mode=='cu':
            pass
        elif options.mode=='cn':
            pass
    f.close()

class Read:
    def __init__(self,sn,s1,e1,sp1,s2,e2,sp2,st,m,seq):
        self.sequence_name = sn
        self.left1 = s1
        self.right1 = e1
        self.span1 = sp1
        self.left2 = s2
        self.right2 = e2
        self.span2 = sp2
        self.strand = st
        self.motif = m
        self.sequence = seq

    def __str__(self):
        return '>%s:%d:%d:%d:%d:%d:%d:%s:%s\n%s'%(self.sequence_name,self.left1,self.right1,self.span1,self.left2,self.right2,self.span2,self.strand,self.motif,self.sequence.tostring())

def adjUni(exons,strand,sequence,f):
    totalLength = 0
    #exonStarts is used to store the relative position of the exons when all exons are concatenated.
    exonStarts=list()
    for i,exon in enumerate(exons):
        exonStarts.append(totalLength)
        totalLength += exon[1]-exon[0]+1
    numOfReads = options.coverage*totalLength/options.readLength + 1
    reads=[]
    exonStarts.reverse()
    for i in range(numOfReads):
        tempIndex = random.randint(0,totalLength-1)
        tempEnd = tempIndex+options.readLength-1
        if tempEnd >= totalLength:
            i -= 1
            continue
        for j,exonstart in enumerate(exonStarts):
            #If the two ends are on different exons.
            if tempIndex >= exonstart and j>0 and tempEnd>=exonStarts[j-1]:
                tempdis1 = tempIndex-exonstart
                k=0
                for k in range(j):
                    if tempEnd>=exonStarts[k]:
                        break
                tempdis2 = tempEnd-exonStarts[k]
                exon1 = exons[len(exons)-1-j]
                exon2 = exons[len(exons)-1-k]
                if exon1[1]-exon1[0]-tempdis1+1>1 and tempdis2+1>1:
                    motif='%s%s%s%s'%(sequence[exon1[1]-1],sequence[exon1[1]],sequence[exon2[0]],sequence[exon2[0]+1])
                else:
                    motif=''
                readseq = (sequence[exon1[0]+tempdis1:exon1[1]+1]+sequence[exon2[0]:exon2[0]+tempdis2+1]).seq
                #note that the motif is for both the reverse and forward strand are based on the motif on the corresponding forward strand.
                if strand=='-':
                    readseq=readseq.reverse_complement()
                # All coordinates are based on the forward strand and left to right. Thus for those on the reverse strand, the numbers should be read
                # in the reverse direction.
                reads.append(Read(sequence.name,exon1[0]+tempdis1,exon1[1],exon1[1]-exon1[0]-tempdis1+1,exon2[0],exon2[0]+tempdis2,tempdis2+1,strand,motif,readseq))
                break
            elif tempIndex >= exonstart:
                tempdis1 = tempIndex - exonstart
                tempdis2 = tempEnd - exonstart
                currExon = exons[len(exons)-1-j]
                motif=''
                readseq = sequence[currExon[0]+tempdis1:currExon[0]+tempdis2+1].seq
                if strand=='-':
                    readseq=readseq.reverse_complement()
                reads.append(Read(sequence.name,currExon[0]+tempdis1,currExon[0]+tempdis2,options.readLength,0,0,0,strand,motif,readseq))
                break
    reads.sort(key=lambda r:(r.motif,r.left1,r.left2))
    for r in reads:
        f.write(str(r)+'\n')
    return numOfReads

if __name__ == "__main__":
    try:
        global options
        options=parser.parse_args()
        #f=open('test.txt','w')
        refGene = getRefGene(options.refGeneSource)
        generateReads(refGene)
    except IOError, msg:
        parser.error(str(msg))
