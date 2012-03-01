#!/usr/python
import os,sys
from Bio import SeqIO
import readBed
import argparse

def 

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Get the splicing signals given sets of junction in bed format. The output file will be inputName+'.modified.bed'")
	parser.add_argument('genome',nargs=1,help='genome fasta file')
	parser.add_argument('bedFiles',nargs='+',help='bed files containing junctions')

	args = parser.parse_args()
	genomeGen = SeqIO.parse(args.genome,'fasta')
	junctions = []
	for bedFile in args.bedFiles:
		tempJuncs = {}
		reader = readBed.BEDReader(bedFile)
		for r in reader:
			if r['chrom'] in tempJuncs:
				tempJuncs[r['chrom']].append(r)
			else:
				tempJuncs[r['chrom']] = [r,]
		junctions.append(tempJuncs)
	chromNames = []
	for chrom in genomeGen:
		chromNames.append(chrom.name)
		for juncSet in junctions:
			for junc in juncSet[chrom.name]:
				blockSizes = [int(i) for i in junc['blockSizes'].replace(',','\t').strip().split('\t')]
				blockStarts = [int(i) for i in junc['blockStarts'].replace(',','\t').strip().split('\t')]
				chromStart = int(junc['chromStart'])
				chromEnd = int(junc['chromEnd'])
				juncStart = chromStart + blockSizes[0]
				juncEnd = chromStart + blockStarts[1]
				leftSig = str(chrom[(juncStart+1):(juncStart+3)].seq)
				rightSig = str(chrom[(juncEnd-2):(juncEnd)].seq)
				sig = leftSig+'_'+rightSig
				junc['name'] = ":".join([junc['name'],sig])

	for bedFile, juncSet in zip(args.bedFiles,junctions):
		out = readBed.BEDWriter(bedFile+'.modified.bed')
		for chromName in chromNames:
			for junc in juncSet[chromName]:
				out.write(junc)
