import os,sys
import argparse
import re
if __name__=='__main__':
	parser=argparse.ArgumentParser(description="Change the header of simulated fasta file for easy computation")

	inputFile = open("../simulated_reads_Test1.fa")
	refFile = open("../simulated_reads_Test1.cig")
	outFile = open("../simulated_reads_Test1.out.fa",'w')
	count = 0
	for r in inputFile:
		outLine = r
		if r.startswith('>'):
			match = re.match('^>(?P<seqName>seq\.\d+)(?P<type>[a|b]).+',r)
			hot = 1
			if match.group('type') in 'bB':
				hot = 2
			cigTokens = refFile.next().strip().split('\t')
			chrom = cigTokens[1]
			exons = cigTokens[4].split(', ')
			tempOutLine = ":".join(":".join([chrom,exon]) for exon in exons)
			outLine = "%s:%s/%d\n"%(match.group('seqName'),tempOutLine,hot)
		outFile.write(outLine)
	inputFile.close()
	refFile.close()
	outFile.close()