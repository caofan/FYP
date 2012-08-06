import os,sys
import argparse

def cig2bed(filename,juncOnly=False):
	"""docstring for cig2bed"""
	f = open(filename)
	out = open(filename.replace('cig','out.bed'),'w')
	count = 0
	for r in f:
		name,chrom,chromStart,cigar,pos,strand=r.split('\t')[0:6]
		chromStart = int(chromStart) - 1
		chromEnd = 0
		blockPoses = pos.split(', ')
		blockCount = len(blockPoses)
		if blockCount <= 1:
			continue
		blockStarts = []
		blockSizes = []
		for i in range(blockCount):
			start,end = [int(p) for p in blockPoses[i].split('-')]
			start -= 1
			blockSizes.append(end-start)
			blockStarts.append(start - chromStart)
			if i == blockCount - 1:
				chromEnd = end
		blockSizeStr = ",".join(str(i) for i in blockSizes)
		blockStartStr = ",".join(str(i) for i in blockStarts)
		out.write("%s\t%d\t%d\t%s\t%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n"
				%(chrom,chromStart,chromEnd,name,0,strand,chromStart,chromEnd,"255,0,0,0",blockCount,blockSizeStr,blockStartStr))
		count += 1
		if count%20000 == 0:
			print count, " records written"
	print "\nTotal ",count, " records written."
	f.close()
	out.close()


if __name__=="__main__":
	#if len(sys.argv) < 2:
	#	print "Please supply the input file name"
	#	sys.exit(1)
	#	cig2bed(filename)
	parser = argparse.ArgumentParser(description="Covert cig files to bed format.")
	parser.add_argument("inputs",nargs="+",help="The cig files to be converted")
	parser.add_argument("-j", "--juncOnly",dest="juncOnly",action="store_true",default=False,help="Specify whether to output only records with junctions.")
	#parser.add_argument("-o","--outDir",dest="outDir",help="The output directory.")
	try:
		options = parser.parse_args()
	except IOError, msg:
		parser.error(str(msg))
	for filename in options.inputs:
		print "="*12,filename,"="*12
		cig2bed(filename,options.juncOnly)
		print "\n\n"
