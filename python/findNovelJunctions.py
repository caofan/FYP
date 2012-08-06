#!/usr/bin/python
import os,sys
sys.path.append('.')
import readBed
import readFlat
import platform
import re



'''
Generage dictionary of the junctions in a dict structure. The keys
are "chr_start_end" strings.
'''
def generateRefDict(reference):
	refDict = set()
	count = 0
	juncCount = 0
	chromPat = re.compile(r'^chr((\d+)|X|Y)$')
	if reference.endswith('txt'):
		ref = readFlat.FlatReader(reference)
		print 'ref'
		for index,feat in enumerate(ref):
			if not chromPat.match(feat['chrom']):
				continue
			count+=1
			blockCount = int(feat['exonCount'])
			starts = feat['exonStarts'].strip(',').split(',')
			ends = feat['exonEnds'].strip(',').split(',')
			for i in range(blockCount-1):
				tempkey = (feat['chrom'],int(ends[i]),int(starts[i+1]))
				if tempkey not in refDict:
					juncCount += 1
					refDict.add(tempkey)
			#if count%10000 == 0:
			#	print count, ' records processed'

	else:
		ref = readBed.BEDReader(reference)
		print 'bed'
		for index,feat in enumerate(ref):
			#feat = line.split('\t')
			if not chromPat.match(feat['chrom']):
				continue
			count += 1
			blockCount = int(feat['blockCount'])
			starts = feat['blockStarts'].strip(',').split(',')
			sizes = feat['blockSizes'].strip(',').split(',')
			for i in range(blockCount - 1):
				tempstart = int(feat['chromStart']) + int(starts[i]) + int(sizes[i])
				tempend = int(feat['chromStart'])+int(starts[i+1])
				tempkey = (feat['chrom'],tempstart,tempend)
				if tempkey not in refDict:
					juncCount += 1
					refDict.add(tempkey)
					#refDict[tempkey] = [feat,]
			#if count%10000 == 0:
			#	print count, ' records processed'

	#try:
	#	feats = ref.read().strip().split('\n')
	#except IOError:
	#	exit('Empty reference file.')

	#print len(feats)
	print "Total of ", count, " records processed."
	#f = open(reference+'.index','w')
	#for record in refDict.items():
	#	f.write(record[0].split('_')[0] + '\t' + record[0].split('_')[1] +'\n')
	#f.close()
	print 'Total of ', juncCount, ' junctions obtained'
	print 'Finish processing %s'%(reference)
	return refDict


def getSideDict(refDict):
	"""Split the refdict into two sizes to look for junctions that only have one side in exon."""
	leftDict = set()
	rightDict = set()
	for r in refDict:
		leftDict.add((r[0],r[1]))
		rightDict.add((r[0],r[2]))
	return leftDict, rightDict

def toBedRecord(key,feat):
	"""Extract information to create a bed record for output."""
	return {'chrom':key[0],'chromStart':key[1]-3,'chromEnd':key[2]+3,'name':feat['name'],'score':feat['score'],'strand':feat['strand'],'thickStart':key[1]-3,'thickEnd':key[2]+3,'itemRgb':feat['itemRgb'],'blockCount':2,'blockSizes':'3,3','blockStarts':'0,%d'%(key[2]-key[1]+3,)}

def printStat(unique, total, name):
	"""Print the stats for a given type of junctions"""
	print "%d/%d unique/total %s junctions detected."%(unique,total,name)

'''
Find the novel junctions in the inputfile. If the entries in the
input can contain several junctions, then this will use the generateDict()
to generate a dictionary similar to the refDict.
'''
def findNovelJunctions(refDict,inputfiles,mismatch,oneside=False,leftDict=set(),rightDict=set()):
	print 'Now generating junctions---------------------'
	for inputfile in inputfiles:
		print "\n=============",inputfile,"============="
		infile = readBed.BEDReader(inputfile)

		sameJunctions = set()
		newJunctions = set()
		misaJunctions = set()
		onesideJunctions = set()
		#record the number of junctions in each category.
		#this could include multiple count of a single junction due to the repetition in input.
		sameCount = 0
		newCount = 0
		misaCount = 0
		onesideCount = 0
		newwriter = readBed.BEDWriter(inputfile+'.new.bed')
		samewriter = readBed.BEDWriter(inputfile+'.same.bed')
		if mismatch != 0:
			misawriter = readBed.BEDWriter(inputfile+'.misa.bed')
		if oneside:
			onesideWriter = readBed.BEDWriter(inputfile+'.oneside.bed')


		for feat in infile:
			if int(feat['blockCount']) < 2:
				continue
			chromStart = int(feat['chromStart'])
			blockCount = int(feat['blockCount'])
			blockStarts = [int(i) for i in feat['blockStarts'].strip(',').split(',')]
			blockSizes = [int(i) for i in feat['blockSizes'].strip(',').split(',')]
			for i in range(blockCount - 1):
				tempStart = chromStart + blockStarts[i] + blockSizes[i]
				tempEnd = chromStart + blockStarts[i+1]
				tempKeys = []
				firstkey = (feat['chrom'],tempStart,tempEnd)
				exist = -1
				if firstkey in refDict:
					sameCount += 1
					if firstkey not in sameJunctions:
						sameJunctions.add(firstkey)
						samewriter.writerow(toBedRecord(firstkey,feat))
				else:
					misa = False
					for ml in range(-mismatch,mismatch+1):
						for mr in range(-mismatch,mismatch+1):
							if (feat['chrom'],tempStart+ml,tempEnd+mr) in refDict:
								misa = True
								break
						if misa:
							break
					if misa:
						misaCount += 1
						if firstkey not in misaJunctions:
							misaJunctions.add(firstkey)
							if mismatch != 0:
								misawriter.writerow(toBedRecord(firstkey,feat))
					else:
						newCount += 1
						if firstkey not in newJunctions:
							newJunctions.add(firstkey)
							newwriter.writerow(toBedRecord(firstkey,feat))

					if oneside and ((firstkey[0],firstkey[1]) in leftDict or (firstkey[0],firstkey[2]) in rightDict):
						onesideCount += 1
						if firstkey not in onesideJunctions:
							onesideJunctions.add(firstkey)
							onesideWriter.writerow(toBedRecord(firstkey,feat))


		print sameCount + misaCount + newCount, " junctions in total detected"
		print len(sameJunctions) + len(newJunctions) + len(misaJunctions), " unique junctions detected"
		printStat(len(newJunctions),newCount,'novel')
		printStat(len(sameJunctions),sameCount, 'same')
		if mismatch != 0:
			printStat(len(misaJunctions),misaCount,'misa')
		if oneside:
			printStat(len(onesideJunctions),onesideCount,'oneside')


if __name__=='__main__':
	if platform.python_version().startswith('2.7'):
		import argparse
		parser = argparse.ArgumentParser(description="Find the novel junctions.")
		parser.add_argument('reference',help='The reference file should be in bed format or refflat format')
		parser.add_argument('inputfiles',nargs='+',help='The input file that contains the junctions to compare with the reference')
		parser.add_argument('-m','--mismatch',dest='mismatch',type=int,default=0,help='The position difference that can be allowed for each anchor.')
		parser.add_argument('--oneside', dest='oneside',action='store_true',default=False,help='Whether to output junctions with only one side validated in the reference.')
		try:
			options = parser.parse_args()
		except IOError, msg:
			parser.error(str(msg))
		refDict = generateRefDict(options.reference)
		inputfiles = options.inputfiles
	else:
		from optparse import OptionParser
		parser = OptionParser("usage: %prog [options] reference query1 [query2 query3 ...]")
		#parser.add_option('-r','--reference',dest='reference',help='The reference file should be in bed format or refflat format')
		#parser.add_option('-i','--input',dest='inputfile',help='The input file that contains the junctions to compare with the reference')
		parser.add_option('-m','--mismatch',dest='mismatch',type=int,default=0,help='The position difference that can be allowed for each anchor.')
		parser.add_option('--oneside', dest='oneside',action='store_true',default=False,help='Whether to output junctions with only one side validated in the reference.')
		(options,args) = parser.parse_args()
		'''
		if options.reference is None:
			print 'The reference file is missing\n'
			parser.print_help()
			exit(-1)
		if options.inputfile is None:
			print 'The input file is missing\n'
			parser.print_help()
			exit(-1)
		'''
		if len(args) < 2:
			print "Tow few arguments.\n"
			parser.print_help()
			exit(-1)
		refDict = generateRefDict(args[0])
		inputfiles = args[1:]
	if options.oneside:
		leftDict,rightDict=getSideDict(refDict)
		findNovelJunctions(refDict,inputfiles,options.mismatch,options.oneside,leftDict,rightDict)
	else:
		findNovelJunctions(refDict,inputfiles,options.mismatch)

