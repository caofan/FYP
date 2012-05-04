import os,sys
import readBed
import re

def removeRandom(records):
	result = []
	pat = re.compile(r'^chr(\d+|X|Y|M)$')
	for r in records:
		if pat.match(r['chrom']):
			result.append(r)
		else:
			print r['chrom']
			continue
	return result


class ExonNode:
	"""docstring for ExonNode"""
	def __init__(self,tempKey):
		self.key = tempKey
		self.nexts = set()
		self.prevCount = 0
		self.reached = 0

	def addNext(self,next):
		initLength = len(self.nexts)
		self.nexts.add(next)
		return len(self.nexts)-initLength > 0

	def addPrev(self):
		self.prevCount+=1

	def getPrev(self):
		return self.prevCount

	def getNexts(self):
		return self.nexts

	def setReached(self):
		if self.prevCount > 0:
			self.reached += 1

	def proceed(self):
		assert self.reached <= self.prevCount
		return self.reached == self.prevCount

	def hasNext(self):
		return len(self.nexts) > 0

	def getKey(self):
		"""docstring for getKey"""
		return self.key

def loadExons(records):
	currName = ''
	exons = {}
	tempCount = 0
	for r in records:
		if currName == '' or currName != r['name']:
			currName = r['name']
			tempCount = 1
		else:
			tempCount += 1
		blockSizes = [int(bs) for bs in r['blockSizes'].replace(',','\t').strip().split('\t')]
		blockStarts = [int(bs) for bs in r['blockStarts'].replace(',','\t').strip().split('\t')]
		lastExonKey = '$'
		for bsi, bst in zip(blockSizes, blockStarts):
			tempStart = int(r['chromStart'])+bst
			tempEnd = tempStart + bsi
			tempKey = ':'.join([r['chrom'],str(tempStart),str(tempEnd),r['strand']])

			if tempKey not in exons:
				exons[tempKey] = ExonNode(tempKey)
			if lastExonKey != '$':
				if exons[lastExonKey].addNext(tempKey):
					exons[tempKey].addPrev()
			lastExonKey = tempKey
		#if lastExonKey != '' and '<' not in exons[lastExonKey]:
		#	exons[lastExonKey].append('<')
	return exons

def writeExon(out,exon,exons):
	"""Recursively write exons to file. The writing follow the strategy that a
	node has to wait all its preceding nodes reach it before preceding.
	Arguments:
	"""
	exonChrom,exonStart,exonEnd,exonStrand = exon.split(':')
	exonSize = str(int(exonEnd) - int(exonStart))
	exonNode = exons[exon]
	#print exonNode.getKey()
	exonNode.setReached()
	if exonNode.proceed():
		if exonNode.hasNext():
			nextExons = list(exonNode.getNexts())
			nextExons.sort(key=lambda k:(int(k.split(':')[1]),int(k.split(':')[2])))
			for ne in nextExons:
				exonName = ','.join([exon,ne])
				out.write('\t'.join([exonChrom,exonStart,exonEnd,exonName,'0',exonStrand,exonStart,exonEnd,'0','1',exonSize,'0\n']))
				writeExon(out,ne,exons)
		else:
			exonName = ','.join([exon,'<'])
			out.write('\t'.join([exonChrom,exonStart,exonEnd,exonName,'0',exonStrand,exonStart,exonEnd,'0','1',exonSize,'0\n']))
		exons.pop(exon)
	else:
		exonName = exon
		out.write('\t'.join([exonChrom,exonStart,exonEnd,exonName,'0',exonStrand,exonStart,exonEnd,'0','1',exonSize,'0\n']))



# The writing of the exons should follow the strategy that a node has to wait
# all its preceding exons reach before proceed.
def writeToFile(exons, out):
	exonList = exons.keys()
	exonList.sort(key=lambda
			k:(k.split(':')[0],int(k.split(':')[1]),int(k.split(':')[2])))
	for exon in exonList:
		if exon in exons:
			writeExon(out,exon,exons)
	out.close()


if __name__=='__main__':
	if len(sys.argv) != 3 && len(sys.argv) != 5:
		print "Usage: python getExons.py bedInput bedOutput <fastaIn> <fastaOut>\n"
		sys.exit(1)
	records = readBed.BEDReader(sys.argv[1])
	outfile = open(sys.argv[2],'w')
	result = removeRandom(records)
	exons = loadExons(result)
	writeToFile(exons, outfile)
	if len(sys.argv) == 5:
		os.system("fastaFromBed -fi %s -bed %s -fo %s -name"%(sys.argv[3],sys.argv[2],sys.argv[4]))
