import os,sys
import re
from collections import defaultdict

if __name__=='__main__':
	if len(sys.argv) < 2:
		print "Usage: python getReadStats.py cig_file"
		print "Please supply reads file in cig format"
		sys.exit(1)
	infile = open(sys.argv[1])
	juncsInReads = defaultdict(int)
	minExon = defaultdict(int)
	for r in infile:
		cigar = r.strip().split('\t')[3]
		juncCount = cigar.count('N')
		if juncCount > 0:
			juncsInReads
			min=1000
			iter = re.finditer("((?P<m>\d+)M(\d+I)*\d+N)|((?P<n>\d+)M(\d+I)*$)",cigar,flags=re.IGNORECASE)
			for m in iter:
				if m.group('m') and int(m.group('m')) < min:
					min = int(m.group('m'))
				if m.group('n') and int(m.group('n')) < min:
					min = int(m.group('n'))
			juncsInReads[juncCount] += 1
			if min >= 40:
				print r
			minExon[min] += 1
			
	for key in juncsInReads:
		print key,': ',juncsInReads[key]
	
	print "\n\n"

	for key in minExon:
		print key,": ",minExon[key]
