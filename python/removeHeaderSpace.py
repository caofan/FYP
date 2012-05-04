import os,sys
import argparse

def process(filename,typeSize):
	f = open(filename)
	currFile = filename.split('/')[-1]
	hot = (currFile.split('.')[-2]).split('_')[-1]
	out = open(filename+'.out','w')
	i = 0
	currHeader = ""
	for l in f:
		outStr = l.strip()
		if i%typeSize == 0:
			items = l.strip().split()
			outStr = '_'.join(items)
			outStr = '/'.join([outStr,hot])
			currHeader = outStr[1:]
		elif i%typeSize == 2:
			outStr = "+" + currHeader
		out.write(outStr+'\n')
		i+=1
	f.close()
	out.close()

def processFile(filename):
	try:
		extension = filename.split('.')[-1]
	except IndexError:
		print "File type error"
		return
	if extension in ['fasta','fa']:
		process(filename,2)
	elif extension in ['fastq']:
		process(filename,4)

if __name__=='__main__':
	parser = argparse.ArgumentParser(description="Remove whitespace in fastq/fasta headers.")
	parser.add_argument('inputs',nargs='+',help='input files')
	#parser.add_argument('type',nargs=1,choices=["fastq","fasta"],help='type of input')
	args = parser.parse_args()
	for filename in args.inputs:
		processFile(filename)
