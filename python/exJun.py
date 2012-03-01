'''
This program extracts the introns that are found between the junctions.
'''

import os,sys
sys.path.append('.')
import readBed
import csv
import argparse
import re
parser = argparse.ArgumentParser(description="Options")
parser.add_argument('-i','--inputdir',dest='inputdir')
parser.add_argument('-o','--outputdir',dest='outdir')
parser.add_argument('-c','--cutoff',dest='cutoff',help='Output the junctions that spans an intron larger than the cutoff. If 0, then no output is given.',type=int, default=0)
def extractIntron(filename):
    reader = readBed.BEDReader(filename)
    try:
        line = reader.next()
    except StopIteration:
        line = None
    output = []
    exist = dict()
    lengths = dict()
    largejunctions = []
    #if not line['blockStarts']:
    #    line = reader.next()
    i = 0
    print options.cutoff
    while line:
        i+=1
        blocksize = int(line['blockSizes'].split(',')[0])
        tempstart = int(line['chromStart']) + blocksize
        tempend = int(line['blockStarts'].split(',')[-1])+int(line['chromStart'])
        templength = tempend - tempstart
        if options.cutoff > 0 and templength > options.cutoff:
            print type(line)
            largejunctions.append(line)
        if templength in lengths.keys():
            lengths[templength] += 1
        else:
            lengths[templength] = 1
        tempout = (line['chrom'],tempstart,tempend,line['name'],line['score'],line['strand'],tempstart,tempend,line['itemRgb'],1,tempend-tempstart,0)

        try:
            exist[str(tempstart)+' '+str(tempend)]
        except KeyError:
            output.append(tempout)
            exist[str(tempstart)+' '+str(tempend)]=1
        try:
            line = reader.next()
        except StopIteration:
            break
    print 'Writing ',os.path.basename(filename),'-----------------------------------------'
    writer = csv.writer(open(os.path.join(options.outdir,os.path.basename(filename)+'.intron'),'wb'),dialect='bed')
    writer.writerows(output)
    f = open(os.path.join(options.outdir,os.path.basename(filename)+'.length'),'w')
    fi = open(os.path.join(options.outdir,os.path.basename(filename)+'.length.raw'),'w')
    items = lengths.items()
    items.sort(key=lambda k:(k[0]))
    for item in items:
        f.write(str(item[0]) + '\t' + str(item[1])+'\n')
        for i in range(item[1]):
            fi.write(str(item[0]))
            fi.write('\n')
    f.close()
    fi.close()
    if options.cutoff > 0:
        largewriter = readBed.BEDWriter(os.path.join(options.outdir,os.path.basename(filename)+'.large'))
        largewriter.writerows(largejunctions)


if __name__=="__main__":
    try:
        global options
        options = parser.parse_args()
        if not (options.inputdir and options.outdir):
            exit('Please specify the input and output directories')
    except IOError, msg:
        parser.error(str(msg))
    inpath = options.inputdir
    for filename in os.walk(inpath).next()[2]:
        if filename.endswith('genefile'):
            continue
        elif not re.match(r'^.+\..+', filename):
            continue
        else:
            extractIntron(os.path.join(inpath,filename))
