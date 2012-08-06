import os
import re
import sys
sys.path.append('.')
import csv
import readBed
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

parser = argparse.ArgumentParser(description="Options")
parser.add_argument('-r', dest='recursive', action='store_true', default=False, help="Recursively look at the subdirectories")
parser.add_argument('-p', '--path', dest='directory', required=True, help="The directory where the files are located.")
parser.add_argument('-l', '--length', dest='readlength', type=int, default=0)
parser.add_argument('-o', '--output', dest='outputdir', default='.',help="The directory to store the output files.")
if __name__ == "__main__":
    try:
        global options
        options = parser.parse_args()
    except IOError, msg:
        parser.error(str(msg))
    counter = 0
    for currdir in os.walk(options.directory):
        if (not options.recursive) and counter > 0:
            break
        counter += 1
        for filename in currdir[2]:
            if re.match('.+\.bed$', filename):
                print 'Processiong ', filename, '-'*20
                reader = readBed.BEDReader(os.path.join(currdir[0],filename))
                fasta = SeqIO.parse(open(os.path.join(currdir[0],filename.split('.')[0] + '.fasta')), 'fasta')
                records = []
                fastahandle = open(os.path.join(options.outputdir,filename.split('.')[0]+'.fasta'),'w')
                output = []
                for line in reader:
                    currSeqRecord = fasta.next()
                    if line['chrom'] == 'polyA':
                        pass
                    elif options.readlength != 0 and len(currSeqRecord.seq) < options.readlength:
                        pass
                    else:
                        tempout = (line['chrom'], line['chromStart'], line['chromEnd'], line['name'], line['score'], line['strand'], line['thickStart'], line['thickEnd'], line['itemRgb'], line['blockCount'], line['blockSizes'], line['blockStarts'])
                        output.append(tempout)
                        FastaWriter(fastahandle,wrap=80).write_file([currSeqRecord,])
                fastahandle.close()
                writer = csv.writer(open(os.path.join(options.outputdir, filename), 'wb'), dialect='bed')
                writer.writerows(output)
