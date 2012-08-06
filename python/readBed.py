#code from http://www.moosechips.com/2011/01/python-chip-seq-bed-file-reader/
#/usr/bin/python

import csv
class CommentedFileReader:
    """
    Helper class for file reading.
    Skips lines starting with '#'

    tsv_file = csv.reader(CommentedFileReader("inputfile.txt"),
                      delimiter='\t')
    for row in tsv_file:
        print row[2] # prints column 3 of each line
    """
    def __init__(self, f, commentstring="#"):
        self.f = open(f, 'rU')
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring) or line.startswith('track'):
            line = self.f.next()
        return line
    def __iter__(self):
        return self

csv.register_dialect('bed', delimiter = '\t',quoting = csv.QUOTE_NONE,skipinitialspace = True)

class BEDReader(csv.DictReader):
    """
    Read BED files into a DictReader.
    See BEDReader.FIELDS for field names

    Example:
        bed = BEDReader("file.bed")
        for line in bed:
            # print the chromStart
            print(line['chromStart'])
    """
    FIELDS = ('chrom', 'chromStart', 'chromEnd',
              'name', 'score', 'strand',
              'thickStart', 'thickEnd',
              'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts')

    def __init__(self, filename):
        csv.DictReader.__init__(self, CommentedFileReader(filename), dialect='bed',
                                fieldnames=self.FIELDS)


class BEDWriter(csv.DictWriter):
    FIELDS = ('chrom', 'chromStart', 'chromEnd',
              'name', 'score', 'strand',
              'thickStart', 'thickEnd',
              'itemRgb',
              'blockCount', 'blockSizes', 'blockStarts')

    def __init__(self, filename):
        csv.DictWriter.__init__(self, open(filename,'wb'), dialect='bed',fieldnames=self.FIELDS,extrasaction='ignore',quoting=csv.QUOTE_NONE)


if __name__ == "__main__":
    bed = BEDReader("4/rawbin.refgene")
    for line in bed:
        print(line['chromStart'])
