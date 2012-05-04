import csv
class CommentedFileReader:
    def __init__(self,f,commentstring='#'):
        self.f=open(f,'rU')
        self.commentstring = commentstring
    def next(self):
        line = self.f.next()
        while line.startswith(self.commentstring):
            line = self.f.next()
        return line
    def __iter__(self):
        return self

csv.register_dialect('refFlat',delimiter='\t',
                     quoting=csv.QUOTE_NONE,
                     skipinitialspace = True)

class FlatReader(csv.DictReader):
    FIELDS=('geneName','name','chrom',
            'strand','txStart','txEnd',
            'cdsStart','cdsEnd','exonCount',
            'exonStarts','exonEnds')

    def __init__(self,filename):
        csv.DictReader.__init__(self,CommentedFileReader(filename),
                                dialect='refFlat',fieldnames=self.FIELDS)

class FlatWriter(csv.DictWriter):
    FIELDS=('geneName','name','chrom',
            'strand','txStart','txEnd',
            'cdsStart','cdsEnd','exonCount',
            'exonStarts','exonEnds')

    def __init__(self,filename):
        csv.DictWriter.__init__(self,open(filename,'wb'),
                                dialect='refFlat',fieldnames=self.FIELDS)

        
