import sys
import argparse
sys.path.append('.')
import readBed
import string
import pylab as plt


def loadons(reffile):
    '''
    Load exons and introns from the annotation file.
    Return two dicts, exons and introns, where the keys
    are the sizes and the values are sets of intons/exons haveing the corresponding size.
    '''
    handler = readBed.BEDReader(reffile)
    introns = {}
    exons = {}
    for r in handler:
        n = int(r['blockCount'])
        starts = [int(t) for t in r['blockStarts'].strip(string.whitespace+",").split(',')]
        sizes = [int(t) for t in r['blockSizes'].strip(string.whitespace+",").split(',')]
        chromStart = int(r['chromStart'])
        #Load introns
        for i in range(n):
            tempStart = chromStart + starts[i]
            tempEnd = tempStart + sizes[i]
            tempkey = '_'.join([r['chrom'],str(tempStart),str(tempEnd)])
            if sizes[i] not in exons:
                exons[sizes[i]] = set([tempkey,])
            else:
                exons[sizes[i]].add(tempkey)
        #Load introns
        for i in range(n-1):
            tempStart = chromStart + starts[i] + sizes[i]
            tempEnd = chromStart + starts[i+1]
            tempkey = '_'.join([r['chrom'],str(tempStart),str(tempEnd)])
            intronSize = tempEnd - tempStart
            if intronSize not in introns:
                introns[intronSize] = set([tempkey,])
            else:
                introns[intronSize].add(tempkey)
    return exons, introns

def plot(data,descr):
    '''
    Plot the distribution of the input data.
    '''
    n,bins,patches = plt.hist(data,bins = 10000, normed=True,rwidth=0.8,histtype='bar',label=['%s Length'%descr,])
    plt.close()
    plt.figure(1)
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.title('%s Size Distribution'%descr)
    plt.legend()
    bcenters = 0.5*(bins[1:] + bins[:-1])
    plt.plot(bcenters,n)
    plt.show()

def getNumSmallerThan(data,s):
    '''
    Return the number of items in data with size smaller than the size as specified in s.
    The argument data is a dict where the keys are the sizes and values are sets/lists.
    '''
    keys = data.keys()
    keys.sort()
    n = 0
    for k in keys:
        if k <= s:
            n += len(data[k])
    return n


def getNumLargerThan(data,s):
    '''
    Return the number of items in data with size larger than the size as specified in s.
    The argument data is a dict where the keys are the sizes and values are sets/lists.
    '''
    temp = {}
    for i in data:
        temp[-i] = data[i]
    return getNumSmallerThan(temp,-s)

if __name__=="__main__":
    parser = argparse.ArgumentParser(description="Generate some statistics.")
    parser.add_argument('ref',metavar='R',help="Annotation file.")
    parser.add_argument('-s',dest='s',nargs='+',type=int,help="Smaller or equal than.")
    parser.add_argument('-l',dest='l',nargs='+',type=int,help="Larger or equal than.")
    args = parser.parse_args()
    exons,introns = loadons(args.ref)
    et = getNumSmallerThan(exons,sys.maxint)
    it = getNumSmallerThan(introns,sys.maxint)
    print "Total exons: %d\nTotal introns: %d"%(et,it)
    if args.s != None:
        print "Smaller than:"
        for n in args.s:
            es = getNumSmallerThan(exons,n)
            ins = getNumSmallerThan(introns,n)
            print '\t%d\t%d\t%f\t%d\t%f'%(n,es,es*1.0/et,ins,ins*1.0/it)
    if args.l != None:
        print "Larger than:"
        for n in args.l:
            el = getNumLargerThan(exons,n)
            il = getNumLargerThan(introns,n)
            print '\t%d\t%d\t%f\t%d\t%f'%(n,el,el*1.0/et,il,il*1.0/it)




