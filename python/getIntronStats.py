import sys
sys.path.append('.')
import readBed
import pylab as plt


refgeneFile = '../mm9.refFlat.bed'
outFile = '../intronSizes.txt'
if __name__=='__main__':
    introns = {}
    exons = {}
    handler = readBed.BEDReader(refgeneFile)
    for r in handler:
        n = int(r['blockCount'])
        starts = [int(t) for t in r['blockStarts'].replace(',','\t').strip().split('\t')]
        sizes = [int(t) for t in r['blockSizes'].replace(',','\t').strip().split('\t')]
        chromStart = int(r['chromStart'])
        for i in range(n-1):
            tempStart = chromStart + starts[i] + sizes[i]
            tempEnd = chromStart + starts[i+1]
            key = '_'.join([r['chrom'],str(tempStart),str(tempEnd)])
            intronSize = tempEnd - tempStart
            if key not in introns:
                introns[key] = intronSize
    intronSizes = introns.values()
    intronSizes.sort(key=lambda k:(k,))
    out = open(outFile,'w')
    for i in intronSizes:
        out.write(str(i) + '\n')
    out.close()
    n,bins,patches = plt.hist(intronSizes,bins = 10000, normed=True,rwidth=0.8,histtype='bar',label=['Intron Length',])
    plt.close()
    plt.figure(1)
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.title('Intron Size Distribution')
    plt.legend()
    bcenters = 0.5*(bins[1:] + bins[:-1])
    plt.plot(bcenters,n)
    plt.show()
