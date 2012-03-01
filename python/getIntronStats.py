import os,sys
sys.path.append('.')
import readBed
import pylab as plt
from scipy import arange, optimize, special


refgeneFile = '../mm9.refFlat.bed'
outFile = '../intronSizes.txt'
rawFile = '../intronSizes_raw.txt'
if __name__=='__main__':
    introns = {}
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
    out = open(rawFile,'w')
    for i in intronSizes:
        out.write(str(i) + '\n')
    out.close()
    
    n,bins,patches = plt.hist(intronSizes,bins = 10000, normed=True,rwidth=0.8,histtype='bar',label=['Intron Length',],cumulative=True)
    n1,bins1,patches1 = plt.hist(intronSizes,bins = 10000,rwidth=0.8,histtype='bar',label=['Intron Length',])
    print len(n),' ',len(bins)
    out = open(outFile,'w')
    for i in range(len(n1)):
        out.write("%d\t%d\t%d\t%f\n"%(bins1[i],bins1[i+1],n1[i],n[i]))
    out.close()
    plt.close()
    
    plt.figure(1)
    plt.title('Intron Size Distribution')
    #plt.hist(intronSizes,bins = 500, normed=True,rwidth=0.8,histtype='bar',label=['Intron Length',],cumulative=True)
    plt.legend()
    bcenters = 0.5*(bins[1:] + bins[:-1])
    plt.plot(bcenters,n)
    plt.xlabel('Length')
    plt.ylabel('Frequency')
    plt.xlim((0,160000))
    plt.savefig("intronStats.png")
    #plt.show()