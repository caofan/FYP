#!/usr/bin/python

import os,sys
sys.path.append('.')
import readBed
exclude = 'chr13'
stats = {}
RESIDUE = 10
GAPS = [400000,200000,100000,50000,20000,]
outdir = 'supposedVSfound/rawbin/simulated-testResidue'
supposedDir = 'supposedJuncs'
outFile = 'residue_dist_stats.txt'
if __name__=='__main__':
    STATS = {}
    supposedTotal = 0
    statsOut = open(outFile,'w')
    for filename in os.walk(supposedDir).next()[2]:
        if not filename.startswith('chr13'):
            tempHanlder = open(os.path.join(supposedDir,filename))
            for r in tempHanlder:
                supposedTotal += 1
    while RESIDUE > 4:
        for gap in GAPS:
            handler = open('rawbin/rnaseq-test/src/rnaseq.cpp')
            out = open('rnaseq.temp','w')
            for r in handler:
                if r.startswith('const int RESIDUE'):
                    r = 'const int RESIDUE=%d;\n'%(RESIDUE,)
                elif r.startswith('const int MAXMULTIGAP'):
                    r = 'const int MAXMULTIGAP=%d;\n'%(gap,)
                out.write(r)
            out.close()
            handler.close()
            os.system('rm rawbin/rnaseq-test/src/rnaseq.cpp')
            os.system('mv rnaseq.temp rawbin/rnaseq-test/src/rnaseq.cpp')
            os.system('make -C rawbin/rnaseq-test')
            os.system('python runRawbin_test.py')
            os.system('python ../python/runFoundvsSupposed.py')
            totalNew = 0
            totalSame = 0
            for filename in os.walk(outdir).next()[2]:
                if filename.endswith('same.bed'):
                    tempHandler = open(os.path.join(outdir,filename))
                    for r in tempHanlder:
                        totalSame += 1
                elif filename.endswith('new.bed'):
                    tempHandler = open(os.path.join(outdir,filename))
                    for r in tempHanlder:
                        totalNew += 1
            STATS['_'.join([str(RESIDUE),str(gap),])] = supposed
            statsOut.write('%d\t%d\t%d\t%d\t%d\n'%(RESIDUE,gap,supposedTotal,totalSame,totalNew,))
        RESIDUE-=1
    statsOut.close()
