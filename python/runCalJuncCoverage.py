import os,sys
import calJuncCoverage
rootdir = "/home/Ken/caofan/outputs"
juncdir = "/home/Ken/caofan/supposedVSfound/rawbin/simulated-19SEP-17SEP-test-12"

outputdir = "/home/Ken/caofan/supposedVSfound/rawbin/simulated-19SEP-17SEP-test-12-calOut"
if not os.path.exists(outputdir):
    os.system("mkdir -p %s"%(outputdir,))
index=[]
index = range(1,20)
index.append("X")
index.append("Y")

signals = {}

BASE = {"A":"T","C":"G","G":"C","T":"A"}


#Rawbin, SpliceMap, Tophat, Simulated
noSupport = [0,0,0,0]
noCov = [0,0,0,0]
noSignal = [0,0,0,0]
AB = [0,0,0,0]
BC = [0,0,0,0]
CA = [0,0,0,0]
ABC = [0,0,0,0]


for i in index:
    currJuncFile = os.path.join(juncdir,'supposed.rawbin.chr%s.fasta.new.bed'%(str(i),))
    os.system("echo %s"%("chr"+str(i),))
    rawbinReadFile = os.path.join(rootdir,"rawbin/simulated-19SEP-17SEP-test-12/chr%s.fasta/Mapped_modified.bam"%(str(i),))
    spliceReadFile = os.path.join(rootdir,"splicemap/simulated/chr%s.fasta/good_hits.bam"%(str(i),))
    tophatReadFile = os.path.join(rootdir,"tophat/simulated/chr%s.fasta/accepted_hits.bam"%(str(i),))
    simulatedReadFile = os.path.join("/home/Ken/caofan/inputs/simulated_bam/chr%s.bam"%(str(i),))

    rawbinoutfile = os.path.join(outputdir,"rawbin_chr%s.txt"%(str(i),))
    # Deal with Rawbin reads first.
    os.system("echo rawbin")
    os.system("python26 calJuncCoverage.py -s -j %s -r %s -a %d -c %d -o %s"%(currJuncFile, rawbinReadFile,10,15,rawbinoutfile))
    # Deal with Splicemap reads
    spliceoutfile = os.path.join(outputdir,"splice_chr%s.txt"%(str(i),))
    os.system("echo splicemap")
    os.system("python26 calJuncCoverage.py -s -j %s -r %s -a %d -c %d -o %s"%(currJuncFile, spliceReadFile,10,15,spliceoutfile))
    tophatoutfile = os.path.join(outputdir,"tophat_chr%s.txt"%(str(i),))
    os.system("echo tophat")
    os.system("python26 calJuncCoverage.py -s -j %s -r %s -a %d -c %d -o %s"%(currJuncFile, tophatReadFile,10,15,tophatoutfile))
    simulatedoutfile = os.path.join(outputdir,"simulated_chr%s.txt"%(str(i),))
    os.system("echo original_reads")
    os.system("python26 calJuncCoverage.py -s -j %s -r %s -a %d -c %d -o %s"%(currJuncFile, simulatedReadFile,0,15,simulatedoutfile))
    tempjunctions,tempsignals = calJuncCoverage.getJunctions(currJuncFile,True)
    for sig in tempsignals:
        sigRev = ""
        for base in sig:
            sigRev = BASE[base] + sigRev
        if sig in signals:
            signals[sig] += tempsignals[sig]
        elif sigRev in signals:
            signals[sigRev] += tempsignals[sig]
        else:
            signals[sig] = tempsignals[sig]
    os.system("echo ")

print "\n----------------------Summary---------------------------"
for key in signals:
    print key,' ',signals[key]


inputfile = 'run'+os.path.basename(outputdir) + '.txt'
inputfile = inputfile.replace('runc','runC')
f = open(inputfile)
r = f.next()
while True:
    if r.startswith('chr'):
        f.next()
        f.next()
        rawbinStats = [int(key) for key in f.next().split('\t')]
        f.next()
        f.next()
        spliceStats = [int(key) for key in f.next().split('\t')]
        f.next()
        f.next()
        tophatStats = [int(key) for key in f.next().split('\t')]
        f.next()
        f.next()
        simulatedStats = [int(key) for key in f.next().split('\t')]

        noSupport[0] += rawbinStats[0]
        noCov[0] += rawbinStats[1]
        noSignal[0] += rawbinStats[2]
        AB[0] += rawbinStats[3]
        BC[0] += rawbinStats[4]
        CA[0] += rawbinStats[5]
        ABC[0] += rawbinStats[6]

        noSupport[1] += spliceStats[0]
        noCov[1] += spliceStats[1]
        noSignal[1] += spliceStats[2]
        AB[1] += spliceStats[3]
        BC[1] += spliceStats[4]
        CA[1] += spliceStats[5]
        ABC[1] += spliceStats[6]

        noSupport[2] += tophatStats[0]
        noCov[2] += tophatStats[1]
        noSignal[2] += tophatStats[2]
        AB[2] += tophatStats[3]
        BC[2] += tophatStats[4]
        CA[2] += tophatStats[5]
        ABC[2] += tophatStats[6]


        noSupport[3] += simulatedStats[0]
        noCov[3] += simulatedStats[1]
        noSignal[3] += simulatedStats[2]
        AB[3] += simulatedStats[3]
        BC[3] += simulatedStats[4]
        CA[3] += simulatedStats[5]
        ABC[3] += simulatedStats[6]
    try:
        r = f.next()
    except StopIteration:
        break
f.close()
print '\n'
print "Mapper\tNo_Support\tNoCov\tNoSignal\tAB\tBC\tCA\tABC"
print 'Rawbin\t%d\t%d\t%d\t%d\t%d\t%d\t%d'%(noSupport[0],noCov[0],noSignal[0],AB[0],BC[0],CA[0],ABC[0],)
print 'SpliceMap\t%d\t%d\t%d\t%d\t%d\t%d\t%d'%(noSupport[1],noCov[1],noSignal[1],AB[1],BC[1],CA[1],ABC[1],)
print 'Tophat\t%d\t%d\t%d\t%d\t%d\t%d\t%d'%(noSupport[2],noCov[2],noSignal[2],AB[2],BC[2],CA[2],ABC[2],)
print 'Simulated\t%d\t%d\t%d\t%d\t%d\t%d\t%d'%(noSupport[3],noCov[3],noSignal[3],AB[3],BC[3],CA[3],ABC[3],)
