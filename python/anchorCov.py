import os,sys
import re
import csv
def within(pos1,pos2,mis):
    if pos1-mis<=pos2 and pos1+mis>=pos2:
        return True
    return False

if __name__=='__main__':
    mis = 10
    base = 'rsup.nun.bed'
    infile = '../'+base+'.same.bed'
    infile2 = '../'+base+'.new.bed'
    reader1 = csv.reader(open(infile,'rb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    reader2 = csv.reader(open(infile2,'rb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    covFile = '../../coverage.bed'
    covreader = csv.reader(open(covFile,'rb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    coverage = {}
    for row in covreader:
        if row[0] in coverage:
            coverage[row[0]].append(row)
        else:
            coverage[row[0]] = [row,]
    tjuncs = []
    fjuncs = []
    for row in reader1:
        tjuncs.append(row)
    for row in reader2:
        fjuncs.append(row)
    corr={}
    for i in range(len(tjuncs)):
        for j in range(len(fjuncs)):
            tj = tjuncs[i]
            fj = fjuncs[j]
            if tj[0] == fj[0]:
                if within(int(tj[1]),int(fj[1]),mis) and within(int(tj[2]),int(fj[2]),mis):
                    if i in corr:
                        corr[i].append(j)
                    else:
                        corr[i]=[j,]
    f = open('../'+base+'.ancove','w')
    for key in corr:
        #the interval is this [start,end) where both are 0-based
        currjunc = tjuncs[key]
        start  = int(currjunc[1])+3
        end = int(currjunc[1])+int(currjunc[-1].split(',')[1])
        leftCov = 0
        rightCov = 0
        fstarts = []
        fends = []
        for i in corr[key]:
            fstarts.append(int(fjuncs[i][1])+3)
            fends.append(int(fjuncs[i][1])+int(fjuncs[i][-1].split(',')[1]))
        fleftCovs = []
        frightCovs = []
        for cov in coverage[currjunc[0]]:
            if int(cov[1]) <= start and int(cov[2]) > start:
                leftCov = int(cov[-1])
            if int(cov[1]) <= end and int(cov[2]) > end:
                rightCov = int(cov[-1])
            for i,j in zip(fstarts,fends):
                if int(cov[1]) <= i and int(cov[2])>i:
                    print 'hahahahahahahahahahahahahhaahahahahhahahahahahaha'
                    fleftCovs.append(int(cov[-1]))
                if int(cov[1]) <= j and int(cov[2])>j:
                    frightCovs.append(int(cov[-1]))

        f.write(str(leftCov)+'\t'+str(rightCov))
        for i in fleftCovs:
            f.write('\t'+str(i))
        for i in frightCovs:
            f.write('\t'+str(i))
        f.write('\n')
        print leftCov,' ',rightCov
    f.close()
