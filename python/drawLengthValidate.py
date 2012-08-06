import sys,os
sys.path.append('.')
import readBed
import numpy as np
from pylab import *


if __name__=='__main__':
    inputdir = '../outputs/rawbin.splicemap.4'
    refSameFile = os.path.join(inputdir,'rawbin.refgene.same.bed')
    refNewFile = os.path.join(inputdir,'rawbin.refgene.new.bed')
    estSameFile = os.path.join(inputdir,'rawbin.est.same.bed')
    estNewFile = os.path.join(inputdir,'rawbin.est.new.bed')

    #Against refgene
    sameHandle = readBed.BEDReader(refSameFile)
    newHandle = readBed.BEDReader(refNewFile)
    sameJunc = []
    newJunc = []
    for i in sameHandle:
        sameJunc.append(i)
    for i in newHandle:
        newJunc.append(i)
    sameLen = []
    newLen = []
    allLen = []
    for i in sameJunc:
        tempLength = int(i['blockStarts'].split(',')[1])-int(i['blockSizes'].split(',')[0])
        sameLen.append(tempLength)
        allLen.append((tempLength,True))
    for i in newJunc:
        tempLength = int(i['blockStarts'].split(',')[1])-int(i['blockSizes'].split(',')[0])
        newLen.append(tempLength)
        allLen.append((tempLength,False))
    allLen.sort(key=lambda k:(k[0],))
    print len(allLen)
    groupby100x = []
    groupby100y = []
    i = 0
    print 'test'
    groupSize = 50
    while i < len(allLen):
        tempList = allLen[i:(i+groupSize)]
        tempTotal = 0
        tempCountTrue=0
        for l in tempList:
            tempTotal += l[0]
            if l[1]:
                tempCountTrue+=1
        tempMean = float(tempTotal)/groupSize
        tempRate = tempCountTrue/float(groupSize)
        groupby100x.append(tempMean)
        groupby100y.append(tempRate)
        i+=groupSize
    figure(1)
    plot(groupby100x,groupby100y,'go')
    xlabel('Length (bp)')
    ylabel('Validation rate')
    figure(2)
    bins = np.arange(0,400001,5000)
    sn,sbins,spatches = hist(sameLen,bins=bins)
    nn,nbins,npatches = hist(newLen,bins=bins)
    figure(3)
    valRate = []
    for i in range(len(sn)):
        if sn[i]+nn[i]<20:
            valRate.append(-0.1)
        else:
            valRate.append(float(sn[i])/(sn[i]+nn[i]))
    for i in range(len(valRate)):
        if valRate[i] < 0:
            if i != 0 and i != len(valRate)-1:
                valRate[i] = 0.5*(valRate[i-1]+valRate[i+1])
    bcenters = 0.5*(bins[1:]+bins[:-1])
    plot(bcenters,valRate,'r^')
    show()
