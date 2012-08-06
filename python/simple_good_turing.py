import os,sys
import re
import csv
import numpy as np
from pylab import *

def smooth(data):
    Z = {}
    x_values = []
    y_values = []
    rawx = []
    rawy = []
    for i in range(len(data)):
        item = data[i]
        if i != 0 and i != len(data)-1:
            Z[item[0]] = 2.0*item[1]/(data[i+1][0]-data[i-1][0])
            x_values.append(np.log10(item[0]))
            y_values.append(np.log10(Z[item[0]]))
            rawx.append(item[0])
            rawy.append(Z[item[0]])
        elif i == 0:
            Z[item[0]] = 2.0*item[1]/(data[i+1][0]-0)
            x_values.append(np.log10(item[0]))
            y_values.append(np.log10(item[1]))
            rawx.append(item[0])
            rawy.append(item[1])
        elif i == len(data)-1:
            Z[item[0]] = 2.0*item[1]/(2*item[0]-2*data[i-1][0])
            rawx.append(item[0])
            rawy.append(Z[item[0]])
    #print x_values
    #plot
    #subplot(121)
    #scatter(x_values,y_values,marker='o',c='r')
    m = np.polyfit(x_values,y_values,5)
    #yfit = np.polyval(m,x_values)
    #plot(x_values,yfit,'b--')
    #subplot(122)
    #scatter(rawx,rawy)
    #m1 = np.polyfit(rawx,rawy,20)
    #yfit1 = np.polyval(m1,rawx)
    #plot(rawx,yfit1,'r--')
    #show()
    return Z,m


def normalize(data,output,rawDataOut,isSecond):
    for i in range(shape(data)[1]):
        col = data[:,i]
        print np.shape(col)
        countDict = {}
        zeroCount = 0
        for v in col:
            if v == 0:
                zeroCount += 1
            elif v in countDict:
                countDict[v] += 1
            else:
                countDict[v] = 1
        N = 0
        for index in countDict:
            print index
            N += index*countDict[index]
        print N
        print sum(col)
        P0 = countDict[1] / float(N)

        transDict = {}
        simpleTransDict = {}
        countItems = countDict.items()
        countItems.sort(key=lambda k:(k[0]))
        Z,m = smooth(countItems)
        calculateX = True
        for index in range(len(countItems)-1):
            item = countItems[index]
            r = item[0]
            nr = float(item[1])
            nr_1 = float(countItems[index+1][1])
            sr = 10**np.polyval(m,np.log10(r))
            sr_1 = 10**np.polyval(m,np.log10(r+1))
            xr = 0
            if calculateX:
                xr = (r+1)*nr_1/nr
            yr = (r+1)*sr_1/sr
            difference = abs(xr-yr)
            if not calculateX:
                transDict[item] = int(round(yr))
                simpleTransDict[item[0]] = int(round(yr))
            elif calculateX and difference > 1.96*np.sqrt((r+1)**2*nr_1*(1+nr_1/nr)/nr**2):
                transDict[item] = int(round(xr))
                simpleTransDict[item[0]] = int(round(xr))
            else:
                calculateX = False
                print index
                transDict[item] = int(round(yr))
                simpleTransDict[item[0]] = int(round(yr))
        transItems = transDict.items()
        transItems.sort(key=lambda k:(k[1]))
        newCountDict = {}
        for item in transItems:
            if item[1] in newCountDict:
                newCountDict[item[1]] += item[0][1]
            else:
                newCountDict[item[1]] = item[0][1]

        for j in range(len(col)):
            if col[j] == 0:
                if isSecond:
                    output[j][i+5] = 0
                    rawDataOut[j][i+5] = col[j]
                else:
                    output[j][i] = 0
                    rawDataOut[j][i] = col[j]
            elif not (col[j] in simpleTransDict):
                if isSecond:
                    output[j][i+5] = col[j]
                    rawDataOut[j][i+5] = col[j]
                else:
                    output[j][i] = col[j]
                    rawDataOut[j][i] = col[j]
            else:
                if isSecond:
                    output[j][i+5] = simpleTransDict[col[j]]
                    rawDataOut[j][i+5] = col[j]
                else:
                    output[j][i] = simpleTransDict[col[j]]
                    rawDataOut[j][i] = col[j]

if __name__=='__main__':
    inputdir = ''
    ouptutdir = ''
    data = np.genfromtxt('/home/caofan/Downloads/simulated_Data.txt')
    #data1 = np.genfromtxt('../simulate/cond1.8.0002.txt')
    data1 = data[:,0:5]
    data2 = data[:,5:10]
    print np.shape(data1)
    #data2 = np.genfromtxt('../simulate/cond2.8.0002.txt')
    output = np.zeros([np.shape(data1)[0],10])
    rawDataOut = np.zeros([np.shape(data1)[0],10])
    normalize(data1,output,rawDataOut,False)

    normalize(data2,output,rawDataOut,True)

    out = open("../a.txt","w")
    for r in output:
        for c in r:
            out.write(str(int(c))+"\t")
        out.write("\n")
    out.close()
    out = open("../a_rawOut.txt","w")
    for r in rawDataOut:
        for c in r:
            out.write(str(int(c))+"\t")
        out.write("\n")
    out.close()
