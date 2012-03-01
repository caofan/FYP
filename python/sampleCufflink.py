import os,sys
import numpy.random as npr

def generateReads(numReads,outFileName):
    fromEachFile = numReads / 10
    out = open(outFileName+"_1.fastq", "w")
    out2 = open(outFileName+"_2.fastq", "w")
    for i in range(45,55):
        filename = "SRR0379"+str(i)+"_1.fastq"
        filename2 = "SRR0379"+str(i)+"_2.fastq"
        infile = open(filename)
        readsInFile = 0
        for r in infile:
            readsInFile+=1
        infile.close()
        readsInFile /= 4
        numbers = set()
        while len(numbers) < fromEachFile:
            tempNumbers = npr.random_integers(1,readsInFile+1, fromEachFile-len(numbers)+100)
            numbers = numbers.union(set(tempNumbers))
        numbers = list(numbers)
        numbers.sort();
        infile1 = open(filename)
        infile2 = open(filename2)
        count = 0
        numTrack = 0
        rowsWritten = 0
        lastRowsWritten = 0
        for r,k in zip(infile1,infile2):
            if count/4 +1 == numbers[numTrack]:
                out.write(r)
                out2.write(k)
                rowsWritten+=1
            count += 1
            if rowsWritten != lastRowsWritten and rowsWritten % 4 == 0:
                lastRowsWritten = rowsWritten
                numTrack+=1
            if numTrack >= len(numbers):
                break
        print rowsWritten, ' from ', "SRR0379"+str(i)
    out.close()
    out2.close()

if __name__=="__main__":
    generateReads(20164140/2,"cuff.paired.small")
    generateReads(20164140,"cuff.paired.large")


