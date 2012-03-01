import os,sys
import readBed

def removeRandom(records):
	result = []
	for r in records:
		if 'random' in r['chrom']:
			print r['chrom']
			continue
		else:
			result.append(r)
	return result

def loadExons(records):
    currName = ''
    exons = {}
    tempCount = 0
    for r in records:
        if currName == '' or currName != r['name']:
            currName = r['name']
            tempCount = 1
        else:
            tempCount += 1
        blockSizes = [int(bs) for bs in r['blockSizes'].replace(',','\t').strip().split('\t')]
        blockStarts = [int(bs) for bs in r['blockStarts'].replace(',','\t').strip().split('\t')]
        lastExonKey = ''
        for bsi, bst in zip(blockSizes, blockStarts):
            tempStart = int(r['chromStart'])+bst
            tempEnd = tempStart + bsi
            tempKey = ':'.join([r['chrom'],str(tempStart),str(tempEnd),r['strand']])
            '''
            if tempKey in exons:
                exons[tempKey].append(':'.join([r['name'],str(tempCount)]))
            else:
                exons[tempKey] = [':'.join([r['name'],str(tempCount)])]
            '''
            if tempKey not in exons:
                exons[tempKey] = [tempKey,]
            if lastExonKey != "" and tempKey not in exons[lastExonKey]:
                exons[lastExonKey].append(tempKey)
            lastExonKey = tempKey
        if lastExonKey != '' and '<' not in exons[lastExonKey]:
            exons[lastExonKey].append('<')
    return exons

def writeToFile(exons, out):
    exonList = exons.items()
    exonList.sort(key=lambda k:(k[0].split(':')[0],k[0].split(':')[1],k[0].split(':')[2]))
    for exon in exonList:
        exonChrom,exonStart,exonEnd,exonStrand = exon[0].split(':')
        exonSize = str(int(exonEnd)-int(exonStart))
        exonName = ','.join(exon[1])
        #exonName += ','+exon[0]
        out.write(exonChrom+'\t'+exonStart+'\t'+exonEnd+'\t'+exonName+'\t0\t'+exonStrand+'\t'+exonStart+'\t'+exonEnd+'\t0\t1\t'+exonSize+'\t0\n')
    out.close()


if __name__=='__main__':
    records = readBed.BEDReader("../mm9.refFlat.bed")
    outfile = open("../mm9.refFlat.exons.new.bed",'w')
    result = removeRandom(records)
    exons = loadExons(result)
    writeToFile(exons, outfile)
