import os,sys
import readFlat

if __name__=='__main__':
    inputfile = '../mm9.refFlat.txt'
    output = '../mm9.genomeBoundaries.txt'
    reader = readFlat.FlatReader(inputfile)
    transcripts = []
    for line in reader:
        if len(line) <= 1:
            break
        if line['chrom'].endswith('random'):
            continue
        starts = []
        ends = []
        for token in line['exonStarts'].replace(',',' ').strip().split(' '):
            starts.append(int(token))
        for token in line['exonEnds'].replace(',',' ').strip().split(' '):
            ends.append(int(token))
        line['exonStarts'] = starts
        line['exonEnds'] = ends
        line['txStart'] = int(line['txStart'])
        line['txEnd'] = int(line['txEnd'])
        line['cdsStart'] = int(line['cdsStart'])
        line['cdsEnd'] = int(line['cdsEnd'])
        line['exonCount'] = int(line['exonCount'])
        transcripts.append(line)
    transcripts.sort(key=lambda k:(k['chrom'],k['txStart'],k['txEnd']))
    genes = []
    currentLeft = transcripts[0]['txStart']
    currentRight = transcripts[0]['txEnd']
    currentStrand = transcripts[0]['strand']
    currentChrom = transcripts[0]['chrom']
    currentGenename = transcripts[0]['geneName']
    for i in range(1,len(transcripts)):
        tx = transcripts[i]
        if tx['chrom'] == currentChrom and tx['strand'] == currentStrand and tx['txStart']<currentRight:
            if currentRight < tx['txEnd']:
                currentRight = tx['txEnd']
        else:
            gene = {}
            gene['start'] = currentLeft
            gene['end'] = currentRight
            gene['strand'] = currentStrand
            gene['chrom'] = currentChrom
            genes.append(gene)
            currentChrom = tx['chrom']
            currentLeft = tx['txStart']
            currentRight = tx['txEnd']
            currentStrand = tx['strand']
            currentGenename = tx['geneName']
    f = open(output,'w')
    genes.sort(key=lambda k:(k['chrom'],k['start'],k['end']))
    for gene in genes:
        f.write(gene['chrom'] + '\t' + gene['strand'] + '\t' + str(gene['start']) + '\t' + str(gene['end']) + '\n')
    f.close()
