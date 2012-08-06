#!/usr/bin/python
import os,sys
sys.path.append('.')
import readBed
import readFlat
import platform
import re



'''
Generage dictionary of the junctions in a dict structure. The keys
are "chr_start_end" strings.
'''
def generateRefDict(reference):
    refDict = dict()
    count = 0
    if reference.endswith('txt'):
        ref = readFlat.FlatReader(reference)
        print 'ref'
        for index,feat in enumerate(ref):
            if not re.match(r'^chr((\d+)|X|Y)$',feat['chrom']):
                continue
            count+=1
            blockCount = int(feat['exonCount'])
            starts = feat['exonStarts'].replace(',','\t').strip().split('\t')
            ends = feat['exonEnds'].replace(',','\t').strip().split('\t')
            for i in range(blockCount-1):
                tempstart = int(ends[i])
                tempend = int(starts[i+1])
                tempkey = feat['chrom']+'_'+str(tempstart)+'_'+str(tempend)
                if tempkey in refDict:
                    refDict[tempkey].add(index)
                else:
                    refDict[tempkey] = set([index,])
            if count%10000 == 0:
                print count, ' records processed'

    else:
        ref = readBed.BEDReader(reference)
        print 'bed'
        for index,feat in enumerate(ref):
            #feat = line.split('\t')
            if not re.match(r'chr((\d+)|X|Y)$',feat['chrom']):
                continue
            count += 1
            blockCount = int(feat['blockCount'])
            starts = feat['blockStarts'].replace(',','\t').strip().split('\t')
            sizes = feat['blockSizes'].replace(',','\t').strip().split('\t')
            for i in range(blockCount - 1):
                tempstart = int(feat['chromStart']) + int(starts[i]) + int(sizes[i])
                tempend = int(feat['chromStart'])+int(starts[i+1])
                tempkey = feat['chrom']+"_"+str(tempstart)+'_'+str(tempend)
                if tempkey in refDict:
                    refDict[tempkey].add(index)
                    #refDict[tempkey].append(feat)
                else:
                    refDict[tempkey] = set([index,])
                    #refDict[tempkey] = [feat,]
            if count%10000 == 0:
                print count, ' records processed'

    #try:
    #    feats = ref.read().strip().split('\n')
    #except IOError:
    #    exit('Empty reference file.')

    #print len(feats)
    print "Total of ", count, " records processed."
    #f = open(reference+'.index','w')
    #for record in refDict.items():
    #    f.write(record[0].split('_')[0] + '\t' + record[0].split('_')[1] +'\n')
    #f.close()
    print 'Finish processing dictionary'
    return refDict



'''
Find the novel junctions in the inputfile. If the entries in the
input can contain several junctions, then this will use the generateDict()
to generate a dictionary similar to the refDict.
'''
def findNovelJunctions(refDict,inputfiles,mismatch,multijunc):
    print 'Now generating junctions---------------------'
    for inputfile in inputfiles:
        infile = readBed.BEDReader(inputfile)
        currJuncs = []
        for r in infile:
            currJuncs.append(r)
        if multijunc:
            inputContent = currJuncs
            currJuncs = generateRefDict(inputfile)

        sameJunctions = []
        newJunctions = []
        misajunctions = []
        for feat in currJuncs:
            if not multijunc:
                tempstart = int(feat['chromStart'])+int(feat['blockStarts'].split(',')[0]) + int(feat['blockSizes'].split(',')[0])
                tempend = int(feat['chromStart'])+int(feat['blockStarts'].split(',')[1])
                chrom = feat['chrom']
            else:
                tempstart = int(feat.split('_')[1])
                tempend = int(feat.split('_')[2])
                chrom = feat.split('_')[0]
            tempkeys = []
            firstkey = chrom+'_'+str(tempstart)+'_'+str(tempend)
            tempkeys.append(firstkey)
            for i in range(-mismatch,mismatch+1):
                for j in range(-mismatch,mismatch+1):
                    if i == 0 and j == 0:
                        continue
                    else:
                        tempkeys.append(chrom+'_'+str(tempstart+i)+'_'+str(tempend+j))
            exists = False
            #Generate the keys that within the user allowed difference.
            for tempkey in tempkeys:
                if tempkey in refDict:
                    exists = True
                    if tempkey != firstkey:
                        misajunctions.append(feat)
                    break
            if exists:
                if multijunc:
                    juncI = currJuncs[feat].pop()
                    currj = inputContent[juncI]
                    juncName = currj['name']
                    tempJunc = {'chrom':chrom,'chromStart':tempstart-3,'chromEnd':tempend+3,'name':currj['name'],'score':currj['score'],'strand':currj['strand'],'thickStart':tempstart-3,'thickEnd':tempend+3,'itemRgb':currj['itemRgb'],'blockCount':2,'blockSizes':'3,3','blockStarts':'0,%s'%(tempend-tempstart+3,)}
                    sameJunctions.append(tempJunc)
                else:
                    sameJunctions.append(feat)
            else:
                if multijunc:
                    juncI = currJuncs[feat].pop()
                    currj = inputContent[juncI]
                    juncName = currj['name']
                    tempJunc = {'chrom':chrom,'chromStart':tempstart-3,'chromEnd':tempend+3,'name':currj['name'],'score':currj['score'],'strand':currj['strand'],'thickStart':tempstart-3,'thickEnd':tempend+3,'itemRgb':currj['itemRgb'],'blockCount':2,'blockSizes':'3,3','blockStarts':'0,%s'%(tempend-tempstart+3,)}
                    newJunctions.append(tempJunc)
                else:
                    newJunctions.append(feat)

        print len(newJunctions)+len(sameJunctions), " junctions in total detected"
        print len(newJunctions), " novel junctions detected"
        print len(sameJunctions), " old junctions detected"
        newwriter = readBed.BEDWriter(inputfile+'.new.bed')
        samewriter = readBed.BEDWriter(inputfile+'.same.bed')
        newwriter.writerows(newJunctions)
        samewriter.writerows(sameJunctions)
        if mismatch > 0:
            misawriter = readBed.BEDWriter(inputfile+'.misa.bed')
            misawriter.writerows(misajunctions)


if __name__=='__main__':
    if platform.python_version().startswith('2.7'):
        import argparse
        parser = argparse.ArgumentParser(description="Find the novel junctions.")
        parser.add_argument('reference',help='The reference file should be in bed format or refflat format')
        parser.add_argument('inputfiles',nargs='+',help='The input file that contains the junctions to compare with the reference')
        parser.add_argument('-m','--mismatch',dest='mismatch',type=int,default=0,help='The position difference that can be allowed for each anchor.')
        parser.add_argument('-l','--multijunc',dest='multijunc',action='store_true',default=False,help='If each entry contains several junctions, this should be set to be true.')
        try:
            options = parser.parse_args()
        except IOError, msg:
            parser.error(str(msg))
        refDict = generateRefDict(options.reference)
        findNovelJunctions(refDict,options.inputfiles,options.mismatch,options.multijunc)
    else:
        from optparse import OptionParser
        parser = OptionParser("usage: %prog [options] reference query1 [query2 query3 ...]")
        #parser.add_option('-r','--reference',dest='reference',help='The reference file should be in bed format or refflat format')
        #parser.add_option('-i','--input',dest='inputfile',help='The input file that contains the junctions to compare with the reference')
        parser.add_option('-m','--mismatch',dest='mismatch',type=int,default=0,help='The position difference that can be allowed for each anchor.')
        parser.add_option('-l','--multijunc',dest='multijunc',action='store_true',default=False,help='If each entry contains several junction.')
        (options,args) = parser.parse_args()
        '''
        if options.reference is None:
            print 'The reference file is missing\n'
            parser.print_help()
            exit(-1)
        if options.inputfile is None:
            print 'The input file is missing\n'
            parser.print_help()
            exit(-1)
        '''
        if len(args) < 2:
            print "Tow few arguments.\n"
            parser.print_help()
            exit(-1)
        refDict = generateRefDict(args[0])
        findNovelJunctions(refDict,args[1:],options.mismatch,options.multijunc)

