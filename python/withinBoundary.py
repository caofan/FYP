import os,sys
import readBed

if __name__=='__main__':
    reffile = '../mm9.genomeBoundaries.txt'
    juncfile = '../outputs/rawbin.splicemap.4/stats/rawbin.est.new.bed.large'
    outputdir = '../outputs/rawbin.splicemap.4/large'
    genes = []
    f = open(reffile)
    lines = f.read().strip().split('\n')
    for line in lines:
        tokens = line.split('\t')
        genes.append([tokens[0],tokens[1],int(tokens[2]),int(tokens[3])])
    f.close()
    reader=readBed.BEDReader(juncfile)
    within = []
    halfin = []
    nonein = []
    diffgene = []
    for line in reader:
        #This specifies which gene they are located. If it is -1, then
        #it is in the intergenic region.
        leftin = -1
        rightin = -1
        starts = line['blockStarts'].replace(',',' ').strip().split(' ')
        sizes = line['blockSizes'].replace(',',' ').strip().split(' ')
        juncStart = int(line['chromStart'])+int(starts[0])+int(sizes[0])
        juncEnd = int(line['chromStart'])+int(starts[1])
        for i,gene in enumerate(genes):
            if line['chrom'] == gene[0] and line['strand']==gene[1] and juncStart >= gene[2] and juncStart <= gene[3]:
                if leftin != -1:
                    print "left ",leftin
                leftin = i
            if line['chrom'] == gene[0] and line['strand']==gene[1] and juncEnd >= gene[2] and juncEnd <= gene[3]:
                if rightin != -1:
                    print "right ",rightin
                rightin = i
        if leftin != -1 and leftin == rightin:
            within.append(line)
        elif (leftin!= -1 and rightin == -1) or (rightin != -1 and leftin == -1):
            halfin.append(line)
        elif leftin != -1 and rightin != -1 and leftin != rightin:
            diffgene.append(line)
        else:
            nonein.append(line)

    print "within ", len(within)
    print "halfin ", len(halfin)
    print "nonein ", len(nonein)
    print "diffgene ", len(diffgene)
    writer  = readBed.BEDWriter(os.path.join(outputdir,'within.bed'))
    writer.writerows(within)
    writer  = readBed.BEDWriter(os.path.join(outputdir,'halfin.bed'))
    writer.writerows(halfin)
    writer  = readBed.BEDWriter(os.path.join(outputdir,'nonein.bed'))
    writer.writerows(nonein)
    writer  = readBed.BEDWriter(os.path.join(outputdir,'diffgene.bed'))
    writer.writerows(diffgene)

