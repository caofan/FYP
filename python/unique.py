import os,sys
import csv

def within(pos1,pos2,mis):
    if pos1-mis<=pos2 and pos1+mis>=pos2:
        return True
    return False

if __name__=='__main__':
    base = 'rsup'
    infile = '../'+base+'.bed'
    #allows +/- 10bp misalignment
    diff = 10
    reader = csv.reader(open(infile,'rb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    juncs = []
    unique = []
    nonunique=[]
    for row in reader:
        juncs.append(row)
    for junc in juncs:
        count = 0
        for row in juncs:
            if junc[0] == row[0]:
                if within(int(row[1]),int(junc[1]),diff) and within(int(row[2]),int(junc[2]),diff):
                    count += 1
        if count > 1:
            nonunique.append(junc)
        else:
            unique.append(junc)

    writer = csv.writer(open('../'+base+'.uni.bed','wb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    writer.writerows(unique)
    writer = csv.writer(open('../'+base+'.nun.bed','wb'),delimiter='\t',quoting=csv.QUOTE_NONE)
    writer.writerows(nonunique)
