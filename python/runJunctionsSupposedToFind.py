import os,sys
import re

inputdir = '../mm9-chromosome-single/polyAfiltered/'
refgenefile = '../mm9.refFlat.bed'

indir = os.walk(inputdir).next()

for filename in indir[2]:
    match =  re.match(r'^chr(?P<id>\d+|X|Y).bed$',filename)
    if  match:
        print "processing ", filename
        infilePath = os.path.join(indir[0],filename)
        os.system("python findNovelJunctions.py -r %s -i %s -l"%(refgenefile,infilePath,))
        os.system("mv %s ../supposedJuncs/chr%s_supposed.bed"%(infilePath+".same.bed",match.group("id"),))
        os.system("rm %s"%(infilePath+".new.bed"))
