import os,sys
inputdir = '../splicemap/cat/output'
outdir = '../splicemap/cat/output'

indir = os.walk(inputdir).next()

for filename in indir[2]:
    path = os.path.join(indir[0],filename)
    os.system('python findNovelJunctions.py -r ../outputs/cufflink.400000/rawbin.4.noann/junctions.bed -i ' +path+' -m 0')
