#!/usr/bin/python
import os,sys

rawbinDir = "../outputs/rawbin/simulated"
spliceMapDir = "../outputs/splicemap/simulated"

rawbinOutDir = "../novelOut/rawbin/simulated"
spliceMapOutDir = "../novelOut/splicemap/simulated"

# the type of reference can be refgene or EST. Needs to be specified
# whenever the reference is changed.
reference = "../mm9.refFlat.bed"
ref_type = "refgene"

# For each chromosome, compare refgene/EST against it.
indirHandle = os.walk(rawbinDir)
# Skip the current directory
indirHandle.next()
for currDir in indirHandle:
    currChrom = os.path.basename(currDir[0])
    filename = os.path.join(currDir[0],"junctions.bed")
    os.system("python findNovelJunctions.py -i %s -r %s -l"%(reference,filename))
    outfilePrefix = os.path.join(rawbinOutDir,ref_type+".rawbin."+currChrom)
    os.system("mv %s %s"%(reference+".same.bed",outfilePrefix+".same.bed"))
    os.system("mv %s %s"%(reference+".new.bed",outfilePrefix+".new.bed"))


indirHandle = os.walk(spliceMapDir)
#Skip the current directory
indirHandle.next()
for currDir in indirHandle:
    currChrom = os.path.basename(currDir[0])
    filename = os.path.join(currDir[0],"junction_color.bed")
    os.system("python findNovelJunctions.py -i %s -r %s -l"%(reference,filename))
    outfilePrefix = os.path.join(spliceMapOutDir,ref_type+".splicemap."+currChrom)
    os.system("mv %s %s"%(reference+".same.bed",outfilePrefix+".same.bed"))
    os.system("mv %s %s"%(reference+".new.bed",outfilePrefix+".new.bed"))

