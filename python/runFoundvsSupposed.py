#!/usr/bin/python26
import os,sys
import re
i=600
rawbinDir = "/home/Ken/caofan/outputs/rawbin/simulated-17OCT-trans"
spliceMapDir = "/home/Ken/caofan/outputs/splicemap/simulated"
tophatDir = '/home/Ken/caofan/outputs/tophat/simulated'

rawbinOutDir = "/home/Ken/caofan/supposedVSfound/rawbin/simulated-17OCT-trans"
spliceMapOutDir = "/home/Ken/caofan/supposedVSfound/splicemap/simulated"
tophatOutDir = "/home/Ken/caofan/supposedVSfound/tophat/simulated"
# the type of reference can be refgene or EST. Needs to be specified
# whenever the reference is changed.
referencedir = "/home/Ken/caofan/supposedJuncs"
ref_type = "supposed"

# For each chromosome, compare refgene/EST against it.
indirHandle = os.walk(rawbinDir)
# Skip the current directory
indirHandle.next()
print "---------------------------------------------Rawbin-------------------"
if not os.path.exists(rawbinOutDir):
    #os.mkdir(rawbinOutDir)
    os.system('mkdir -p '+ rawbinOutDir)
for currDir in indirHandle:
    currChrom = os.path.basename(currDir[0])
    chromMatch = re.match(r'^(?P<chrom>chr(\d+|X|Y)).',currChrom)
    print currChrom
    filename = os.path.join(currDir[0],"junctions.bed")
    reference = os.path.join(referencedir,chromMatch.group('chrom')+"_supposed.bed")
    os.system("python26 /home/Ken/caofan/python/findNovelJunctions.py -r %s -i %s -l"%(reference,filename))
    outfilePrefix = os.path.join(rawbinOutDir,ref_type+".rawbin."+currChrom)
    os.system("mv %s %s"%(filename+".same.bed",outfilePrefix+".same.bed"))
    os.system("mv %s %s"%(filename+".new.bed",outfilePrefix+".new.bed"))
    print "\n"

'''

print "\n-----------------------------------------SpliceMap-----------------"
indirHandle = os.walk(spliceMapDir)
#Skip the current directory
indirHandle.next()
for currDir in indirHandle:
    currChrom = os.path.basename(currDir[0])
    if not re.match(r"^chr(\d+|X|Y)\.fasta$",currChrom):
        continue
    currChrom = currChrom.split(".")[0]
    filename = os.path.join(currDir[0],"junction_color.bed")
    reference = os.path.join(referencedir,currChrom+"_supposed.bed")
    os.system("python26 findNovelJunctions.py -r %s -i %s -l"%(reference,filename))
    outfilePrefix = os.path.join(spliceMapOutDir,ref_type+".splicemap."+currChrom)
    os.system("mv %s %s"%(filename+".same.bed",outfilePrefix+".same.bed"))
    os.system("mv %s %s"%(filename+".new.bed",outfilePrefix+".new.bed"))
    print "\n"

print "\n-----------------------------------------Tophat---------------------"
indirHandle = os.walk(tophatDir)
#Skip the current directory
indirHandle.next()
if not os.path.exists(tophatOutDir):
    os.system("mkdir -p %s"%(tophatOutDir,))
for currDir in indirHandle:
    currChrom = os.path.basename(currDir[0])
    if not re.match(r"^chr(\d+|X|Y)\.fasta$",currChrom):
        continue
    currChrom = currChrom.split(".")[0]
    filename = os.path.join(currDir[0],"junctions.bed")
    reference = os.path.join(referencedir,currChrom+"_supposed.bed")
    os.system("python26 findNovelJunctions.py -r %s -i %s -l"%(reference,filename))
    outfilePrefix = os.path.join(tophatOutDir,ref_type+".tophat."+currChrom)
    os.system("mv %s %s"%(filename+".same.bed",outfilePrefix+".same.bed"))
    os.system("mv %s %s"%(filename+".new.bed",outfilePrefix+".new.bed"))
    print "\n"
'''
