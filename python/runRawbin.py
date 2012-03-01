#!/usr/bin/python
import os

inputdir = "/data1/caofan/inputs/simulated"
indir = os.walk(inputdir).next()[2]

outdir = "/data1/caofan/outputs/rawbin/simulated"
if not os.path.exists(outdir):
    os.system("mkdir "+outdir)

for filename in indir:
    print filename
    inputfile= os.path.join(inputdir,filename)
    os.system("./rawbin/src/rnaseq -g mm9/mm9 -q "+inputfile + " -G 400000 -N 2")
    outputdir = os.path.join(outdir, filename)
    if not os.path.exists(outputdir):
        os.system("mkdir "+outputdir)
    print "Copying ", filename
    os.system("cp Mapped.sam coverage.bed junctions.bed " + outputdir)
    print "Finished ", filename
