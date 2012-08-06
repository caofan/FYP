#!/usr/bin/python
import os

inputdir = "/data1/caofan/inputs/simulated"
indir = os.walk(inputdir).next()[2]

cfgfile = "/data1/caofan/FYP/SpliceMap3352_example_linux-64/run.cfg"

for filename in indir:
    lines = []
    reader = open(cfgfile)
    isreadFile = False
    for l in reader:
        if isreadFile:
            lines.append(os.path.join(inputdir,filename)+"\n")
            isreadFile = False
        elif l.startswith("> reads_list1"):
            lines.append(l)
            isreadFile = True
        elif l.startswith("out_path"):
            outpath = "/data1/caofan/outputs/splicemap/simulated/"+filename
            lines.append("out_path = "+outpath+"\n")
            if ! os.path.exists(outpath):
                os.system("mkdir " + outpath)
        elif l.startswith("temp_path"):
            temppath = "/data1/caofan/outputs/splicemap/simulated/"+filename+"/tmp"
            lines.append("temp_path = " + temppath +"\n")
            if ! os.path.exists(temppath):
                os.system("mkdir " + temppath)
        else:
            lines.append(l)
    reader.close()
    writer = open(cfgfile,'w')
    for l in lines:
        writer.write(l)
    writer.close()
    os.system("./FYP/SpliceMap3352_example_linux-64/bin/runSpliceMap "+cfgfile)
