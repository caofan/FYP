import os,sys
import re

readLength = 76
rootdir = os.walk('.')
rootdir.next()
for indir in rootdir:
    print indir[0]
    f = open(os.path.join(indir[0],"Mapped.sam"))
    outfile = os.path.join(indir[0],"Mapped_modified.sam")
    out = open(outfile,'w')
    out.write('''@HD\tVN:1.0\tSO:coordinate\n@SQ\tSN:chr1\tLN:197195432\n@SQ\tSN:chr10\tLN:129993255\n@SQ\tSN:chr11\tLN:121843856\n@SQ\tSN:chr12\tLN:121257530\n@SQ\tSN:chr13\tLN:120284312\n@SQ\tSN:chr14\tLN:125194864\n@SQ\tSN:chr15\tLN:103494974\n@SQ\tSN:chr16\tLN:98319150\n@SQ\tSN:chr17\tLN:95272651\n@SQ\tSN:chr18\tLN:90772031\n@SQ\tSN:chr19\tLN:61342430\n@SQ\tSN:chr2\tLN:181748087\n@SQ\tSN:chr3\tLN:159599783\n@SQ\tSN:chr4\tLN:155630120\n@SQ\tSN:chr5\tLN:152537259\n@SQ\tSN:chr6\tLN:149517037\n@SQ\tSN:chr7\tLN:152524553\n@SQ\tSN:chr8\tLN:131738871\n@SQ\tSN:chr9\tLN:124076172\n@SQ\tSN:chrM\tLN:16299\n@SQ\tSN:chrX\tLN:166650296\n@SQ\tSN:chrY\tLN:15902555\n''')
    for r in f:
        tokens = r.strip().split("\t")
        total = 0
	readLength = int(tokens[-2].split(":")[-1])
        if re.match(r"^\d+M$",tokens[5]):
            match = re.match(r"^(?P<m>\d+)M$",tokens[5])
            if int(match.group("m")) < readLength:
                remain = readLength-int(match.group('m'))
                tokens[5] += str(remain)+"X"
            #tokens[10] = int(match.group("m"))*"I"
            tokens[10] = readLength*"I"
            total = readLength
        elif re.match(r"^(?P<m1>\d+)M(?P<n>\d+)N(?P<m2>\d+)M",tokens[5]):
	    tokens[10] = readLength*"I"
            match=re.match(r"^(?P<m1>\d+)M(?P<n>\d+)N(?P<m2>\d+)M",tokens[5])
            total = int(match.group("m1"))+int(match.group("m2"))
            if total > readLength:
                remain = readLength - int(match.group('m1'))
                tokens[5] = match.group("m1")+"M"+str(int(match.group("n"))+total-readLength)+"N"+str(remain)+"M"
                tokens.append("XO:i:%d"%(total-readLength,))
        if total < readLength:
            continue
        outline = ""
        for token in tokens:
            outline += token + "\t"
        outline = outline.strip()
        out.write(outline+"\n")
    out.close()
    f.close()
    bamfile = os.path.join(indir[0],"Mapped_modified.bam")
    os.system("samtools view -Sb %s -o %s"%(outfile,bamfile,))
    os.system("samtools sort %s %s"%(bamfile,os.path.join(indir[0],"Mapped_modified",)))
    os.system("samtools index %s"%(bamfile,))

