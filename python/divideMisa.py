#!/usr/bin/env python
import os,sys
sys.path.append('.')
import re
import readBed

comp = {}
comp['A'] = 'T'
comp['C'] = 'G'
comp['G'] = 'C'
comp['T'] = 'A'

def rc(signal):
    sigl = list(signal)
    sigl.reverse()
    rsig = ''.join(sigl)
    rcsig= []
    for r in rsig:
        try:
            rcsig.append(comp[r.upper()])
        except KeyError:
            rcsig.append(r)
    return ''.join(rcsig)

if __name__=='__main__':
    infile = '../outputs/cufflink.400000/rawbin.4.noann/cat/rsup.nun.bed.same.bed'
    reader = readBed.BEDReader(infile)
    junc = {}
    for r in reader:
        match = re.match(r'JUN.+\_(?P<signal>[A-Z]{2}-[A-Z]{2})\_.+',r['name'])
        if match:
            sig = match.group('signal')
            if sig in junc:
                junc[sig].append(r)
            elif rc(sig) in junc:
                junc[rc(sig)].append(r)
            else:
                junc[sig] = [r,]
    for key in junc:
        outfile = infile+'.'+key+'.bed'
        writer = readBed.BEDWriter(outfile)
        writer.writerows(junc[key])
