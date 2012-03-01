import os,sys
inputdir = '../splicemap/cat/output'
outdir = '../rawbin/all_fa'
multisplitdir = '../rawbin/all_multireads' # the directory to store the multi-split reads
bamdir = '../rawbin/all_bam'
nonehit = '../rawbin/all_nonehit'

indir = os.walk(inputdir).next()

for filename in indir[2]:
    if filename.endswith('.bed'):
        #and not os.path.exists(os.path.join(outdir,filename+'.splitreads.bam')):
        print filename
        path = os.path.join(indir[0],filename)
        cmd = 'python extractSplitReads.py -o ' + outdir + ' '
        #os.system('python extractSplitReads.py -o ' + outdir + ' ' +path+' ../splicemap/good_hits.bam')
        if multisplitdir != '':
            cmd += '-m ' + multisplitdir + ' '
        if bamdir != '':
            cmd += '-b ' + bamdir + ' '
        if nonehit != '':
            cmd += '-n ' + nonehit + ' '
        cmd += path +' ../splicemap/good_hits.bam'
        #print cmd
        os.system(cmd)

