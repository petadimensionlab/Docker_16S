 # required : sra-stat, fastq-dump #
import os, shutil, glob , sys

args = sys.argv

cdir = os.getcwd()
data_dir = '/condir/'+args[1]
pfastqdump_dir = '/sratool/bin/'

is_pairend = args[2]
is_gzip = 'yes'

thread_num = 8

## fastq dump ##
os.chdir(data_dir)
files = glob.glob('*.sra')
for f in files:
    if is_pairend=='yes':
        #cmd = 'export PATH=$PATH:'+pfastqdump_dir+' && '+pfastqdump_dir+'./pfastq-dump -t %d -split-3 -O %s %s' % (thread_num,data_dir,f)
        cmd = 'export PATH=$PATH:'+pfastqdump_dir+' && '+pfastqdump_dir+'./fastq-dump -split-files -O %s %s' % (data_dir,f)
        os.system(cmd)
    else:
        #cmd = 'export PATH=$PATH:'+pfastqdump_dir+' && '+pfastqdump_dir+'./pfastq-dump -t %d -O %s %s' % (thread_num,data_dir,f)
        cmd = 'export PATH=$PATH:'+pfastqdump_dir+' && '+pfastqdump_dir+'./fastq-dump -O %s %s' % (data_dir,f)
        os.system(cmd)

## gzip ##
if is_gzip=='yes':
    os.chdir(data_dir)
    files = glob.glob('*.fastq')
    for f in files:
        cmd = 'gzip %s' % (f)
        os.system(cmd)