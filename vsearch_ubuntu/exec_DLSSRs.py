# http://blog.amelieff.jp/?eid=231191
import os, shutil, glob , sys

args = sys.argv

cdir = os.getcwd()
data_dir = '/condir/'+args[1]
#prefetch_dir = '/Users/snakaoka/ws/apps/sratoolkit/bin/'
prefetch_dir = '/sratool/bin/'
DL_dir = '/root/ncbi/public/sra/'

## convert sra to fastq ##
os.chdir(data_dir)
cmd = prefetch_dir+'./prefetch --option-file '+os.path.join(data_dir,'SRR_Acc_List.txt')
os.system(cmd)
os.chdir(cdir)

## move all data to an assigned folder ##
os.chdir(DL_dir)
files = glob.glob('*.sra')
for f in files:
    tmp = os.path.join(data_dir,f)
    shutil.move(f,tmp)