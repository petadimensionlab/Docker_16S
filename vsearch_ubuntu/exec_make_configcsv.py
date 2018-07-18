import os , sys

args = sys.argv

data_dir = '/condir/'+args[1]
is_pairend = args[2]

os.chdir(data_dir)
fr = open('SRR_Acc_List.txt','r').readlines()
fw = open('config.csv','w')

if is_pairend=='yes':
    R1 = '_1.fastq.gz'
    R2 = '_2.fastq.gz'
    for line in fr:
        Acc = line.replace('\n','')
        R1name = Acc+R1
        R2name = Acc+R2
        lst = [Acc,R1name,R2name]
        saveline = ','.join(lst)
        fw.write(saveline+'\n')
    #fw.write('\n')
    fw.close()
else:
    R = '.fastq.gz'
    for line in fr:
        Acc = line.replace('\n','')
        Rname = Acc+R
        lst = [Acc,Rname]
        saveline = ','.join(lst)
        fw.write(saveline+'\n')
    #fw.write('\n')
    fw.close()