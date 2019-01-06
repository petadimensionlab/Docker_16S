import os, re, glob, shutil ,sys
from Bio import SeqIO
args = sys.argv

def append(input_fn,output_handle,fmt='fastq',sep='.',sample_name=None):
    ## Taken from MICCA append function ##
    def append_sample_name(records,sample_name):
        for record in records:
            record.id = '{0};sample={1}'.format(record.id,sample_name)
            try:
                description = record.description.split(None,1)[1]
            except IndexError:
                description = ''
            record.description = description
            yield record
    if sample_name is None:
        sample_name = os.path.basename(input_fn).split(sep)[0]
        sample_name = sample_name.replace('_merged','')
    sample_name_nows = re.sub('\s+','_',sample_name)
    records_in = SeqIO.parse(input_fn,fmt)
    records_out = append_sample_name(records_in,sample_name_nows)
    SeqIO.write(records_out,output_handle,fmt)

#### configurations ####
is_pairend = args[2]
is_filecheck = 'yes'
is_fastqc = 'yes'
is_trim = 'yes'
is_append = 'yes'
is_qualfilter = 'yes'
is_unique = 'yes'
is_sort = 'yes'
is_chimeracheck = 'yes'
is_uchime_ref = 'yes'
is_cluster = 'yes'
is_taxonomy = 'yes'
is_denovopicking = 'yes'
is_mapping = 'yes'
is_deletefiles = 'yes'

## softwares ##
vsearch = '/opt/conda/pkgs/vsearch-2.7.0-1/bin/vsearch'
fastqc = '/opt/conda/pkgs/fastqc-0.11.7-4/bin/fastqc'
Trimmomatic = '/opt/conda/pkgs/trimmomatic-0.36-5/share/trimmomatic-0.36-5/trimmomatic.jar'
mothur = '/opt/conda/bin/mothur'
## input / output file ##
input_dir = '/condir/'+args[1]
output_dir = input_dir
config_file_name = 'config.csv'

## reference file ##
chimera_ref = '/ref/chimeras/rRNA16S.gold.NAST_ALIGNED.fasta'
# GreenGenes #
gg97fasta = '/ref/gg_13_5_otus/rep_set/97_otus.fasta'
gg97tax = '/ref/gg_13_5_otus/taxonomy/97_otu_taxonomy.txt'
# Silva #
silva_fasta = '/ref/silva_128/rep_set/97_otus_16S.fasta'
silva_tax = '/ref/silva_128/taxonomy/consensus_taxonomy_all_levels.txt'
# RDP #
rdp_fasta = '/ref/trainset16_022016.rdp/trainset16_022016.rdp.fasta'
rdp_tax = '/ref/trainset16_022016.rdp/trainset16_022016.rdp.tax'

fastalist = { 'gg':gg97fasta,'silva':silva_fasta,'rdp':rdp_fasta }
taxlist = { 'gg':gg97tax,'silva':silva_tax,'rdp':rdp_tax }
header = { 'gg':'gg','silva':'silva','rdp':'rdp' }

reflist = ['gg','silva','rdp']
ref = reflist[0]

## parameters ##
maxdiffs = 5 #10
minmergelen = 1
maxmergelen = 600
minseqlength = 50
maxee_rate = 0.5/100 #%
minsize = 2
thr = 0.97
threadsnum = 4

## file names ##
appended_file_name = 'appended.fastq'
filtered_file_name = 'filtered.fasta'
unique_file_name = 'unique.fasta'
sorted_file_name = 'sorted.fasta'
nonchimera_file_name = 'nonchimera.fasta'
centroids_file_name = 'centroids.fasta'
denovootu_file_name = header[ref]+'denovootu.txt'
hits_file_name = header[ref]+'hits.txt'
refotuname_file_name = header[ref]+'otu_name.txt'
refotuID_file_name = header[ref]+'otu.txt'

## file integrity check ##
if is_filecheck=='yes':
    if is_pairend=='yes':
        fwds = set(glob.glob('^*_1.fastq.gz'))
        revs = set(glob.glob('^*_2.fastq.gz'))
        ffiles = rfiles = set()
        config_file = os.path.join(input_dir,config_file_name)
        fr = open(config_file,'r').readlines()
        for line in fr:
            line = line.replace('\n','')
            tmp = line.split(',')
            sID = tmp[0]
            ffiles.add(tmp[1])
            rfiles.add(tmp[2])
        fwddiff = list(fwds-ffiles)
        revdiff = list(revs-rfiles)
        if len(fwddiff)>0 or len(revdiff)>0:
            msg = 'Some FASTQ files are lacking. Quit computation.'
            print( msg )
            quit()
        else:
            msg = 'All FASTQ files exist. Continue computation.'
            print( msg )
    else:
        fwds = set(glob.glob('*.fastq.gz'))
        ffiles = set()
        config_file = os.path.join(input_dir,config_file_name)
        fr = open(config_file,'r').readlines()
        for line in fr:
            line = line.replace('\n','')
            tmp = line.split(',')
            sID = tmp[0]
            ffiles.add(tmp[1])
        fwddiff = list(fwds-ffiles)
        if len(fwddiff)>0:
            msg = 'Some FASTQ files are lacking. Quit computation.'
            print( msg )
            quit()
        else:
            msg = 'All FASTQ files exist. Continue computation.'
            print( msg )

## Quality Check by fastqc ##
if is_fastqc=='yes':
    fastqc_dir = os.path.join(input_dir,'fastqc')
    if not os.path.exists(fastqc_dir):
        os.mkdir(fastqc_dir)
    if is_pairend=='yes':
        config_file = os.path.join(input_dir,config_file_name)
        fr = open(config_file,'r').readlines()
        for line in fr:
            line = line.replace('\n','')
            tmp = line.split(',')
            sID = tmp[0]
            fwdfile = os.path.join(input_dir,tmp[1])
            revfile = os.path.join(input_dir,tmp[2])
            cmd = '%s %s' % (fastqc,fwdfile)
            os.system(cmd)
            cmd = '%s %s' % (fastqc,revfile)
            os.system(cmd)
    else:
        config_file = os.path.join(input_dir,config_file_name)
        fr = open(config_file,'r').readlines()
        for line in fr:
            line = line.replace('\n','')
            tmp = line.split(',')
            sID = tmp[0]
            fn = tmp[1]
            cmd = '%s %s' % (fastqc,fn)
            os.system(cmd)
    qcfiles = glob.glob('*_fastqc*')
    for qcfile in qcfiles:
        before = os.path.join(input_dir,qcfile)
        after = os.path.join(fastqc_dir,qcfile)
        shutil.move(before,after)

## merge pair-end files ##
if is_pairend=='yes':
    input_fns = []
    config_file = os.path.join(input_dir,config_file_name)
    fr = open(config_file,'r').readlines()
    for line in fr:
        line = line.replace('\n','')
        tmp = line.split(',')
        sID = tmp[0]
        fwdfile = os.path.join(input_dir,tmp[1])
        revfile = os.path.join(input_dir,tmp[2])
        paired_fwdfile = os.path.join(input_dir,'paired_'+tmp[1])
        paired_revfile = os.path.join(input_dir,'paired_'+tmp[2])
        unpaired_fwdfile = os.path.join(input_dir,'unpaired_'+tmp[1])
        unpaired_revfile = os.path.join(input_dir,'unpaired_'+tmp[2])
        tmpfile = os.path.join(output_dir,sID+'_merged.fasta')
        input_fns.append(tmpfile)
        ## trimming by Trimmomatic ##
        if is_trim=='yes':
            cmd = 'java -jar %s PE %s %s %s %s %s %s SLIDINGWINDOW:4:15 MINLEN:%d' % (Trimmomatic,fwdfile,revfile,paired_fwdfile,unpaired_fwdfile,paired_revfile,unpaired_revfile,minseqlength)
            os.system(cmd)
            cmd = '%s -fastq_mergepairs %s --reverse %s -fastq_maxdiffs %s -fastq_minmergelen %s -fastq_maxmergelen %s -fastqout %s' % (vsearch,paired_fwdfile,paired_revfile,str(maxdiffs),str(minmergelen),str(maxmergelen),tmpfile)
            os.system(cmd)
        else:
            cmd = '%s -fastq_mergepairs %s --reverse %s -fastq_maxdiffs %s -fastq_minmergelen %s -fastq_maxmergelen %s -fastqout %s' % (vsearch,fwdfile,revfile,str(maxdiffs),str(minmergelen),str(maxmergelen),tmpfile)
            os.system(cmd)
else:
    input_fns = []
    config_file = os.path.join(input_dir,config_file_name)
    fr = open(config_file,'r').readlines()
    for line in fr:
        line = line.replace('\n','')
        tmp = line.split(',')
        sID = tmp[0]
        fn = tmp[1]
        #input_file = os.path.join(output_dir,sID+'.fastq.gz')
        input_file = os.path.join(output_dir,fn)
        ## trimming by Trimmomatic ##
        if is_trim=='yes':
            tmpfile = os.path.join(output_dir,sID+'_trimmed.fastq')
            cmd = 'java -jar %s SE %s %s SLIDINGWINDOW:4:15 MINLEN:%d' % (Trimmomatic,input_file,tmpfile,minseqlength)
            os.system(cmd)
            input_fns.append(tmpfile)
        else:
            cmd = 'gunzip %s' % (input_file)
            #os.system(cmd)
            #input_fn = input_file.replace('.gz','')
            #input_fns.append(input_fn)

## append file ##
appended_file = os.path.join(output_dir,appended_file_name)
if is_append=='yes':
    if is_pairend=='yes':
        ## append merged-pair fastq files ##
        with open(appended_file,'w') as output_handle:
            for input_fn in input_fns:
                append(input_fn,output_handle,'fastq','.')
        output_handle.close()
        ## remove tmporary files ##
        tmpfastqfiles = glob.glob('*_merged.fastq')
        for item in tmpfastqfiles:
            os.remove(item)
    else:
        ## append fastq files ##
        with open(appended_file,'w') as output_handle:
            for input_fn in input_fns:
                append(input_fn,output_handle,'fastq','.')
        output_handle.close()
        if is_trim!='yes':
            for input_fn in input_fns:
                cmd = 'gzip -f %s' % (input_fn)
                os.system(cmd)
else:
    msg = 'Skip appending.'
    print( msg )

## filtering ##
filtered_file = os.path.join(output_dir,filtered_file_name)
if is_qualfilter=='yes':
    cmd = '%s --fastq_filter %s --fastaout %s --fastq_maxee_rate %s' % (vsearch,appended_file,filtered_file,str(maxee_rate))
    os.system(cmd)
else:
    msg = 'Skip filtering.'
    print( msg )

## get unique sequences ##
unique_file = os.path.join(output_dir,unique_file_name)
if is_unique=='yes':
    cmd = '%s --derep_fulllength %s --output %s --sizeout --minseqlength %s' % (vsearch,filtered_file,unique_file,str(minseqlength))
    os.system(cmd)
else:
    msg = 'Skip unifying.'
    print( msg )

## sorted sequences ##
sorted_file = os.path.join(output_dir,sorted_file_name)
if is_sort=='yes':
    cmd = '%s --sortbysize %s --output %s --minsize %s' % (vsearch,unique_file,sorted_file,str(minsize))
    os.system(cmd)
    #print(cmd)
else:
    msg = 'Skip sorting.'
    print( msg )

## chimera detection ##
nonchimera_file = os.path.join(output_dir,nonchimera_file_name)
if is_chimeracheck=='yes':
    if is_uchime_ref=='yes':
        cmd = '%s \"#chimera.uchime(fasta=%s,reference=%s,processors=%s); remove.seqs(fasta=current,accnos=current);\"' % (mothur,sorted_file,chimera_ref,str(threadsnum))
        #print(cmd)
        os.system(cmd)
        tmp_name = sorted_file_name.replace('.fasta','')
        tmp_name = tmp_name+'.pick.fasta'
        tmpfile = os.path.join(output_dir,tmp_name)
        cmd = 'mv %s %s' % (tmpfile,nonchimera_file)
        os.system(cmd)
        #print(cmd)
    else:
        cmd = '%s --uchime_denovo %s --nonchimeras %s' % (vsearch,sorted_file,nonchimera_file)
        os.system(cmd)
else:
    msg = 'Skip chimera check.'
    print( msg )

## cluster ##
centroids_file = os.path.join(output_dir,centroids_file_name)
if is_cluster=='yes':
    #cmd = '%s --cluster_fast %s --id %s --centroids %s' % (vsearch,sorted_file,str(thr),centroids_file)
    cmd = '%s --cluster_smallmem %s --id %s --consout %s --usersort' % (vsearch,nonchimera_file,str(thr),centroids_file)
    os.system(cmd)
else:
    msg = 'Skip clustering.'
    print( msg )

## Assign taxonomy ##
hits_file = os.path.join(output_dir,hits_file_name)
if is_taxonomy=='yes':
    cmd = '%s --usearch_global %s --db %s --strand both --id %s --threads %s --userfields \"query+target+id\" --userout %s' % (vsearch,centroids_file,fastalist[ref],str(thr),str(threadsnum),hits_file)
    os.system(cmd)
else:
    msg = 'Skip taxonomy assignment.'
    print( msg )

## de novo picking ##
denovootu_file = os.path.join(output_dir,denovootu_file_name)
if is_denovopicking=='yes':
    cmd = '%s --usearch_global %s --db %s --strand both --id %s --threads %s --otutabout %s' % (vsearch,filtered_file,centroids_file,str(thr),str(threadsnum),denovootu_file)
    os.system(cmd)
else:
    msg = 'Skip denovo picking.'
    print( msg )

## Mapping to GreenGene database ##
if is_mapping=='yes':
    refid2refname = {}
    frref = open(taxlist[ref],'r').readlines()
    for line in frref:
        line = line.replace('\n','')
        tmp = line.split('\t')
        refid2refname[tmp[0]] = tmp[1]
    ## make a mapping from "user defined sequence ID" to "greengene ID"
    myid2refid = {}
    myid2refname = {}
    frh = open(hits_file,'r').readlines()
    for line in frh:
        tmp = line.split('\t') # refID: tmp[1]
        refid = tmp[1]
        refname = refid2refname[refid]
        myids = tmp[0].split(';') # myID : tmp[0]
        myid2refid[myids[0]] = refid
        myid2refname[myids[0]] = refname

    frd = open(denovootu_file,'r').readlines()
    header = frd[0] # restore header information
    frd.pop(0) # remove the header line
    ## taxonomy name version ##
    refotu_name_file = os.path.join(output_dir,refotuname_file_name)
    fw = open(refotu_name_file,'w')
    fw.write(header)
    for line in frd:
        tmp = line.split('\t') # tmp[0]: myID
        try:
            tmp[0] = myid2refname[tmp[0]]
            fw.write('\t'.join(tmp))
        except:
            continue
    fw.close()
    ## taxonomy ID version ##
    refotu_ID_file = os.path.join(output_dir,refotuID_file_name)
    fw = open(refotu_ID_file,'w')
    fw.write(header)
    for line in frd:
        tmp = line.split('\t') # tmp[0]: myID
        try:
            tmp[0] = myid2refid[tmp[0]]
            fw.write('\t'.join(tmp))
        except:
            continue
    fw.close()
else:
    msg = 'Skip taxonomy mapping.'
    print( msg )

## delete intermediate files ##
if is_deletefiles=='yes':
    os.chdir(input_dir)
    mothurlog_file = glob.glob('mothur*logfile')
    mothurpick_file = glob.glob('sorted.pick.fasta')
    paired_files = glob.glob('paired*.fastq.gz')
    unpaired_files = glob.glob('unpaired*.fastq.gz')
    merged_files = glob.glob('*_merged.fasta')
    trimmed_files = glob.glob('*_trimmed.fasta')
    uchime_files = glob.glob('*uchime*')
    rm_files = ['appended.fastq','unique.fasta','nonchimera.fasta']+mothurlog_file+uchime_files+mothurpick_file
    if is_pairend=='yes':
        rm_files = rm_files+paired_files+unpaired_files+merged_files
    else:
        rm_files = rm_files+trimmed_files
    for file in rm_files:
        cmd = 'rm %s' % (file)
        os.system(cmd)

