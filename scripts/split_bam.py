import os
import pysam

sample='testA'

allread=open('readname_'+sample+'.txt','r').read().split('\n')[:-1]
namesdel=set([c.split('\t')[0] for c in allread if 'DEL-small' in c])
nameldel=set([c.split('\t')[0] for c in allread if 'DEL-large' in c])
namesins=set([c.split('\t')[0] for c in allread if 'INS-small' in c])
namelins=set([c.split('\t')[0] for c in allread if 'INS-large' in c])
nameinv=set([c.split('\t')[0] for c in allread if 'INV' in c])
namedup=set([c.split('\t')[0] for c in allread if 'DUP' in c])
namehdr=set([c.split('\t')[0] for c in allread if 'HDR' in c])
namefull=set([c.split('\t')[0] for c in allread if 'Full-AAV-integration' in c])
nameunmodify=set([c.split('\t')[0] for c in allread if 'Unmodified' in c])
namesnp=set([c.split('\t')[0] for c in allread if 'SNP' in c])

f=pysam.AlignmentFile('filtered_'+sample+'.bam','rb')
f1=pysam.AlignmentFile(sample+'_split_unmodified.bam','wb',template=f)
f2=pysam.AlignmentFile(sample+'_split_indel_large_del.bam','wb',template=f)
f3=pysam.AlignmentFile(sample+'_split_indel_small_del.bam','wb',template=f)
f4=pysam.AlignmentFile(sample+'_split_indel_large_ins.bam','wb',template=f)
f5=pysam.AlignmentFile(sample+'_split_indel_small_ins.bam','wb',template=f)
f6=pysam.AlignmentFile(sample+'_split_indel_inv.bam','wb',template=f)
f7=pysam.AlignmentFile(sample+'_split_indel_dup.bam','wb',template=f)
f8=pysam.AlignmentFile(sample+'_split_snp.bam','wb',template=f)
f9=pysam.AlignmentFile(sample+'_split_aav_hdr.bam','wb',template=f)
f10=pysam.AlignmentFile(sample+'_split_aav_full.bam','wb',template=f)


allread=f.fetch()
for read in allread:
    if read.query_name in namesdel:
        f3.write(read)
    if read.query_name in nameldel:
        f2.write(read)
    if read.query_name in namesins:
        f5.write(read)
    if read.query_name in namelins:
        f4.write(read)
    if read.query_name in nameinv:
        f6.write(read)
    if read.query_name in namedup:
        f7.write(read)
    if read.query_name in namehdr:
        f9.write(read)
    if read.query_name in namefull:
        f10.write(read)
    if read.query_name in namesnp:
        f8.write(read)
    if read.query_name in nameunmodify:
        f1.write(read)
f.close()
f1.close();f2.close();f3.close();f4.close();f5.close();f6.close();f7.close();f9.close();f10.close()
f8.close()

os.system('samtools index '+sample+'_split_unmodified.bam')
os.system('samtools index '+sample+'_split_indel_large_del.bam')
os.system('samtools index '+sample+'_split_indel_small_del.bam')
os.system('samtools index '+sample+'_split_indel_small_ins.bam')
os.system('samtools index '+sample+'_split_indel_large_ins.bam')
os.system('samtools index '+sample+'_split_indel_inv.bam')
os.system('samtools index '+sample+'_split_indel_dup.bam')
os.system('samtools index '+sample+'_split_aav_hdr.bam')
os.system('samtools index '+sample+'_split_aav_full.bam')
os.system('samtools index '+sample+'_split_snp.bam')
