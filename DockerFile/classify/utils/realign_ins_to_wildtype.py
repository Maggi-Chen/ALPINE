import pysam
import os

def realign_ins_wt(filename, wt_reference):
    read_class = open("readname_"+filename+".txt","r").read().split('\n')[:-1]
    ins_readname = {}
    for c in read_class:
        if "INS-large" in c:
            ins_readname[c.split('\t')[0]] = c.split('\t')[2]
    ins_map_info={}

    os.system("minimap2 -ax map-hifi -t 8 --secondary=no " + wt_reference + " temp_"+filename+"_insseq.fa > temp_insseq_wt_reference_"+filename+".sam")
    os.system("samtools sort -@ 8 temp_insseq_wt_reference_"+filename+".sam -o InsSeq_WT_reference_"+filename+".bam")
    os.system("samtools index InsSeq_WT_reference_"+filename+".bam")
    os.system("samtools view -f 4 InsSeq_WT_reference_"+filename+".bam | samtools fastq - -0 InsSeq_WT_reference_"+filename+"_MM2_unmapped.fastq")
    os.system("bwa index " + wt_reference)
    os.system("bwa mem -t 8 " + wt_reference + " InsSeq_WT_reference_"+filename+"_MM2_unmapped.fastq > InsSeq_WT_reference_"+filename+"_BWA.sam")

    f=pysam.AlignmentFile("temp_insseq_wt_reference_"+filename+".sam", "r")
    allreads=f.fetch()
    for read in allreads:
        readname = read.query_name.split('_cigarins_')[0]
        if read.is_supplementary or readname not in ins_readname or readname in ins_map_info:
            continue
        if read.flag == 4:
            continue
        mapped_length = read.reference_end-read.reference_start
        mapped_ratio = mapped_length/float(ins_readname[readname])
        ins_map_info[readname] = "\t".join([read.reference_name, str(read.reference_start), 
            str(read.reference_end), str(mapped_length), ins_readname[readname], str(mapped_ratio)])

    f=pysam.AlignmentFile("InsSeq_WT_reference_"+filename+"_BWA.sam", "r")
    allreads=f.fetch()
    for read in allreads:
        readname = read.query_name.split('_cigarins_')[0]
        if read.is_supplementary or readname not in ins_readname or readname in ins_map_info:
            continue
        if read.flag == 4:
            continue
        mapped_length = read.reference_end-read.reference_start
        mapped_ratio = mapped_length/float(ins_readname[readname])
        ins_map_info[readname] = "\t".join([read.reference_name, str(read.reference_start),
            str(read.reference_end), str(mapped_length), ins_readname[readname], str(mapped_ratio)])

    f=open('large_ins_mapping_info_'+filename+'.txt','w')
    columns = ['readname', 'reference', 'reference_start', 'reference_end', 'mapped_length', 'ins_length', 'mapped_ratio']
    
    f.write('\t'.join(columns)+'\n')
    for read in ins_map_info:
        f.write(read+'\t'+ins_map_info[read]+'\n')
    f.close()
