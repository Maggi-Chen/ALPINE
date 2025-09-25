#!/usr/bin/env python3
# coding: utf-8

# import modules
import argparse
import os
import time
import pysam
from datetime import datetime
from typing import NamedTuple, Optional

from utils.verify_and_update_hdr_vartag import verify_update_hdr_tag
from utils.verify_and_update_itr_vartag import verify_itr_seq
from utils import debreak_detect


def read_aav_config(config_path):
    
    """
    Read in AAV vector config file and return a directory of AAV info.
    config file should contain 5 columns:
    AAV_ID, AAV_name, left_itr_end_pos, right_itr_begin_pos, insertstart, insert_end
    """

    allvector = open(config_path,'r').read().split('\n')
    allvector = [c for c in allvector if c!='']
    header=allvector[0][1:].split('\t')
    allvector = [c for c in allvector if c[0]!="#"] # Remove header lines if there is any

    aavinfo={}
    for line in allvector:
        line=line.split('\t')
        if line[0] not in aavinfo:
            aavinfo[line[0]]={}
        aavinfo[line[0]][line[1]]={}
        for i in range(len(header)-2):
            if header[i+2]=='Ref_Type':
                aavinfo[line[0]][line[1]][header[i+2]]=line[i+2]
            else:
                if line[i+2]!='NA':
                    aavinfo[line[0]][line[1]][header[i+2]]=int(line[i+2])
                else:
                    aavinfo[line[0]][line[1]][header[i+2]]=None
    return aavinfo


def process_read_targetsite(allread, 
                            chrom, 
                            pos, 
                            filename,
                            read_tag,
                            allindel,
                            covering_bp):
    """
    Detect structural variant events within target regions from sequence alignment records.
    combine read-level and segment-level analyses to identify SV and small variants.
    """
    f = open('temp_'+filename+'_insseq.fa','a')
    segmentreads = {}
    for read in allread:
        if read.is_secondary:
            continue
        read_tag[read.query_name] = ''
        if read.reference_start < pos-10 and read.reference_end > pos+10:
            covering_bp.add(read.query_name)
        if True:
            readinfo = debreak_detect.cigardeletion(
                read.flag, chrom, read.reference_start, read.cigarstring, 1, 100000)
            indelinfo = readinfo[0]
            if read.has_tag('SA'):
                if read.is_supplementary:
                    simpleinfo = [read.query_name, read.flag, chrom, read.reference_start, read.reference_end,
                                  readinfo[2], read.mapping_quality, '']
                else:
                    simpleinfo = [read.query_name, read.flag, chrom, read.reference_start, read.reference_end,
                                  readinfo[2], read.mapping_quality, read.query_sequence]
                if read.query_name not in segmentreads:
                    segmentreads[read.query_name] = [simpleinfo]
                else:
                    segmentreads[read.query_name] += [simpleinfo]
            goodindel = []
            insid = 1
            for sv in indelinfo:
                if 'D-cigar' in sv:
                    if min(pos+20, sv[1]+sv[2])-max(sv[1], pos-20) > 0:
                        goodindel += [sv]
                if 'Unmodified-with-SNP' in sv:
                    if pos-20 <= sv[1] <= pos+20:
                        goodindel += [sv]
                if 'I-cigar' in sv:
                    if pos-20 < sv[1] < pos+20:
                        goodindel += [sv]
                        if sv[2] > 50:
                            insertseq = read.query_sequence[sv[4]:sv[4]+sv[2]]
                            f.write('>'+read.query_name+'_cigarins_' +
                                    str(insid)+'\n'+insertseq+'\n')

            if read.query_name not in allindel:
                allindel[read.query_name] = goodindel
            else:
                allindel[read.query_name] += goodindel
    f.close()

    # detect sv from segment
    for readname in segmentreads:
        readgroup = segmentreads[readname]
        if len(readgroup) < 2:
            continue
        segmentsv = debreak_detect.segmentdeletion_tumor(readgroup, 50, 100000)
        goodsegment = []
        for sv in segmentsv:
            sv = sv.split('\t')
            if 'D-segment' in sv or 'DUP-segment' in sv or 'INV-segment' in sv:
                if min(pos+10, int(sv[1])+int(sv[2]))-max(int(sv[1]), pos-10) > 0:
                    goodsegment += [[sv[0], int(sv[1]), int(sv[2]), sv[3]]]
            if 'I-segment' in sv:
                if pos-10 < int(sv[1]) < pos+10:
                    goodsegment += [[sv[0], int(sv[1]), int(sv[2]), sv[3]]]
        if readname in allindel:
            allindel[readname] += goodsegment
        else:
            allindel[readname] = goodsegment
    f = open('allindel_readname_' + filename + '.txt', 'w')
    for c in allindel:
        f.write(c+'\t'+str(allindel[c])+'\n')
    f.close()

    return (read_tag, allindel, covering_bp)


def process_read_aav_hdr(aavsite, 
                     aav_tag,
                     aav_id,
                     ref_name,
                     aavinfo,
                     aavlen,
                     itrlen,
                     filename):
    gene_start=aavinfo[aav_id][ref_name]['gene_start']
    gene_end=aavinfo[aav_id][ref_name]['gene_end']

    f=open('temp_'+filename+'_realignWT.fa','a')

    for read in aavsite:
        if read.is_secondary:
            continue
        align_tuple=read.get_aligned_pairs()
        (itr1,gene,itr2) = get_itr_hdr_length(align_tuple,None,None,None,None,gene_start,gene_end)
        if gene >10:
            if read.query_name not in aav_tag:
                aav_tag[read.query_name] = ['HDR-'+aav_id]
            elif 'HDR-'+aav_id not in aav_tag[read.query_name]:
                aav_tag[read.query_name] += ['HDR-'+aav_id]
            if read.query_name not in aavlen:
                aavlen[read.query_name]={}
            if aav_id not in aavlen[read.query_name] or gene > aavlen[read.query_name][aav_id]:
                aavlen[read.query_name][aav_id] = gene
        else:
            if not read.is_supplementary:
                f.write('>'+read.query_name+'\n'+read.query_sequence+'\n')
    f.close()
    return (aav_tag, aavlen, itrlen)

def get_itr_hdr_length(align_tuple,five_prime_itr_start,five_prime_itr_end,three_prime_itr_start,three_prime_itr_end,gene_start,gene_end):
    itr1_spos=None;itr1_epos=None
    itr2_spos=None;itr2_epos=None
    gene_spos=None;gene_epos=None
    if five_prime_itr_start==None and five_prime_itr_end==None:
        five_prime_itr_start=-1
        five_prime_itr_end=-1

    if three_prime_itr_start==None and three_prime_itr_end==None:
        three_prime_itr_start=99999999
        three_prime_itr_end=99999999

    for pair in align_tuple:
        if pair[1]==None or pair[0]==None:
            continue
        if pair[1]<five_prime_itr_start:
            continue
        if five_prime_itr_start<=pair[1]<=five_prime_itr_end:
            if itr1_spos == None:
                itr1_spos=pair[0]
            if itr1_epos==None or pair[1]>itr1_epos:
                itr1_epos=pair[0]
        elif gene_start<=pair[1]<=gene_end:
            if gene_spos == None:
                gene_spos=pair[0]
            if gene_epos ==None or pair[1]>gene_epos:
                gene_epos=pair[0]
        elif three_prime_itr_start<=pair[1]<=three_prime_itr_end:
            if itr2_spos == None:
                itr2_spos=pair[0]
            if itr2_epos==None or pair[1]>itr2_epos:
                itr2_epos=pair[0]
        if pair[1]>three_prime_itr_end:
            break

    if itr1_spos==None and itr1_epos==None:
        itr1=0
    elif itr1_spos==None or itr1_epos==None:
        itr1=0
    else:
        itr1=itr1_epos-itr1_spos

    if itr2_spos==None and itr2_epos==None:
        itr2=0
    elif itr2_spos==None or itr2_epos==None:
        itr2=0
    else:
        itr2=itr2_epos-itr2_spos

    if gene_spos==None and gene_epos==None:
        gene=0
    elif gene_spos==None or gene_epos==None:
        gene=0
    else:
        gene=gene_epos-gene_spos

    return (itr1,gene,itr2)

def process_read_aav_itr(aavsite, 
                     aav_tag,
                     aav_id,
                     ref_name,
                     aavinfo,
                     aavlen,
                     itrlen,
                     filename):

    five_prime_itr_start=aavinfo[aav_id][ref_name]['ITR1_start']
    five_prime_itr_end=aavinfo[aav_id][ref_name]['ITR1_end']
    three_prime_itr_start=aavinfo[aav_id][ref_name]['ITR2_start']
    three_prime_itr_end=aavinfo[aav_id][ref_name]['ITR2_end']
    gene_start=aavinfo[aav_id][ref_name]['gene_start']
    gene_end=aavinfo[aav_id][ref_name]['gene_end']

    f=open('temp_'+filename+'_realignWT.fa','a')

    for read in aavsite:
        if read.is_secondary:
            continue

        align_tuple=read.get_aligned_pairs()
        (itr1,gene,itr2) = get_itr_hdr_length(align_tuple,five_prime_itr_start,five_prime_itr_end,three_prime_itr_start,three_prime_itr_end,gene_start,gene_end)

        if itr1>10 or itr2>10:
            if read.query_name not in aav_tag:
                aav_tag[read.query_name] = ['ITR_integration-'+aav_id]
            elif 'ITR_integration-'+aav_id not in aav_tag[read.query_name]:
                aav_tag[read.query_name] += ['ITR_integration-'+aav_id]
            if read.query_name not in itrlen:
                itrlen[read.query_name]={}
            aavinslen=itr1+itr2+gene
            if aav_id not in itrlen[read.query_name] or itrlen[read.query_name][aav_id] < aavinslen:
                itrlen[read.query_name][aav_id] = aavinslen
        elif gene >10:
            if read.query_name not in aav_tag:
                aav_tag[read.query_name] = ['HDR-'+aav_id]
            elif 'HDR-'+aav_id not in aav_tag[read.query_name]:
                aav_tag[read.query_name] += ['HDR-'+aav_id]
            if read.query_name not in aavlen:
                aavlen[read.query_name]={}
            if aav_id not in aavlen[read.query_name] or gene > aavlen[read.query_name][aav_id]:
                aavlen[read.query_name][aav_id] = gene
        else:
            if not read.is_supplementary:
                f.write('>'+read.query_name+'\n'+read.query_sequence+'\n')
    f.close()
    return (aav_tag, aavlen, itrlen)


def realign_to_wt(filename,wtchrom,wtpos,read_tag,allindel,covering_bp,reffile):
    refseq=open(reffile,'r').read().split('>')[1:]
    for chrom in refseq:
        if wtchrom==chrom.split('\n')[0].split(' ')[0]:
            f=open('temp_'+filename+'_wtref.fa','w')
            f.write('>'+wtchrom+'\n'+''.join(chrom.split('\n')[1:])+'\n')
            f.close()
            break
    os.system('minimap2 -ax map-hifi -t 8 --secondary=no temp_'+filename+'_wtref.fa temp_' +
              filename+'_realignWT.fa |samtools sort -o temp_'+filename+'_realignWT.bam')
    os.system('samtools index temp_'+filename+'_realignWT.bam')

    f=pysam.AlignmentFile('temp_'+filename+'_realignWT.bam','rb')
    realigned=f.fetch(wtchrom)
    (read_tag, allindel, covering_bp) = process_read_targetsite(realigned,
                            wtchrom,
                            wtpos,
                            filename,
                            read_tag,
                            allindel,
                            covering_bp)
    return (read_tag, allindel, covering_bp)


def svsize(a):  # sorting a list based on second element
    return a[2]


def find_aav_insertion(samname,aavinfo,aavlen,aav_tag,itrlen):
    try:
        f = pysam.AlignmentFile(samname, 'r')
    except:
        return (aav_tag, aavlen, itrlen)

    allread = f.fetch()

    vector_names={}
    for aav_id in aavinfo:
        for refname in aavinfo[aav_id]:
            vector_names[refname]=aav_id

    goodaav =0
    fullaav = 0
    for read in allread:
        readname = read.query_name.split('_cigarins')[0]
        if read.flag == 4 or read.reference_name not in vector_names:
            continue
        aav_id = vector_names[read.reference_name]
        align_tuple=read.get_aligned_pairs()
        (itr1,gene,itr2) = get_itr_hdr_length(align_tuple,
                                        aavinfo[aav_id][read.reference_name]['ITR1_start'],
                                        aavinfo[aav_id][read.reference_name]['ITR1_end'],
                                        aavinfo[aav_id][read.reference_name]['ITR2_start'],
                                        aavinfo[aav_id][read.reference_name]['ITR2_end'],
                                        aavinfo[aav_id][read.reference_name]['gene_start'],
                                        aavinfo[aav_id][read.reference_name]['gene_end'])

        if itr1>10 or itr2>10:
            fullaav+=1
            if readname not in aav_tag:
                aav_tag[readname] = ['ITR_integration-'+aav_id]
            elif 'ITR_integration-'+aav_id not in aav_tag[readname]:
                aav_tag[readname] += ['ITR_integration-'+aav_id]
            aavinslen = itr1+gene+itr2
            if readname not in itrlen:
                itrlen[readname]={}
            if aav_id not in itrlen[readname] or itrlen[readname][aav_id] < aavinslen:
                itrlen[readname][aav_id] = aavinslen
        else:
            goodaav+=1
            if readname not in aav_tag:
                aav_tag[readname] = ['HDR-'+aav_id]
            elif 'HDR-'+aav_id not in aav_tag[readname]:
                aav_tag[readname] += ['HDR-'+aav_id]

            aavinslen = gene
            if readname not in aavlen:
                aavlen[readname]={}
            if aav_id not in aavlen[readname]  or aavlen[readname][aav_id] < aavinslen:
                aavlen[readname][aav_id] = aavinslen

    print('rescue ins aav: ', goodaav, fullaav)
    return (aav_tag, aavlen, itrlen)

 
def find_variant_clip(temp_unmodified, 
                      bampath, 
                      filename, 
                      fasta, 
                      chrom, 
                      pos,
                      aav_tag, 
                      aavinfo,
                      aavlen,
                      itrlen):
    vector_names={}
    for aav_id in aavinfo:
        for refname in aavinfo[aav_id]:
            vector_names[refname]=aav_id

    primaryalign = {}
    # extract clipped sequence for unmodified reads with >100bp clip on either side
    f = pysam.AlignmentFile(bampath, 'rb')
    allread = f.fetch(chrom)
    tempfile = open('temp_'+filename+'_clippedseq.fa', 'w')
    for read in allread:
        if read.query_name not in temp_unmodified \
                or read.flag > 16 or read.has_tag('SA'):
            continue
        readcigar = read.cigartuples
        ifclip = False
        cigarinfo = [0, 0, 0]
        if readcigar[0][0] == 4:
            cigarinfo[0] = readcigar[0][1]
            if readcigar[0][1] > 100:
                tempfile.write('>'+read.query_name+'_leftclip\n' +
                               read.query_sequence[:readcigar[0][1]]+'\n')
                ifclip = True
        if readcigar[-1][0] == 4:
            cigarinfo[2] = readcigar[-1][1]
            if readcigar[-1][1] > 100:
                tempfile.write('>'+read.query_name+'_rightclip\n' +
                               read.query_sequence[-readcigar[0][1]:]+'\n')
                ifclip = True
        cigarinfo[1] = read.query_alignment_length
        if ifclip:
            primaryalign[read.query_name] = [[read.query_name, read.flag, read.reference_name,
                                              read.reference_start, read.reference_end, cigarinfo, read.mapping_quality]]
    tempfile.close()

    # align clip sequence & detect variant or aav integration
    os.system('minimap2 -ax map-hifi -t 8 --secondary=no ' + fasta + ' temp_' +
              filename+'_clippedseq.fa > temp_'+filename+'_clippedseq.sam')
    try:
        f = pysam.AlignmentFile('temp_'+filename+'_clippedseq.sam', 'r')
    except:
        return (aav_tag, [], {}, aavlen, itrlen)
    allread = f.fetch()
    readtag = {}
    goodaav = []
    realunmodi = []
    fullaav = []
    variant = []
    svinfo = {}
    for read in allread:
        readname = read.query_name.split('_rightclip')[0].split('_leftclip')[0]
        if read.flag == 4:
            readtag[readname] = 'Unmodify'
            continue
        if readname not in readtag and (read.reference_name not in vector_names and read.reference_name !=chrom):
            readtag[readname] = 'Unmodify'
            continue
        # contain aav sequence
        if read.reference_name in vector_names:
            aav_id = vector_names[read.reference_name]
            align_tuple=read.get_aligned_pairs()
            (itr1,gene,itr2) = get_itr_hdr_length(align_tuple,
                                        aavinfo[aav_id][read.reference_name]['ITR1_start'],
                                        aavinfo[aav_id][read.reference_name]['ITR1_end'],
                                        aavinfo[aav_id][read.reference_name]['ITR2_start'],
                                        aavinfo[aav_id][read.reference_name]['ITR2_end'],
                                        aavinfo[aav_id][read.reference_name]['gene_start'],
                                        aavinfo[aav_id][read.reference_name]['gene_end'])

            if itr1>10 or itr2>10:
                if readname not in aav_tag:
                    aav_tag[readname] = ['ITR_integration-'+aav_id]
                elif 'ITR_integration-'+aav_id not in aav_tag[readname]:
                    aav_tag[readname] += ['ITR_integration-'+aav_id]
                aavinslen = itr1+gene+itr2
                if readname not in itrlen:
                    itrlen[readname]={}
                if aav_id not in itrlen[readname] or itrlen[readname][aav_id] < aavinslen:
                    itrlen[readname][aav_id] = aavinslen
            else:
                if readname not in aav_tag:
                    aav_tag[readname] = ['HDR-'+aav_id]
                elif 'HDR-'+aav_id not in aav_tag[readname]:
                    aav_tag[readname] += ['HDR-'+aav_id]

                aavinslen = gene
                if readname not in aavlen:
                    aavlen[readname]={}
                if aav_id not in aavlen[readname]  or aavlen[readname][aav_id] < aavinslen:
                    aavlen[readname][aav_id] = aavinslen
        # same chrom
        if read.reference_name == chrom:
            prime = primaryalign[readname][0]
            cigarinfo = [0, 0, 0]
            cigarinfo[0] = read.cigartuples[0][1] if read.cigartuples[0][0] in [
                4, 5] else 0
            cigarinfo[0] = read.cigartuples[0][1] if read.cigartuples[0][0] in [
                4, 5] else 0
            if read.is_reverse:
                flag = 2064
                cigarinfo[2] = read.cigartuples[0][1] if read.cigartuples[0][0] in [
                    4, 5] else 0
                cigarinfo[0] = read.cigartuples[-1][1] if read.cigartuples[-1][0] in [4, 5] else 0
            else:
                flag = 2048
                cigarinfo[0] = read.cigartuples[-1][1] if read.cigartuples[-1][0] in [4, 5] else 0
                cigarinfo[2] = read.cigartuples[0][1] if read.cigartuples[0][0] in [
                    4, 5] else 0
            cigarinfo[1] = read.query_alignment_length
            if 'leftclip' in read.query_name:
                cigarinfo[2] += sum(prime[5])-sum(cigarinfo)
            primaryalign[readname] += [[readname, flag, chrom, read.reference_start,
                                        read.reference_end, cigarinfo, read.mapping_quality]]


    for readname in primaryalign:
        if len(primaryalign[readname]) > 1 and (readname not in readtag or readtag[readname] not in ['FULL', 'HDR']):
            readsv = debreak_detect.segmentdeletion_tumor(
                primaryalign[readname], 50, 100000)
            goodsv = []
            for c in readsv:
                c = c.split('\t')
                if c[0] == chrom and min(int(c[1])+int(c[2]), pos+10)-max(pos-10, int(c[1])) >= 0:
                    goodsv += [[c[0], int(c[1]), int(c[2]), c[3]]]
            if len(goodsv) == 0:
                continue
            svinfo[readname] = goodsv
            readtag[readname] = 'SV'

    for readname in readtag:
        if readtag[readname] == 'HDR':
            goodaav += [readname]
        if readtag[readname] == 'SV':
            variant += [readname]
        if readtag[readname] == 'FULL':
            fullaav += [readname]
        if readtag[readname] == 'Unmodify':
            realunmodi += [readname]
    print('rescue ins aav: ', len(goodaav), len(fullaav))
    print('rescue variant: ', len(variant))
    f = open('rescued_clip_variant', 'w')
    for c in svinfo:
        goodsv = svinfo[c]
        goodsv.sort(key=svsize, reverse=True)
        f.write(str(goodsv[0])+'\n')
    f.close()
    return (aav_tag, variant, svinfo, aavlen, itrlen)

def find_longest(calls, lengths):
    aav_id = calls[0].split('-')[-1]
    longest = lengths[aav_id]
    for call in calls:
        if lengths[call.split('-')[-1]] > longest:
            aav_id = call.split('-')[-1]
            longest = lengths[call.split('-')[-1]]

    return (aav_id,longest)

def merge_aav_tag(read_tag, aavinfo, aav_tag, aavlen, itrlen, truncated_hdr_cutoff):
    final_aavlen={}
    fulllength={}
    for aav_id in aavinfo:
        if aav_id not in fulllength:
            fulllength[aav_id]={}
        for refname in aavinfo[aav_id]:
            if aavinfo[aav_id][refname]['Ref_Type']=="HDR":
                fulllength[aav_id]["HDR"]=aavinfo[aav_id][refname]['gene_end']-aavinfo[aav_id][refname]['gene_start']
            elif aavinfo[aav_id][refname]['Ref_Type']=="ITR":
                fulllength[aav_id]["ITR"]=aavinfo[aav_id][refname]['gene_end']-aavinfo[aav_id][refname]['gene_start'] + aavinfo[aav_id][refname]['ITR1_end']-aavinfo[aav_id][refname]['ITR1_start'] + aavinfo[aav_id][refname]['ITR2_end']-aavinfo[aav_id][refname]['ITR2_start'] 
    for aav_id in fulllength:
        if "HDR" not in fulllength[aav_id]:
            for refname in aavinfo[aav_id]:
                if aavinfo[aav_id][refname]['Ref_Type']=="ITR":
                    fulllength[aav_id]["HDR"]=aavinfo[aav_id][refname]['gene_end']-aavinfo[aav_id][refname]['gene_start']

    for readname in aav_tag:
        hdr_calls = [c for c in aav_tag[readname] if 'HDR' in c]
        if len(hdr_calls)>0:
            (aav_id, aav_ins_len) = find_longest(hdr_calls, aavlen[readname])
            if truncated_hdr_cutoff * (fulllength[aav_id]["HDR"]) <= aav_ins_len:
                read_tag[readname] = 'HDR-'+aav_id
            else:
                read_tag[readname] = 'Non-HDR-without-ITR-'+aav_id
            final_aavlen[readname] = aav_ins_len
            continue

        itr_calls = [c for c in aav_tag[readname] if 'ITR' in c]
        if len(itr_calls)>0:
            (aav_id, itr_ins_len) = find_longest(itr_calls, itrlen[readname])
            if truncated_hdr_cutoff * (fulllength[aav_id]["ITR"]) <= itr_ins_len:
                read_tag[readname] = 'Non-HDR-with-ITR-'+aav_id
            else:
                read_tag[readname] = 'Non-HDR-with-ITR-'+aav_id
            final_aavlen[readname] = itr_ins_len
            continue

    return (read_tag, final_aavlen)


def find_unmapped_reads(bampath, read_tag):
    f = pysam.AlignmentFile(bampath, 'rb')
    allread = f.fetch(until_eof=True)
    unmapped = []
    for read in allread:
        if read.query_name not in read_tag and read.query_name not in unmapped:
            unmapped += [read.query_name]
    return unmapped


def assign_read(read_tag,
                aav_tag, 
                allindel, 
                filename, 
                fasta, 
                aavinfo,
                bampath, chrom,
                pos, covering_bp,
                aavlen, 
                itrlen,
                full_hdr_cutoff):
    # realign inserted seq from INS calls - find aav sequence
    os.system('minimap2 -ax map-hifi -t 8 --secondary=no ' + fasta + ' temp_' +
              filename+'_insseq.fa > temp_'+filename+'_insseq.sam')
    (aav_tag, aavlen, itrlen) = find_aav_insertion('temp_'+filename+'_insseq.sam',
                                                    aavinfo,
                                                    aavlen,
                                                    aav_tag,
                                                    itrlen)

    # realign unmapped clip sequence from Unmodified reads - find aav sequence
    temp_unmodified = []
    for readname in read_tag:
        if read_tag[readname] == '' and (readname not in allindel or len(allindel[readname]) == 0):
            temp_unmodified += [readname]
    temp_unmodified = set(temp_unmodified)
    (aav_tag, variant, svinfo, aavlen, itrlen) = find_variant_clip(temp_unmodified, 
                                                                      bampath,
                                                                      filename,
                                                                      fasta, 
                                                                      chrom,
                                                                      pos, 
                                                                      aav_tag,
                                                                      aavinfo,
                                                                      aavlen,
                                                                      itrlen)
    variant = set(variant)

    (read_tag, final_aavlen) = merge_aav_tag(read_tag, aavinfo, aav_tag, aavlen, itrlen, full_hdr_cutoff)

    # write output read-group txt file
    unmodified = 0
    small_del = 0
    small_ins = 0
    large_del = 0
    large_ins = 0
    inv = 0
    dup = 0
    hdr = 0
    truncated = 0
    full = 0
    snp = 0
    unknown=0

    unknown_readnames = find_unmapped_reads(bampath, read_tag)

    f = open('temp_readname_'+filename+'.txt', 'w')
    for readname in read_tag:
        if read_tag[readname] == '':
            if readname not in allindel or len(allindel[readname]) == 0:
                if readname not in variant:
                    if readname in covering_bp:
                        f.write(readname+'\tUnmodified\t\n')
                        unmodified += 1
                    else:
                        unknown+=1
                        f.write(readname+'\tUnclassified\t\n')
                    continue
                allsv = svinfo[readname]
            else:
                allsv = allindel[readname]
            if len(allsv) == 1:
                sv = allsv[0]
            else:
                allsv.sort(key=svsize, reverse=True)
                sv = allsv[0]
            if 'D-' in sv[3]:
                if sv[2] >= 50:
                    large_del += 1
                    f.write(readname+'\tDEL-large\t'+str(sv[2])+'\n')
                else:
                    small_del += 1
                    f.write(readname+'\tDEL-small\t'+str(sv[2])+'\n')
            if 'I-' in sv[3]:
                if sv[2] >= 50:
                    large_ins += 1
                    f.write(readname+'\tINS-large\t'+str(sv[2])+'\n')
                else:
                    small_ins += 1
                    f.write(readname+'\tINS-small\t'+str(sv[2])+'\n')
            if 'INV-' in sv[3]:
                inv += 1
                f.write(readname+'\tINV\t'+str(sv[2])+'\n')
            if 'DUP-' in sv[3]:
                dup += 1
                f.write(readname+'\tDUP\t'+str(sv[2])+'\n')
            if 'Unmodified-with-SNP' in sv[3]:
                snp += 1
                f.write(readname+'\tUnmodified-with-SNP\t'+str(sv[2])+'\n')

        else:
            if 'Non-HDR-without-ITR-' in read_tag[readname]:
                f.write(readname+'\t'+read_tag[readname]+'\t'+str(final_aavlen[readname])+'\n')
                truncated += 1
                continue
            if 'HDR-' in read_tag[readname]:
                hdr += 1
                f.write(readname+'\t'+read_tag[readname]+'\t'+str(final_aavlen[readname])+'\n')
                continue
            if 'ITR-' in read_tag[readname]:
                full += 1
                f.write(readname+'\t'+read_tag[readname]+'\t'+str(final_aavlen[readname])+'\n')
                continue

    for readname in unknown_readnames:
        unknown+=1
        f.write(readname+'\tUnclassified\t\n')

    f.close()

    print('number of reads:')
    print('unmodified\t', unmodified)
    print('DEL-large\t'+str(large_del)+'\nDEL-small\t'+str(small_del)+'\nINS-large\t'+str(large_ins) +
         '\nINS-small\t'+str(small_ins)+'\nINV\t'+str(inv)+'\nDUP\t'+str(dup)+
         '\nHDR\t'+str(hdr)+'\nNon-HDR-without-ITR\t'+str(truncated)+'\nITR\t'+str(full))
    print('Unmodified-with-SNP\t'+str(snp))
    print('Unclassified'+str(unknown))


# run workflow
def classify_bam(filename, 
                 bampath, 
                 fasta, 
                 aavinfo,
                 full_hdr_cutoff):
    f = pysam.AlignmentFile(bampath, 'rb')
    read_tag={}
    aavlen={}
    itrlen={}
    aav_tag={}
    covering_bp = set([])
    allindel={}
    wtchrom=''
    pos=''
    for aav in aavinfo:
        for ref_name in aavinfo[aav]:
            alignment=f.fetch(ref_name)
            if aavinfo[aav][ref_name]['Ref_Type'] == "WT":
                wtchrom=ref_name
                pos=(int(aavinfo[aav][ref_name]['gene_start'])+int(aavinfo[aav][ref_name]['gene_end']))//2
                tempfile = open('temp_'+filename+'_insseq.fa', 'w') # This is to clear previous existing temp file
                (read_tag, allindel, covering_bp) = process_read_targetsite(alignment,
                            ref_name,
                            pos,
                            filename,
                            read_tag,
                            allindel,
                            covering_bp)
            elif aavinfo[aav][ref_name]['Ref_Type'] == "ITR":
                (aav_tag, aavlen, itrlen) = process_read_aav_itr(alignment, 
                                              aav_tag,
                                              aav,
                                              ref_name,
                                              aavinfo,
                                              aavlen,
                                              itrlen,
                                              filename)
            else:
                (aav_tag, aavlen, itrlen) = process_read_aav_hdr(alignment,
                                              aav_tag,
                                              aav,
                                              ref_name,
                                              aavinfo,
                                              aavlen,
                                              itrlen,
                                              filename)
    # realign WT-like reads to WT only reference
    (read_tag, allindel, covering_bp) = realign_to_wt(filename,
                                            wtchrom,
                                            pos,
                                            read_tag,
                                            allindel,
                                            covering_bp,
                                            fasta)

    read_tag = assign_read(read_tag,
                           aav_tag,
                           allindel,
                           filename,
                           fasta,
                           aavinfo,
                           bampath,
                           wtchrom,
                           pos,
                           covering_bp,
                           aavlen,
                           itrlen,
                           full_hdr_cutoff)



# --------------------------------------------------------------------------------------

class Args(NamedTuple):
    """command line arguments"""
    input: str
    fastq: str
    prefix: str
    ref: str
    vect_name: str
    left_itr_end_pos: int
    right_itr_start_pos: int
    ins_start: int
    ins_end: int
    full_hdr_cutoff: int
    five_prime_HA_arm_seq: str
    three_prime_HA_arm_seq: str
    ha_match_ratio: Optional[float] = 0.98
    left_itr_seq: Optional[str] = None
    right_itr_seq: Optional[str] = None
    seed_size: Optional[int] = 15
    perc_identity: Optional[int] = 80


def get_args() -> Args:
    """ Get command-line arguments """
    parser = argparse.ArgumentParser(usage="Identify and classify variant types. \n-Author:\tMaggie Chen, Xing-Huang Gao BMS 2023(R)",
                                     description="Identify full-length, KI and INDELS",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter
                                     )

    parser.add_argument("-i", "--input", type=str, required=True,
                        help="input bam file")
    parser.add_argument("-fq", "--fastq", type=str, required=True,
                        help="filtered fastq file")
    parser.add_argument("-p", "--prefix", type=str, required=True,
                        help="sample prefix")
    parser.add_argument("-r", "--ref", type=str, required=True,
                        help="reference genome fa file including vectors")
    parser.add_argument("--config", type=str, required=True,
                        help="config file for AAV vector")
    parser.add_argument("--full_hdr_cutoff", type=float, required=False, default=0.99,
                        help="length cutoff (percentage) for full-length aav integration")
    parser.add_argument("--five_prime_HA_arm_seq", type=str, required=True,
                        help="five prime HA arm sequence")
    parser.add_argument("--three_prime_HA_arm_seq", type=str, required=True,
                        help="three prime HA arm sequence")
    parser.add_argument("--left_itr_seq",
                        help="left ITR sequence")
    parser.add_argument("--right_itr_seq",
                        help="right ITR sequence")
    parser.add_argument("--ha_match_ratio", type=float,
                        default=0.98, help="HA sequence match ratio")
    parser.add_argument("--seed_size", default=15, help="word size for initial sequence match")
    parser.add_argument("--perc_identity", default=80, help="percentage identity for sequence alignment")
    return parser.parse_args()


def main() -> None:
    args = get_args()
    
    bampath = args.input
    fastq = args.fastq
    filename = args.prefix
    fasta = args.ref

    full_hdr_cutoff = args.full_hdr_cutoff
    left_ha_seq = args.five_prime_HA_arm_seq
    right_ha_seq = args.three_prime_HA_arm_seq
    left_itr_seq = args.left_itr_seq
    right_itr_seq = args.right_itr_seq
    ha_match_ratio = args.ha_match_ratio
    seed_size = args.seed_size
    perc_identity = args.perc_identity

    
    readname_file = f"temp_readname_{filename}.txt"

    aavinfo = read_aav_config(args.config)
    

    # Get current date and time
    current_datetime = datetime.now()
    current_date = current_datetime.date()
    current_time_formatted = current_datetime.strftime("%H:%M:%S")

    output_str = "# Starting..."
    print(f">{output_str:=^200}<\n{current_date}, {current_time_formatted}")
    start_time = time.time()

    classify_bam(filename,
                 bampath,
                 fasta,
                 aavinfo,
                 full_hdr_cutoff)
    
    verify_update_hdr_tag(readname_file, fastq, left_ha_seq,
                          right_ha_seq, f'{filename}_ValidatedHA', 'HDR', ha_match_ratio)
    
    if left_itr_seq is not None and right_itr_seq is not None:
        print("# Verifying ITR sequences ...")
        verify_itr_seq(f'readname_{filename}_ValidatedHA.txt', fastq, filename, left_itr_seq,
                   right_itr_seq, seed_size=seed_size, percent_identity=perc_identity)
    else:
        print("# No ITR sequences provided. Skipping ITR sequence verification.")
        os.rename(f'readname_{filename}_ValidatedHA.txt', f'readname_{filename}.txt')

    
    end_time = time.time()
    print(f'# {filename}: Done. Total elapsed time: \
          {(end_time - start_time)/60:.1f} minutes.')



if __name__ == '__main__':
    main()

