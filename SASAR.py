#!/usr/bin/env python3
"""
# June 2021
# If using this pipeline please cite : XXXXXXXXXX
#--------------------------------------------------------------------------+    
#                                                   
#	SASAR is a meta-assembly tool 
#       to reconcile different long read assemblies without a reference guide.  
#	Example commands:    
#	python SASAR.py in_dir (option -t / -i/ -c/ -m/ -r / -o)                                              
#                                                        
#--------------------------------------------------------------------------+
#                                                      
#	AUTHOR: panpan ZHANG                            
#	CONTACT: njaupanpan@gmail.com                      
#                                                         
#	LICENSE:                                            
# 	GNU General Public License, Version 3               
#	http://www.gnu.org/licenses/gpl.html  
#                                             
#	VERSION: V.1                    
#                                                                                                       
#--------------------------------------------------------------------------+
"""

import argparse
import re
import sys
import os
import pathlib
import gzip
import subprocess
import itertools
import multiprocessing
import pybedtools
import pandas as pd
from pybedtools import BedTool
from functools import reduce
#import tempfile
#import collections
#from operator import index
#import random
#import shutil

__version__ = '1.0'

if str(pd.__version__) == "":
    print("Please install pandas") 

if str(pybedtools.__version__)== "":
    print("Please install pybedtools") 


def get_arguments(args):
    parser = argparse.ArgumentParser(description='Long read assembly reconciliation', add_help=False, 
                            usage="python SASAR.py in_dir")

    required_args = parser.add_argument_group('Positional arguments')
    required_args.add_argument('in_dir', type=str,
                            help='input directory containing all assemblies')
    #required_args.add_argument('out_dir', type=str,
    #                        help='output directory for SASAR-assembly')

    setting_args = parser.add_argument_group('Settings')
    setting_args.add_argument('-t', metavar="INT",type=int, default=get_default_thread(),
                            help='number of CPU threads for whole genome alignment')
    setting_args.add_argument('-i', metavar="INT",type=int,default=95,
                            help="minimum identity confidence score [95]")
    setting_args.add_argument('-c', metavar="INT",type=int,default=95,
                            help="minimum coverage confidence score [95]")
    setting_args.add_argument('-m',metavar="INT", type=int,default=50000,
                            help="minium overlap length [50000]")
    setting_args.add_argument('--repeat_size',metavar="INT", type=int,default=50000,
                            help="repeat size [50000]")

    output_args = parser.add_argument_group("Output options")
    output_args.add_argument("-o", metavar="PATH", type=str, default="SASAR_output", 
                            help="output directory for SASAR-assembly [./SASAR_output]")

    other_args = parser.add_argument_group('Other')
    other_args.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                            help='show this help message and exit')
    other_args.add_argument('-v','--version', action='version',
                            version='SASAR v' + __version__,
                            help="show program's version number and exit")

    args = parser.parse_args(args)
    return args

def main(args=None):
    args = get_arguments(args)
    #parser = get_parser()
    #random.seed(0)
    out_dir = args.o
    if not os.path.isdir(out_dir):
        os.mkdir(out_dir)
    out_dir = os.path.abspath(out_dir) + "/"


    wga=process_wga(args.in_dir,args.o,args.t)
    filter=id_filter(args.i,args.c,args.m,args.repeat_size,args.o)
    #pans=Pan_contig(args.o)
    #ends=Extend_ends(args.i,args.c,args.m,args.o)
    #os.makedirs(args.o, exist_ok=True)
    print('SASAR to {}'.format(args.o))
    print()

def process_wga(in_dir,o,t):
    l = [str(x) for x in sorted(pathlib.Path(in_dir).glob('*'))
                      if x.is_file()]
    #print(l)
    l = [x for x in l if
                      x.endswith('.fasta') or x.endswith('.fasta.gz') or
                      x.endswith('.fna') or x.endswith('.fna.gz') or
                      x.endswith('.fa') or x.endswith('.fa.gz')]  
    #print(l)
    print("---------------------------------------"+'Total {:,} assemblies as input'.format(len(l))+"----------------------------------------------------")
    #cmd ="mkdir " +out_dir 
    #subprocess.call(cmd, shell=True) 
    for x in l:
        if x.endswith('.fasta') or x.endswith('.fna') or x.endswith('.fa'):
            fname=x.split('/')[-1]
            outfa=o+"/NEW_"+str(fname)
            print("---------------------------------------Reheader assembly : "+str(fname)+" --------------------------------------------------------------")
            cmd ="bioawk -cfastx '{ print \">\"$name\"_len\"length($seq);print $seq}' " + x +"  >" +outfa
            #subprocess.call(cmd, shell=True) 
        if x.endswith('.fasta.gz') or x.endswith('.fna.gz') or x.endswith('.fa.gz'):
            fname=x.split('/')[-1].split('.gz')[0]
            outfa=o+"/NEW_"+str(fname)
            print("---------------------------------------Reheader assembly : "+str(fname)+" --------------------------------------------------------------")
            cmd ="bioawk -cfastx '{ print \">\"$name\"_len\"length($seq);print $seq}' " + x +"  >" +outfa
            #subprocess.call(cmd, shell=True) 
    l = [str(x) for x in sorted(pathlib.Path(o).glob('NEW_*'))
                      if x.is_file()]
    print("---------------------------------------Whole genome alignment of every two assemblies----------------------------------------------------")
    cmd='cat ' +' '.join([str(i) for i in l]) + ' > ' +o+'/WGA.fa'
    print(cmd)
    subprocess.call(cmd, shell=True)
    os.chdir(o)
    length = len(l)
    i = 0
    # Iterating 
    while i < length:
        r = l[i].split('/')[-1]
        i += 1
        j=i
        while j < length:
            q = l[j].split('/')[-1]
            j += 1
            print(r,q)
            outpaf="Q"+str(r)+"_"+str(q)+".paf" 
            minimap2_command = ["minimap2", "--cs", "-cxasm20", "-t", str(t),r,q,"-o",outpaf]
            print(minimap2_command)
            #minimap2_out = subprocess.run(minimap2_command, stdout=subprocess.PIPE).stdout.decode()
    cmd = "awk 1 Q*.paf >WGA.paf"  
    #subprocess.call(cmd, shell=True)  


    print("---------------------------------------Get metrics from alignment---------------------------------------------------------------------------")
    paf=open('WGA.paf')    
    outfile=open('WGA.statsF','w') 
    #--------------------------------------------------------------------------------------------------------------------------------------------------+
    #   This script is to get the statistics from paf file, such as Insertion, Deletion, Subsititution 
    #   and Gap-compressed Identity by the same definition from minimap2 developer 
    #   Heng Li's blog: http://lh3.github.io/2018/11/25/on-the-definition-of-sequence-identity
    #		       
    #	Tag  Type                       Description                      
    #
    #	tp   A    Type of aln: P/primary, S/secondary and I,i/inversion 
    #	cm   i    Number of minimizers on the chain                     
    #    s1   i    Chaining score
    #    s2   i    Chaining score of the best secondary chain
    #    NM   i    Total number of mismatches and gaps in the alignment
    #    MD   Z    To generate the ref sequence in the alignment
    #    AS   i    DP alignment score
    #    ms   i    DP score of the max scoring segment in the alignment
    #    nn   i    Number of ambiguous bases in the alignment
    #    ts   A    Transcript strand (splice mode only)
    #    cg   Z    CIGAR string (only in PAF) M:MATCH; I:iNSERTION; D:DELETION
    #    cs   Z    Difference string
    #-------------------------------------------------------------------------------------------------------------------------------------------------+
    #headers = ['refID','rstart', 'rend', 'rlen','direction','queryID', 'qlen','qstart', 'qend' ,'match','block',
    #'blast_iden','gap_compress_iden','divergence','PS',
    #'rcov', 'qcov','ins_max','delt_max','sub','gap_compress'] 
    headers =  ['queryID','qstart','qend','qlen', 'refID','direction','rstart','rend', 'rlen','gap_compress_iden','PS']
    outfile.write('\t'.join(headers))
    outfile.write('\n')
    total=0
    n_primary=0
    #parts = pd.DataFrame([i.split('\t') for i in pafFile]), if install pandas module
    for line in paf:         
        parts = line.strip().split("\t")
        total=total+1
        #get tag "cg" for cigar
        cg=line.strip().split("cg:Z:")[1].split("\t")[0]
        PS=line.strip().split("tp:A:")[1].split("\t")[0]
        DE=line.strip().split("de:f:")[1].split("\t")[0]
        #value=re.findall(r'(\d+)(\w)', cg)
        M= re.compile(r'(\d+)[M]').findall(cg)
        I= re.compile(r'(\d+)[I]').findall(cg)
        D= re.compile(r'(\d+)[D]').findall(cg)
        match_list =list(map(int, M))
        ins_list = list(map(int, I))
        ins_max=max(ins_list,default=0)
        delt_list = list(map(int, D))
        delt_max=max(delt_list,default=0)
        ins_compress = cg.count("I")
        delt_compress = cg.count("D")
        # NM is edit distance: NM = INS + DEL + SUB,
        # but for minimap2 NM = #mismatches + #I + #D + #ambiguous_bases, NM does not count reference skip `N`.
        NM= int(parts[10])- int(parts[9])
        sub=NM-ins_compress-delt_compress
        ins_compress = cg.count("I")
        delt_compress = cg.count("D")
        blast_iden =  100.0 *int((parts[9]))/ int((parts[10]))
        gap_compress = ins_compress +delt_compress
        gap_compress_iden = 100.0 *int((parts[9]))/(int(parts[10])-gap_compress)
        ql = int(parts[3])-int(parts[2])
        qcov = 100.0 *int(ql)/int((parts[1]))
        rl = int(parts[8])-int(parts[7])
        rcov = 100.0 *int(rl)/int((parts[6]))
        resultss = {
		 "queryID": parts[0],
		 "qlen":  int(parts[1]),
		 "qstart": int(parts[2]),
		 "qend": int(parts[3]),
		 "direction": parts[4],
		 "refID": parts[5],
		 "rlen": int(parts[6]),
		 "rstart": int(parts[7]),
		 "rend": int(parts[8]),
		 "match": int(parts[9]),
		 "block": int(parts[10]),
		 "qcov" : qcov,
		 "rcov" : rcov,
         "PS":PS,
         "divergence":DE,
		 "ins_max": ins_max,
		 "delt_max": delt_max,
		 "sub": sub,
		 "gap_compress":gap_compress,
		 "blast_iden": blast_iden,
		 "gap_compress_iden": gap_compress_iden,
         }
        out_row = (str(resultss[x]) for x in headers)
        outfile.write('\t'.join(out_row))
        outfile.write('\n')
    cmd = "awk -v OFS='\t' 'NR>1{print $5,$7,$8,$9,$1,$6,$2,$3,$4,$10,$11}' WGA.statsF >WGA.statsR;awk 1 WGA.statsF WGA.statsR > WGA.stats;rm WGA.statsF WGA.statsR"
    subprocess.call(cmd, shell=True)
    #filter for primary alignment
    #subprocess.call("awk '$11 ~ /P/' WGA.stats > WGA.statsP", shell=True)
    #cat WGA.stats3 |awk '$15P ~ /P/' |awk '{l[$1"\t"$6]++}END{for (x in l) print x,l[x]}'  > WGA.stats3.txt
    #awk -v FS='\t' -v OFS='\t' 'FNR==NR{A[$1 FS $2]=$3;next}{print ($1 FS $6 in A ) ? $0 OFS A[$1 FS $6] : $0 OFS 0}'  WGA.stats3.txt WGA.stats3 > WGA.stats
    cmd = "cat WGA.statsP |awk -v OFS='\t' 'NR>1{print $1\"__\"$5,$0}' |cut -f1,3-12|bedtools sort |bedtools merge -c 4,5,6,9,10 -o distinct,distinct,distinct,distinct,max -d 10000 -s |sed 's/__/\t/g' |cut -f1,3-9 > WGA.stats.sum"
    #subprocess.call(cmd, shell=True)

class dd_list(dict):
    def __missing__(self,k):
        r = self[k] = []
        return r
D = dd_list()

def id_filter(i,c,m,repeat_size,o):
    print("---------------------------------------Contig identity and coverage filter----------------------------------------------------------------------")
    df=pd.read_csv('WGA.stats', skiprows=1,sep='\t',
    names =  ['queryID','qstart','qend','qlen', 'refID','direction','rstart','rend', 'rlen','identity','PS'])
    if df.empty:
        print('DataFrame is empty!')
    else:
        #filter for primary alignment
        df=df[df.PS == "P"]
        #df=df[(df.queryID == "at.sd_utg668_len3242322") | (df.queryID == "at.wt_ctg107_len67542") ]
        #df=df[(df.queryID == "tig00004409_len418920") | (df.queryID == "at.sd_utg3040_len265040 ") ]
        df['queryID'] = df[['queryID','refID']].apply(lambda x : '{}__{}'.format(x[0],x[1]), axis=1)
        #merge overlap and keep strandedness
        x =BedTool.sort(BedTool.from_dataframe(df))
        y =BedTool.merge(x,s=True, c = '4,5,6,9,10',
        o = 'distinct,distinct,distinct,distinct,max',d=10000)
        d=BedTool.to_dataframe(y,disable_auto_names=True, header=None,names = ['queryID','qstart','qend','qlen','refID','direction','rlen','identity'])
        d['queryID']=d['queryID'].str.split('__', expand = True)[0]
        #sum up query contig coverage and keep strandedness
        d['qcov'] =((d.qend-d.qstart)/d.qlen*100).round(3)
        dcov=d.groupby(['queryID','refID','direction'])['qcov'].sum().reset_index(name='coverage').round(3)
        d=pd.merge(d, dcov, on=['refID','queryID','direction']) 
        #d=d[(d.queryID == "tig00004409_len418920") | (d.queryID == "at.sd_utg3040_len265040") ]
        d= d[(d.qlen * d.coverage/100>10000)].sort_values(['coverage'])
        #print(d.sort_values(['coverage']))
        #add ref start and end information
        df=d[['queryID','qstart', 'qend','refID','direction','coverage']]
        df.columns=['refID','rstart', 'rend','queryID','direction','rcoverage']
        d=pd.merge(d, df, on=['queryID','refID','direction']) 
        #reference contig coverage
        #d['rcov'] =(d.rend -d.rstart)/d.rlen*100
        d.to_csv('WGA.merge', header=False, index = False, sep='\t') 
        
        df0=d.copy()
        #print(df)
        #df1=df[['queryID','qstart','qend','qlen','refID','direction','rlen','identity','qcov','coverage','rstart','rend','rcoverage']]
        
        #print(df0)
        #df0.columns=['refID','rstart','rend','rlen','queryID','direction','qlen','identity','qcov','rcoverage','qstart','qend','coverage']
        #dfr=df0[['queryID','qstart','qend','qlen','refID','direction','rlen','identity','qcov','coverage','rstart','rend','rcoverage']]
        #print(dfr)
        #df=pd.concat([df0,dfr])
        #df=dfr.append(df0)
        #print(df)
        print("---------------------------------------Remove Redundancy of contigs---------------------------------------------------------------------------")
        if c: 
            d = d[ (d.qlen < d.rlen ) & (round(d.coverage,0) >= c)].sort_values(['coverage'],ascending=False)
        if i: 
            d = d[ (d.qlen < d.rlen) & (round(d.identity,0) >= i)]
            #print(d)
            d1=d.groupby(['rlen','refID'])['queryID'].agg(list).reset_index().sort_values(['rlen'],ascending=False)
            #print(d1.dtypes)
            D=dict(zip(d1['refID'],d1['queryID']))
            #print(D)
            count=0
            newDic = dict()
            limit = len(list(D))
            lists=[k for k in D.keys()]
            while count < limit:
                v_num=[v for v in D.values()][int(count)]
                k_num=[k for k in D.keys()][int(count)]
                newDic[k_num] = v_num 
                newDic_v=[v for v in newDic.values()] 
                newDic_v_merge=list(itertools.chain.from_iterable(newDic_v))
                c=count
                while c < limit:
                    k_num_next=[k for k in D.keys()][int(c)]
                    #print(k_num_next)
                    if k_num_next  in newDic_v_merge:
                        #print(k_num_next+" is in the last round")           
                        lists.remove( k_num_next)
                    break
                count+=1
                f_out=open('WGA.r1.txt','w')
                f_out.write("\n".join(str(item) for item in lists))
            #subprocess.call("cat WGA.list.txt |tr ' ' '\n' |sort|uniq > WGA.r1.txt", shell=True)   

#def Extend_ends(i,c,m,o):
        cmd = "seqtk subseq WGA.fa WGA.r1.txt > WGA.r1.fa"  
        subprocess.call(cmd, shell=True)   
    #cmd = "grep -wFf WGA.r1.txt WGA.stats.merge > WGA.stats.r1" 
    #subprocess.call(cmd, shell=True) 
    #print("---------------------------------------Extend two ends of contigs---------------------------------------------------------")
        #names = ['queryID','qstart','qend','qlen','refID','direction','rlen',
        #'identity','qcov', 'coverage', 'rstart', 'rend', 'rcoverage']
        #df0 = pd.read_csv('WGA.stats.merge', sep='\t',names=names,header=None)
        #df=df0[(df0.queryID == "tig00004409_len418920") | (df0.queryID == "at.sd_utg3040_len265040") ]
   
        #print(df)
    #if df.empty:
    #    print('DataFrame stats is empty!')
    #else:
        #print(df)
        
        f=open('WGA.r1.txt', 'r')
        selection = [line.strip() for line in f]
        df=df0[(pd.DataFrame(df0.refID.tolist()).isin(selection).any(1)) | (pd.DataFrame(df0.queryID.tolist()).isin(selection).any(1)) ]
        print(df)
        #df=d0
        #print(df)
        #df_sels=df[~df.refID.isin(dfE_ref.refID.tolist())]
        #dsearch=pd.read_csv('WGA.r1.txt',skiprows=0,sep='\t',names = ['refID'])
        #df=pd.merge(df, dsearch, on=['refID']) 

        #print(df.head(50))
        df=df[(df.queryID.astype(str) == "tig00000133_len7134675") |(df.refID.astype(str) == "tig00000133_len7134675") ]
        #df=df[(df.refID.astype(str) == "tig00000133_len7134675") ]
        df['Eq_end']=df.qlen-df.qend
        df['Er_end']=df.rlen-df.rend
        #print(df)

        df1=df.groupby(['refID','queryID','direction'])['Eq_end'].min().reset_index(name='Eq')
        df2=df.groupby(['refID','queryID','direction'])['Er_end'].min().reset_index(name='Er')
        df3=df.groupby(['refID','queryID','direction'])['qstart'].min().reset_index(name='Sq')
        df4=df.groupby(['refID','queryID','direction'])['rstart'].min().reset_index(name='Sr')
        data_frames=[df1,df2,df3,df4]
        df0 = reduce(lambda  left,right: pd.merge(left,right,on=['refID','queryID','direction'],
                                            how='outer'), data_frames)
        df=pd.merge(df, df0, on=['refID','queryID','direction']) 
        #print(df.sort_values(['coverage']))

        def add_end(df):      
            if (df.direction == '+'):
                if repeat_size:
                    if (df.Sq>=df.Sr):
                        if (df.Er<df.Eq) or (df.Eq< df.Er <repeat_size):
                            return "inside"
                    if (df.Sq<=df.Sr) and (df.Er< repeat_size <df.Eq) :
                        s=int(df.Eq - df.Er)
                        return 'END:{}'.format(s) 
                    if (df.Eq < df.Er) and (df.Sr< repeat_size <df.Sq) :
                        s=int(df.Sq - df.Sr)
                        return 'START:{}'.format(s)
            elif (df.direction == '-'): 
                if repeat_size:
                    if (df.Eq>=df.Sr): 
                        if (df.Sq>df.Er) or (df.Sq< df.Er <repeat_size):
                            return "inside"         
                    if (df.Eq<=df.Sr)and (df.Eq <repeat_size)and (df.Er< repeat_size <df.Sq) :              
                        s=int(df.Sq - df.Er) 
                        return 'END:{}'.format(s) 
                    if (df.Sq < df.Er) and (df.Sq <repeat_size) and (df.Sr< repeat_size <df.Eq):
                        s=int(df.Eq - df.Sr)
                        return 'START:{}'.format(s)            
            return "NONE"
        
        df['Ef_len'] = df.apply(add_end, axis=1)
        #print(df)
        df=df[df.Ef_len.str.contains("START|END")]
        df['INFO']=df.Ef_len.str.split(":").str[0]
        df['Ef_len']=df.Ef_len.str.split(":").str[1]
        #print(df)
        #dsearch.columns=dsearch.columns.str.replace('refID', 'queryID')
        #df=pd.merge(df, dsearch, on=['queryID']).sort_values(['rlen'])
        #print(df.head(50))
        #print(df.sort_values(['queryID']))
        #df['diff1']=(df.qend-df.qstart)/(df.rend-df.rstart)
        #df['diff2']=(df.qend-df.qstart)-(df.rend-df.rstart)
        if m:
            df=df[df.qend-df.qstart > m]
            df= df[df.groupby(['refID','INFO'])['Ef_len'].transform(max) == df['Ef_len']]
            print(df)


def get_default_thread():
    return min(multiprocessing.cpu_count(), 16)

if __name__ == '__main__':
    main()
