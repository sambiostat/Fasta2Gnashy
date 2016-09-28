#!usr/bin/python
######################################################################################
#Wrapping script for preCC process, produce gnashy file from fastq file
#Last update: 07/02/2016 --updated _PHIX_
#Author:Jinchun Zhang
#All read should be named in this fashion: 
#	EXP_experimentname_BRD_barcodeseq_IDX_indexsequence_Read1.fastq
#	EXP_experimentname_BRD_barcodeseq_IDX_indexsequence_Read2.fastq
#	EXP_experimentname_BRD_barcodeseq_IDX_indexsequence_Index.fastq	
#For each experiment, output a gnashy txt file
#For each pair of index and barcode, output a Filter_Summary.csv file for Quality COntrol
#	Header : Aggregate, Count
#For version2,it only output a summary.txt for Filtering step & No plasmid mapping step
#Note for Barcode_Primer with length 3, the read1 seq should be 
#		3bp Primer Barcode+TTTACGCAGACTATCTTTCTAGGGTTAA+(TCTAGCTGCATCAGGATCATATCGT) or Genome_seq
#	for Barcode Primer with length 5, the read1 seq should be 
#		5bp Primer Barcode+GCGTCAATTTTACGCAGACTATCTTTCTAGGGTTAA+(TCTAGCTGCATCAGGATCATATCGT) or Genome_seq
########################################################################################

from Bio.Seq import Seq
from Bio import SeqIO
import numpy as np
import pandas as pd
import argparse
import os
import re
import subprocess
#import time

parser = argparse.ArgumentParser(prog='preCallingCards.py', description="Pre Calling Cards Wrap Script, find insertions from raw sequencing data ")
parser.add_argument("--i",dest="inp", type=str,nargs=1, help="type the Input directory, <../../INUPT/> ")
parser.add_argument("--o",dest="out", type=str,nargs=1, help="type the Output directory,  <../../OUTPUT/>")
parser.add_argument("--p",dest="mapping", type=str,nargs=1, help="type the mapping option, <SINGLE> ,<PAIR>")
parser.add_argument("--r",dest="refg", type=str,nargs=1, help="type the reference genome used to mapping,<hg19>,<hg18>,<mm10>,<mm9>")
parser.add_argument("--b",dest="brc", type=str,nargs=1, help="type the barcode file directory, it should be comma seperated with first row as the practical meaning and second row as the sequence of the barcode")
parser.add_argument("--e",dest="exp", type=str,nargs=1, help="The experiment group which corresponds to barcode and index group,<exp_name1, exp_name2,>")
parser.add_argument("--t",dest="tt", type=str,nargs=1, default=['TRUE'],help="The option to keep TTAA for read2, <TRUE> ,<FALSE>, default is TRUE")
parser.add_argument("--5",dest="five", type=str,nargs=1, help="type the length to be trimmed at 5'  ")
parser.add_argument("--3",dest="three", type=str,nargs=1, help="type the length to be trimmed at 3' ")


args=parser.parse_args()
IN_folder=args.inp[0]
OUT_folder=args.out[0]
mapping_opt=args.mapping[0]
bfile=args.brc[0]
exp=(args.exp[0]).split(",")
refg=args.refg[0]
TTAA= args.tt[0]
five=args.five[0]
three=args.three[0]


ref_folder= "/scratch/rmlab/ref/bowtie2_indexes/"
phix_file = "/scratch/rmlab/ref/phiX.fa"
fa_phix= SeqIO.parse(phix_file,"fasta")
phix_seq=fa_phix.next()
plasmid={
'uncut_plasmid_seq':"CCGGTACTCGAGATCCCCCAGCGGAAGAGCGCCCAATACG", #40bp uncut plasmid seq
'taq_plasmid_seq':"TTAACGATAAGTAAAATGTAAAATCACAGG", #30bp 
'csp61_plasmid_seq': "TTTACATATATATTTATTAGACAAGAAAAG", #30bp
'msp1_plasmid_seq': "AAAAAAGACCCGACGATATGATCCTGATGC", #30bp
}


#######################Change the Barcode file  when needed###############################
#read in csv barcode file and change into dictionary
#Input should be exp1, barcode11, index11,index12, index13,index14
#                exp2, barcode21, index21,index22, index23,index24
#Change the input to dictionary, experiment_name is key

def _READB_(filename):
        D={}
        f = open(filename,"r")
        for x in f:
                iarray=x.split("\n")[0].split(",")
                D[iarray[0]]=[iarray[1], iarray[2:]]
        return D

##############################FILTER############################################
#A comparing system between sequences that allowed *% percent mismatch
def _MATCH_ (X,Y,m):
        count=0
	ll=min(len(Y),len(X))
        for i in range(ll):
                if X[i]==Y[i]:
                        count=count+1
        if (ll-count) <= m:
                out=True
        else:
                out=False
        return out


#Chieck if every 25bp of the test sequence matched phix_sequence 
#X should be the test sequence and Y is the phix sequence.
def _PHIX_ (X,phix):
        count=0
        X=str(X)
        Y=str(phix.seq)
        ma=int( (len(X))/25 )+1
        out=False
        for i in range(1, ma):
                matchobj = re.search( X[(i-1)*25 : i*25] , Y)
                if matchobj:
                        count += 1
        if count >= 3:
                out=True
        else:
                phix=phix.reverse_complement()
                Y=str(phix.seq)
                count=0
                for i in range(1, ma):
                        matchobj = re.search( X[(i-1)*25: i*25] , Y)
                        if matchobj:
                                count += 1
                if count>=3:
                        out=True
        return out



# Note:  Read1 (PCR barcode + 24bp LTR seq(TTTACGCAGACTATCTTTCTAG) + TTAA(option) + genomic seq)
def _R1FILTER_(read,b):
        count={ 'Read1_total':0, 'Read1_unknown':0 , 'Read1_plasmid':0, 'Read1_genome_TTAA':0,
             'Read1_barcode':0,  'Read1_genome_Other':0, 'Read1_phix':0
            }
	ID=[]
	G_other=[]
	G_TTAA =[]
        bl=len(b)
        for rec in read:
		count['Read1_total'] += 1
                if str(rec.seq[0:bl]) == b:
			count['Read1_barcode'] += 1
                        test=str(rec.seq[bl:bl+53])
                        if _MATCH_( test , "TTTACGCAGACTATCTTTCTAGGGTTAATCTAGCTGCATCAGGATCATATCGT", 4 ):
                                count['Read1_plasmid'] += 1
                        elif _MATCH_( test[0:28] , 'TTTACGCAGACTATCTTTCTAGGGTTAA',2):
                                count['Read1_genome_TTAA'] += 1
				G_TTAA.append(1)
				G_other.append(0)
				ID.append(rec.id)
                        elif _MATCH_(str(test[0:24]) ,'TTTACGCAGACTATCTTTCTAGGG',1):
                                count['Read1_genome_Other'] += 1
				G_TTAA.append(0)
				G_other.append(1)
				ID.append(rec.id)
                elif _PHIX_(str(rec.seq), phix_seq):
                        count['Read1_phix'] += 1
                else:
                        count['Read1_unknown'] += 1
	out=pd.DataFrame ({"ID":ID, "Genome_Other_Index":G_other, "Genome_TTAA_Index":G_TTAA})
        return count,out

#Read2 has to be filter out from uncut plasmid sequence.
#return a dictionary with 4 keys have sequence id mathces plasmid sepquence.
def _R2FILTER_ (read,plm):
        count={'Read2_Total':0, 'Read2_uncut':0,'Read2_unknown':0,'Read2_phix':0,
             'Read2_taq':0,'Read2_csp61':0,'Read2_msp1':0,'Read2_plasmid_taq':0,
             'Read2_plasmid_csp61':0,'Read2_plasmid_msp1':0
            }
	ID=[]
	msp1=[]
	taq=[]
	csp61=[]
        for rec in read:
		count['Read2_Total'] += 1
                if _MATCH_(str(rec.seq[0:40]), str(plm['uncut_plasmid_seq']),4):
                        count['Read2_uncut'] += 1
                elif str(rec.seq[0:11])=='CCGGTACTCGA':
                        if _MATCH_(str(rec.seq[11:41]),str(plm['taq_plasmid_seq']),3):
                                count['Read2_plasmid_taq'] += 1
                        else:
                                count['Read2_taq'] += 1
				taq.append(1)
				msp1.append(0)
				csp61.append(0)
				ID.append(rec.id) 
                elif str(rec.seq[0:7])=='CCGGTAC':
                        if _MATCH_(str(rec.seq[7:37]),str(plm['csp61_plasmid_seq']),3):
                                count['Read2_plasmid_csp61'] += 1
                        else:
                                count['Read2_csp61'] += 1
				csp61.append(1)
				taq.append(0)
				msp1.append(0)
				ID.append(rec.id)
                elif str(rec.seq[0:5])=='CCGGC':
                        if  _MATCH_(str(rec.seq[5:35]), str(plm['msp1_plasmid_seq']),3):
                                count['Read2_plasmid_msp1'] += 1
                        else:
                                count['Read2_msp1'] += 1
				msp1.append(1)
				csp61.append(0)
				taq.append(0)
				ID.append(rec.id)
                elif _PHIX_(str(rec.seq),phix_seq):
                        count['Read2_phix'] += 1
                else:
                        count['Read2_unknown'] += 1
	out=pd.DataFrame({"ID":ID, "msp1_index":msp1, "taq_index":taq, "csp61_index":csp61})
        return count, out

###########################Index Read###################################
# For each experiment we have multiple index read
def _IndFILTER_ (read, iseq):
	count={'Index_total':0, "Index_match":0, "Index_notmatch":0}
	matchid= []
	for rec in read:
		count['Index_total'] += 1
		if str(rec.seq[0:len(iseq)]) == str(iseq):
			count["Index_match"] += 1
			matchid.append(rec.id)
		else:
			count["Index_notmatch"] += 1
		out = pd.DataFrame({'ID':matchid})
	return count, out

#####################FILTER Function for all the 
def _FILTER_(r1, r2, ir,o, tb, pb, p):
        if not os.path.exists(o):
                os.makedirs(o)
	summary_outfile = o+"Filter_Summary.csv"
        fastq_read1 = SeqIO.parse(r1,"fastq")
        fastq_read2 = SeqIO.parse(r2,"fastq")
	fastq_indexread = SeqIO.parse(ir,"fastq")
        R1_count, read1_out = _R1FILTER_(fastq_read1, pb)
        R2_count, read2_out = _R2FILTER_ (fastq_read2, plasmid)
	Ind_count, iread_out = _IndFILTER_ (fastq_indexread,tb )
	Barcode_Index_Matched = len(pd.merge(read1_out,iread_out,on='ID'))
	summary =  dict(R1_count.items()+ R2_count.items()+Ind_count.items()+[('Barcode_Index_Matched',Barcode_Index_Matched)])
        summary = pd.DataFrame(summary.items(),columns=['Aggregate','Count'])
	summary = summary.sort(['Aggregate'])
	summary.to_csv(summary_outfile, index=False, header=['Aggregate','Count'])
	if p == 'TRUE':
        	outfile1=o+"Filtered_TTAA_cutter_R1.fq"
                outfile2=o+"Filtered_TTAA_cutter_R2.fq"
		r1_out=read1_out[read1_out['Genome_TTAA_Index']==1]
		keep = pd.merge(pd.merge(r1_out,read2_out,on='ID'),iread_out,on ='ID')
	else:
		outfile1=o+"Filtered_genomic_cutter_R1.fq"
		outfile2=o+"Filtered_genomic_cutter_R2.fq"
                keep = pd.merge(pd.merge(read1_out,read2_out,on='ID'),iread_out,on ='ID')
        handle1 =  open(outfile1,'w')
        handle2 =  open(outfile2,'w')
        for rec in SeqIO.parse(r1, "fastq"):
		block = keep[keep['ID']==rec.id]
                if len(block) == 1:
                        rec = rec[24:] #trim before TTAA
                        print >> handle1, rec.format("fastq")
        for rec in SeqIO.parse(r2, "fastq"):
		block = keep[keep['ID'] == rec.id]
                if len(block) == 1:
                        if int(block['taq_index'])== 1:
				rec = rec[8:]
			elif int(block['csp61_index'])==1:
				rec = rec[4:]
                        print >> handle2, rec.format("fastq")

############################Mapping Function##########################################################
def _MAP_ (map_folder, ref, opt, r1, r2,fv,tr):
        if not os.path.exists(map_folder):
                os.makedirs(map_folder)
        sam_file=map_folder+"/Mapped.sam"
        if opt=="SINGLE":
		mapping_string="bowtie2 -x "+ref_folder+ref+"  "+r2+" -S "+ sam_file+" -5 "+ fv +"  -3 "+ tr + "1> "+map_folder+"/stdout.txt "+ "2> "+map_folder+"/stderr_all.txt"
        elif opt=="PAIR":
		mapping_string="bowtie2 -I 0 -X 1000 --dovetail -x "+ ref_folder+ref+" -1 "+ r1+" -2 "+r2+" -S "+ sam_file+" -5 "+ fv +"  -3 "+ tr + "1> "+map_folder+"/stdout.txt "+ "2> "+map_folder+"/stderr_all.txt"
        print " "+ mapping_string
        os.system(mapping_string)
	os.system('grep -v Warning '+map_folder+'/stderr_all.txt > '+map_folder+'/stderr.txt')
        subprocess.call("samtools view "+sam_file+" -S -b -q 22 -o "+map_folder+"/Mapped.bam" ,shell=True) ##-q 10 (mapping quality cut off=10)
        subprocess.call("samtools sort "+map_folder+"/Mapped.bam  -o "+ map_folder+"/Mapped.sorted.bam" ,shell=True)
        subprocess.call("bedtools bamtobed -i "+map_folder+"/Mapped.sorted.bam  > "+ map_folder+"/Mapped.sorted.bed",shell=True)

##########################Used to turn mapping bed to gnashy file#################
def _GNASHY_ (inp):
        inbed=pd.read_csv(inp, sep='\t',header=-1)
        del inbed[4],inbed[5]
        X=inbed.groupby([0,1,2])[3].count().reset_index()
        X[4]=0
        for m in range(len(X)):
                X.loc[m,4]=int(0.5*(X.loc[m,2]+X.loc[m,1]))
                matchobj=re.search("(?<=chr)\d+", X.loc[m,0])
                if matchobj:
                        X.loc[m,0]=matchobj.group(0)
                elif    X.loc[m,0]=='chrX':
                        X.loc[m,0]=23
                elif X.loc[m,0]=='chrY':
                        X.loc[m,0]=24
                else:
                        X.loc[m,0]='DEL'
        total=pd.DataFrame(columns=["Chromosome","Position","Count"])
        total["Chromosome"]=X[0]
        total["Position"]=X[4]
        total["Count"]=X[3]
        total= total[total["Chromosome"] !='DEL']
        return total

###########################################################################################
barcode_dic= _READB_(bfile)
if not os.path.exists(OUT_folder):
        os.makedirs(OUT_folder)

#Loop for each experiment
for e in exp:
	a=0
        barcode_seq=barcode_dic[e][0]
        index_group=barcode_dic[e][1]
        e_output=OUT_folder+"Experiment_"+str(e)+"/"  #One folder for one experiment 
        if not os.path.exists(e_output):
                os.makedirs(e_output)
	gnashy_outfile=e_output+'Experiment_'+str(e)+".gnashy"
        #For each pair of index and barcode, filter reads and mapping back to reference group
	print 'Start analyze experiment: '+e #+' at '+time.strftime('%X %x %Z')
        for i in index_group:
		a=a+1
                prefix="EXP_"+e+"_BRD_"+barcode_seq+"_IDX_"+i ####Specify the prefix for a barcode*index combination
		print
		print  prefix
                read1=IN_folder+prefix+"_Read1.fastq"
                read2=IN_folder+prefix+"_Read2.fastq"
		iread=IN_folder+prefix+"_Index.fastq"
                Ind_Folder=e_output+prefix+"/"     #one folder for one index
                if not os.path.exists(Ind_Folder):
                        os.makedirs(Ind_Folder)
                #Step1 fileter read1 and read2.fq
#                print "		_FILTER_ step starts " #+time.strftime('%X %x %Z')
#                _FILTER_(read1, read2, iread, Ind_Folder+"FILTER_OUTPUT/", i, barcode_seq, TTAA)
                #Step2 MAPPINGQ
                mapout=Ind_Folder+"MAP_OUT"
                print "		_MAP_ Step starts " # +time.strftime('%X %x %Z')
                if TTAA == 'FALSE':
                        _MAP_ (mapout, refg, mapping_opt, Ind_Folder+"FILTER_OUTPUT/Filtered_genomic_cutter_R1.fq", Ind_Folder+"FILTER_OUTPUT/Filtered_genomic_cutter_R2.fq", five, three)
                else :
                        _MAP_(mapout, refg, mapping_opt, Ind_Folder+"FILTER_OUTPUT/Filtered_TTAA_cutter_R1.fq", Ind_Folder+"FILTER_OUTPUT/Filtered_TTAA_cutter_R2.fq", five, three)
		#Step3  merge and output gnashy file for each barcode+index_group pair
                print "		_GNASHY_ step starts "  #+ time.strftime('%X %x %Z')
                gnashyout=_GNASHY_(mapout+"/Mapped.sorted.bed")
                if a ==1:
                        gnashyout.to_csv(gnashy_outfile,index=False,sep='\t',header=["Chromosome","Coordinate","Count"])
                else:
                        gnashyout.to_csv(gnashy_outfile,index=False,sep='\t',header=False,mode="a")

