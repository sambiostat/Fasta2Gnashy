#!usr/bin/python

############################################
#This is the QC scripts for the calling cards experiment which is used to summary the QC information of the step produce Gnahsy from Sequencing 
#Run for each experiment
#Input: path of the Outputfolder from Fasta2Gnashy.py
#different from Version1: it used when bowtie2 were called during FASTA 2 GNASHY step. Since the stderr of mapping is different between BOWTIE and BOWTIE2
#LAST UPDATED: 2016-07-10
###########################################

import numpy as np
import pandas as pd
import argparse
import sys
import re
import subprocess

parser = argparse.ArgumentParser(prog='preCallingCards.py', description="Pre Calling Cards Wrap Script, find insertions from raw sequencing data ")
parser.add_argument("--i",dest="inp", type=str,nargs=1, help="type the Input directory, <../../INUPT/> ")
parser.add_argument("--o",dest="out", type=str,nargs=1, help="type the Output File Name,  <../../QC_summary.out>")
parser.add_argument("--b",dest="brc", type=str,nargs=1, help="type the barcode file directory, it should be comma seperated with first row as the practical meaning and second row as the sequence of the barcode")
parser.add_argument("--e",dest="exp", type=str,nargs=1, help="The experiment of interest,only one is allowed")

args=parser.parse_args()
IN_folder=args.inp[0]
outfilename=args.out[0]
bfile=args.brc[0]
exp=args.exp[0]

def _READB_(filename):
        D=[]
        f = open(filename,"r")
        for x in f:
                iarray=x.split("\n")[0].split(",")
		if str(iarray[0]) == exp:
			ind=iarray[2:]
			for i in ind:
                		D.append ( [iarray[0],iarray[1],i])
	D=pd.DataFrame(D)
	D.columns = ['Experiment','Primer_Barcode','Transposon_Barcode']
        D['Transposon_Reads'] = pd.Series(np.zeros(len(D)))
	D['Fraction_Transposon_Reads'] = pd.Series(np.zeros(len(D)))
	D['Reads_Processed'] = pd.Series(np.zeros(len(D)))
	D['Percent_Aligned'] = pd.Series(np.zeros(len(D)))
	D['Percent_Unaligned'] = pd.Series(np.zeros(len(D)))
        return D

def Analyze_Gnashyfile(gnashyfilename, fhandle):
        fhandle.write( "Experiment analyzed: " + gnashyfilename)
        temp_frame = pd.read_csv(gnashyfilename,delimiter = "\t")
        fhandle.write(" 	Number of transpostions: " + str(len(temp_frame))+"\n")
        fhandle.write(" 	Median number of reads: " + str(np.median(temp_frame['Count']))+"\n")
        num_insertions_lt_2 = float(len(temp_frame[temp_frame['Count']<2]))/len(temp_frame)
        fhandle.write(" 	Transpositions with only 1 read: " + str(num_insertions_lt_2)+"\n")
        num_insertions_lt_6 = float(len(temp_frame[temp_frame['Count']<6]))/len(temp_frame)
        fhandle.write( "         Transpositions with 5 or fewer reads: " + str(num_insertions_lt_6)+"\n")
        num_insertions_lt_11 = float(len(temp_frame[temp_frame['Count']<11]))/len(temp_frame)
        fhandle.write( "         Transpositions with 10 or fewer reads: " + str(num_insertions_lt_11)+"\n")

qc_frame = _READB_(bfile)

for i in qc_frame.index:
	e=qc_frame.ix[i,'Experiment']
	b=qc_frame.ix[i,'Primer_Barcode']
	t=qc_frame.ix[i,'Transposon_Barcode']
	in_name=IN_folder+'Experiment_'+e+'/EXP_'+e+'_BRD_'+b+'_IDX_'+t+'/FILTER_OUTPUT/Filter_Summary.csv'
	infile=pd.read_csv(in_name)
	qc_frame.ix[i, 'Transposon_Reads']=int(infile[infile['Aggregate'] == 'Barcode_Index_Matched']['Count']) 
	in_name=IN_folder+'Experiment_'+e+'/EXP_'+e+'_BRD_'+b+'_IDX_'+t+'/MAP_OUT/stderr.txt'
	infile = open(in_name,'r')
	text = infile.readline()
	matchObj = re.search(r'(\d*) reads',text)
	reads_processed = int(matchObj.group(1))
	text = infile.readline()
	text = infile.readline()
	matchObj = re.search(r'(\d+)\s.(.*)%\)* aligned 0 times',text)
        Percent_aligned = float(matchObj.group(2))
        Percent_not_aligned = float( 100 - Percent_aligned )
        qc_frame.ix[i,'Reads_Processed'] = reads_processed
        qc_frame.ix[i,'Percent_Aligned'] = Percent_aligned
        qc_frame.ix[i,'Percent_Unaligned'] = Percent_not_aligned

for i in qc_frame.index:
	qc_frame.ix[i,'Fraction_Transposon_Reads'] = qc_frame.ix[i,'Transposon_Reads']/(sum(qc_frame['Transposon_Reads']))

qc_frame.to_csv(outfilename,sep = "\t")

filehandle = open(outfilename,"a")
filehandle.write("Mapping summary\n") 
filehandle.write('Experience: ' + exp+"\n")
filehandle.write ("Total fraction of reddads assigned to this experiment: " + str(sum(qc_frame['Fraction_Transposon_Reads']))+'\n' )
in_name=IN_folder+'Experiment_'+e+'/Experiment_'+e+'.gnashy'
Analyze_Gnashyfile(in_name,filehandle)
filehandle.close()
