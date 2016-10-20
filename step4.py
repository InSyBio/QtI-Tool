#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################
#This code is the outcome of the project Innovative Processing of TMT Proteomics. This is a common project of InSyBio LTD and Nestle Institute of Health Sciences. 
#This program implements the step 2 code
####################################

from pyopenms import *
import os, shutil
import random
import math
import copy
import time
import sys
import subprocess
import numpy
import scipy.stats as st
from shutil import copyfile
import platform
import scipy
from collections import defaultdict
from bisect import bisect_left
from subprocess import Popen
from svmutil import *

def generate_final_report(unified_spectrum_list_filename,mapping_spectra_filename,quantitative_values_filename,classes_filename, similarities_mapping_spectra_filename,filtered_files_folder,quantified_search_folder,unified_spectra_search_filename,similarity_threshold,selected_spectra_filename,scaffold_peptide_report_filename,quantified_folder,mode,output_folder,replicates,missing_values_filename,use_inputs_vector_filename):

	unified_spectra = MSExperiment()
	file = MzMLFile()
	file.load(unified_spectrum_list_filename, unified_spectra)
	print("Unified spectra list was successfully parsed!")
	#parse use_inputs_vector_filename
	use_inputs_vector_fid=open(use_inputs_vector_filename,'r')
	use_inputs_vector=list()
	for line1 in use_inputs_vector_fid:
		words=line1.split('\t')
		for i in range(len(words)):
			use_inputs_vector.append(words[i].strip())
	number_of_tmt_channels_used=use_inputs_vector.count("1")
	missing_values_fid2      = open(missing_values_filename,'r')
	missing_values_list=list()
	missing_values_raws=0;
	for line1 in missing_values_fid2:
		words=line1.split('\t')
		missing_values_list.append([])
		for i in range(len(words)):
			missing_values_list[missing_values_raws].append(words[i].strip())
		missing_values_raws+=1
	#parse mapping spectra filename
	mapping_spectra=list()
	filenames=list();
	mapping_fid=open(mapping_spectra_filename,"r")
	number_of_lines=0
	for line1 in mapping_fid:
		if number_of_lines==0:
			words=line1.split("\t")
			for i in range(len(words)):
				filenames.append(words[i].strip())
		else:
			mapping_spectra.append([])
			words=line1.split("\t")
			for i in range(len(words)):
				mapping_spectra[number_of_lines-1].append(words[i].strip())
		number_of_lines=number_of_lines+1
	print("Mapping spectra file was successfully parsed!")
	#parse scaffold report
	scaffold_identified_proteins=list()
	scaffold_peptide_report_fid=open(scaffold_peptide_report_filename,"r")
	number_of_lines=0
	rows_so_far=0
	for line1 in scaffold_peptide_report_fid:
		number_of_lines=number_of_lines+1;
	total_number_of_lines=number_of_lines-1;	
	scaffold_peptide_report_fid=open(scaffold_peptide_report_filename,"r")
	number_of_lines=0
	rows_so_far=0
	for line1 in scaffold_peptide_report_fid:
		if number_of_lines>40 and number_of_lines<total_number_of_lines:
			words=line1.split("\t")
			scaffold_identified_proteins.append([])
			scaffold_identified_proteins[rows_so_far].append(words[2])
			scaffold_identified_proteins[rows_so_far].append(words[6].upper())
			rows_so_far=rows_so_far+1
		number_of_lines=number_of_lines+1
	print("Scaffold peptide report was successfully parsed!")
	#parse quantitative filename	
	quantitative_spectra=list()
	quantitative_fid=open(quantitative_values_filename,"r")
	number_of_lines1=0
	for line1 in quantitative_fid:
		number_of_lines1=number_of_lines1+1
	number_of_lines=0
	
	quantitative_fid=open(quantitative_values_filename,"r")
	for line1 in quantitative_fid:
		if number_of_lines!=(number_of_lines1-1):
			quantitative_spectra.append([])
			words=line1.split("\t")
			for i in range(len(words)-1):
				if use_inputs_vector[i%6]=="1":
					quantitative_spectra[number_of_lines].append(float(words[i]))
		number_of_lines=number_of_lines+1
	print("Quantitative file was successfully parsed!")
	#parse similarities filename
	similarities=list()
	similarities_fid=open(similarities_mapping_spectra_filename,"r")
	number_of_lines=0
	for line1 in similarities_fid:
		if number_of_lines!=0:
			similarities.append([])
			words=line1.split("\t")
			for i in range(len(words)):
				similarities[number_of_lines-1].append(float(words[i]))
		number_of_lines=number_of_lines+1
	print("Similarities file was successfully parsed!")
	#parse unified spectra search filename
	unified_spectra_search_info=list()
	unified_spectra_search_fid=open(unified_spectra_search_filename,"r")
	number_of_lines=0
	queries_so_far=0
	queries_to_spectra=list()
	flag=0;
	for line1 in unified_spectra_search_fid:
		if "decoy_peptides" in line1:
			flag=1;
		string_to_search="q"+str(queries_so_far+1)+"_p1="
		if line1.startswith("scans="):
			data1=line1.split("=")
			queries_to_spectra.append(int(data1[1].strip()))
		if flag==0 and line1.startswith(string_to_search):
			two_parts=line1.split("=")
			four_parts=two_parts[1].split(";")
			if two_parts[1].strip()=="-1":
				unified_spectra_search_info.append([]);
				unified_spectra_search_info[queries_so_far].append(-1)
				unified_spectra_search_info[queries_so_far].append(-1)
				unified_spectra_search_info[queries_so_far].append(-1)
				unified_spectra_search_info[queries_so_far].append(-1)
			else:	
					
				words=four_parts[0].split(",")	
				unified_spectra_search_info.append([]);
				unified_spectra_search_info[queries_so_far].append(words[4])
				unified_spectra_search_info[queries_so_far].append(words[7].replace(".",","))
				four_parts[1].strip()
				different_proteins_full=four_parts[1].strip().split(",")
				different_proteins="";
				peptide_positions="";
				number_of_proteins_so_far=0;
				for protein in different_proteins_full:
					protein_info=protein.split(":")
					pept_pos=""
					for pos in xrange(1,len(protein_info)):
						if pos==len(protein_info)-1:
							pept_pos=pept_pos+protein_info[pos]
						else:
							pept_pos=pept_pos+protein_info[pos]+":"
					if number_of_proteins_so_far==0:
						different_proteins=different_proteins+protein_info[0].replace("\"","")
						peptide_positions=peptide_positions+pept_pos
					else:
						different_proteins=different_proteins+","+protein_info[0].replace("\"","")
						peptide_positions=peptide_positions+","+pept_pos
					number_of_proteins_so_far=number_of_proteins_so_far+1;
				#unified_spectra_search_info[queries_so_far].append(four_parts[1].strip())
				unified_spectra_search_info[queries_so_far].append(different_proteins)
				unified_spectra_search_info[queries_so_far].append(peptide_positions)
			queries_so_far=queries_so_far+1;
	print("Unified spectra search file was successfully parsed!")					

	#parse selected spectra filename
	selected_spectra=list()	
	selected_spectra_fid=open(selected_spectra_filename,"r")
	number_of_lines=0
	for line1 in selected_spectra_fid:
		words=line1.split("\t")
		selected_spectra.append(words[0])
		number_of_lines=number_of_lines+1;
	print("Selected spectra file was successfully parsed!")
	peptides=list()
	for pep in range(len(scaffold_identified_proteins)):
		peptides.append(scaffold_identified_proteins[pep][1])
	final_report_fid=open(output_folder+"final_report_"+str(time.strftime("%Y_%m_%d"))+".tsv","w")
	final_report_fid.write("Spectrum Query Number\tQueries Peak List\tHighest Mascot Searched Peptide\tMascot Score\tAssociated Proteins from Mascot\tPositions of peptides in associated Proteins from Mascot\tAssociated proteins from Scaffold\t")
	for i in range(len(filenames)):
		for j in range(6):
			if use_inputs_vector[j]=="1":
				final_report_fid.write("Quantitative Value of "+ str(filenames[i])+" TMT channel "+str(j)+"\t")
	final_report_fid.write("Differentially Quantified in more than 50% of experiments Flag\n")	
	common_spectra_experiments_1_2=0
	common_spectra_experiments_1_3=0
	common_spectra_experiments_2_3=0
	common_spectra_experiments_1_4=0
	common_spectra_experiments_2_4=0
	common_spectra_experiments_3_4=0
	common_spectra_experiments_1_2_4=0
	common_spectra_experiments_1_3_4=0
	common_spectra_experiments_2_3_4=0
	common_spectra_experiments_1_2_3=0
	common_spectra_experiments_1_2_3_4=0
	spectra_exp1=0
	spectra_exp2=0
	spectra_exp3=0
	spectra_exp4=0
	common_peptides_experiments_1_2=0
	common_peptides_experiments_1_3=0
	common_peptides_experiments_2_3=0
	common_peptides_experiments_1_4=0
	common_peptides_experiments_2_4=0
	common_peptides_experiments_3_4=0
	common_peptides_experiments_1_2_4=0
	common_peptides_experiments_1_3_4=0
	common_peptides_experiments_2_3_4=0
	common_peptides_experiments_1_2_3=0
	common_peptides_experiments_1_2_3_4=0
	peptides_exp1=0
	peptides_exp2=0
	peptides_exp3=0
	peptides_exp4=0;
	common_peptides=list()
	common_peptides_50=list()
	common_peptides_90=list()
	peptides1=list()
	peptides2=list()
	peptides3=list()
	peptides4=list()
	peptides1_2=list()
	peptides2_3=list()
	peptides1_3=list()
	peptides1_4=list()
	peptides2_4=list()
	peptides3_4=list()
	peptides1_2_3=list()
	peptides1_2_3_4=list()
	peptides2_3_4=list()
	peptides1_2_4=list()
	peptides1_3_4=list()
	common_spectra_50_percent=0
	common_spectra_90_percent=0
	common_identified_spectra=0
	common_identified_spectra_unique=0
	dif_quant_spectra_50_percent=0
	dif_quant_spectra_90_percent=0
	dif_quant_peptides_90_percent=list()
	dif_quant_identified_spectra=0
	dif_quant_peptides=list()
	dif_quant_identified_spectra_unique=0

	for i in range(queries_so_far):
		spectra_to_files=list();
		for j in range(len(filenames)):
			spectra_to_files.append(0);
		spectra_to_experiments=list();
		for j in range(len(filenames)/replicates):
			spectra_to_experiments.append(0);
		final_report_fid.write(str(i)+"\t");
		peaks=unified_spectra[i].get_peaks();
		if len(peaks)==0:
			string_to_write="[]";
		else:
			mz_values_old, intensities_old = zip(*peaks)
			mz_values = [float(v) for v in mz_values_old]
			intensities = [float(v) for v in intensities_old]
			string_to_write="["
			for k in range(len(intensities)):
				string_to_write=string_to_write+"("+str(mz_values[k])+":"+str(intensities[k])+")"
			string_to_write=string_to_write+"]"
		final_report_fid.write(string_to_write+"\t")
		try:		
			ind=queries_to_spectra.index(i)
		except ValueError:
			ind=-1;
		if ind!=-1:
			final_report_fid.write(str(unified_spectra_search_info[ind][0])+"\t")
			final_report_fid.write(str(unified_spectra_search_info[ind][1])+"\t")
			final_report_fid.write(str(unified_spectra_search_info[ind][2])+"\t")
			final_report_fid.write(str(unified_spectra_search_info[ind][3])+"\t")
			#write scaffold results
			if (str(unified_spectra_search_info[ind][0])=="-1"):
				final_report_fid.write("-1"+"\t")
			else:
				try:
					index = peptides.index(str(unified_spectra_search_info[ind][0]))
				except ValueError:
					index = -1			
				if index==-1:	
					final_report_fid.write("-1"+"\t")
				else:
					final_report_fid.write(str(scaffold_identified_proteins[index][0])+"\t")
		else:
			final_report_fid.write("-1"+"\t")
			final_report_fid.write("-1"+"\t")
			final_report_fid.write("-1"+"\t")
			final_report_fid.write("-1"+"\t")
			final_report_fid.write("-1"+"\t")
		for j in range(len(filenames)):
			if quantitative_spectra[i][j*6]!=-1:
				spectra_to_files[j]=1;
			for k in range(number_of_tmt_channels_used):
				string_to_write=str(quantitative_spectra[i][j*number_of_tmt_channels_used+k])
				string_to_write=string_to_write.replace(".",",")
				final_report_fid.write(string_to_write+"\t")
		for j in range(len(filenames)/replicates):
			for k in range(replicates):
				if spectra_to_files[j*replicates+k]==1:
					spectra_to_experiments[j]=1;
		sum_spectra_to_experiments=sum(spectra_to_experiments)
		if sum_spectra_to_experiments>=(0.5)*(len(filenames)/replicates):
			common_spectra_50_percent=common_spectra_50_percent+1;
		if sum_spectra_to_experiments>=(0.9)*(len(filenames)/replicates):
			common_spectra_90_percent=common_spectra_90_percent+1;
		if sum_spectra_to_experiments>0 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
			common_identified_spectra=common_identified_spectra+1;
			common_peptides.append(unified_spectra_search_info[ind][0]);
		if sum_spectra_to_experiments>=(0.5)*(len(filenames)/replicates) and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
			common_peptides_50.append(unified_spectra_search_info[ind][0]);
		if sum_spectra_to_experiments>=(0.9)*(len(filenames)/replicates) and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
			common_peptides_90.append(unified_spectra_search_info[ind][0]);
		if len(filenames)/replicates==3:
			if spectra_to_experiments[0]==1:
				spectra_exp1=spectra_exp1+1
			if spectra_to_experiments[1]==1:
				spectra_exp2=spectra_exp2+1
			if spectra_to_experiments[2]==1:
				spectra_exp3=spectra_exp3+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1:
				common_spectra_experiments_1_2=common_spectra_experiments_1_2+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1:
				common_spectra_experiments_1_3=common_spectra_experiments_1_3+1
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1:
				common_spectra_experiments_2_3=common_spectra_experiments_2_3+1
			if spectra_to_experiments[0]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp1=peptides_exp1+1
				peptides1.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp2=peptides_exp2+1
				peptides2.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp3=peptides_exp3+1
				peptides3.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2=common_peptides_experiments_1_2+1
				peptides1_2.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_3=common_peptides_experiments_1_3+1
				peptides1_3.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_2_3=common_peptides_experiments_2_3+1
				peptides2_3.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2_3=common_peptides_experiments_1_2_3+1
				peptides1_2_3.append(unified_spectra_search_info[ind][0])
		if len(filenames)/replicates==4:
			if spectra_to_experiments[0]==1:
				spectra_exp1=spectra_exp1+1
			if spectra_to_experiments[1]==1:
				spectra_exp2=spectra_exp2+1
			if spectra_to_experiments[2]==1:
				spectra_exp3=spectra_exp3+1
			if spectra_to_experiments[3]==1:
				spectra_exp4=spectra_exp4+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1:
				common_spectra_experiments_1_2=common_spectra_experiments_1_2+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1:
				common_spectra_experiments_1_3=common_spectra_experiments_1_3+1
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1:
				common_spectra_experiments_2_3=common_spectra_experiments_2_3+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_1_4=common_spectra_experiments_1_4+1
			if spectra_to_experiments[1]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_2_4=common_spectra_experiments_2_4+1
			if spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_3_4=common_spectra_experiments_3_4+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_1_2_4=common_spectra_experiments_1_2_4+1
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_1_3_4=common_spectra_experiments_1_3_4+1
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1:
				common_spectra_experiments_2_3_4=common_spectra_experiments_2_3_4+1	
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1:
				common_spectra_experiments_1_2_3=common_spectra_experiments_1_2_3+1	
				
			if spectra_to_experiments[0]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp1=peptides_exp1+1
				peptides1.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp2=peptides_exp2+1
				peptides2.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp3=peptides_exp3+1
				peptides3.append(unified_spectra_search_info[ind][0])
				
			if spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				peptides_exp4=peptides_exp4+1
				peptides4.append(unified_spectra_search_info[ind][0])
			
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2=common_peptides_experiments_1_2+1
				peptides1_2.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_3=common_peptides_experiments_1_3+1
				peptides1_3.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_2_3=common_peptides_experiments_2_3+1
				peptides2_3.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2_3=common_peptides_experiments_1_2_3+1
				peptides1_2_3.append(unified_spectra_search_info[ind][0])
				
			if spectra_to_experiments[0]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_4=common_peptides_experiments_1_4+1
				peptides1_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_2_4=common_peptides_experiments_2_4+1
				peptides2_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_3_4=common_peptides_experiments_3_4+1
				peptides3_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2_4=common_peptides_experiments_1_2_4+1
				peptides1_2_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_3_4=common_peptides_experiments_1_3_4+1
				peptides1_3_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_2_3_4=common_peptides_experiments_2_3_4+1
				peptides2_3_4.append(unified_spectra_search_info[ind][0])
			if spectra_to_experiments[0]==1 and spectra_to_experiments[1]==1 and spectra_to_experiments[2]==1 and spectra_to_experiments[3]==1 and (unified_spectra_search_info[ind][0]!="-1" and unified_spectra_search_info[ind][0]!=-1):
				common_peptides_experiments_1_2_3_4=common_peptides_experiments_1_2_3_4+1
				peptides1_2_3_4.append(unified_spectra_search_info[ind][0])
			
		try:
			index = selected_spectra.index(str(i))
		except ValueError:
			index = -1
		if index==-1:	
			final_report_fid.write("-1"+"\n")
		else:
			final_report_fid.write("1"+"\n")
		if index!=-1 and sum_spectra_to_experiments>=(0.5)*(len(filenames)/replicates):
			dif_quant_spectra_50_percent=dif_quant_spectra_50_percent+1			
		if index!=-1 and sum_spectra_to_experiments>=(0.9)*(len(filenames)/replicates):
			dif_quant_spectra_90_percent=dif_quant_spectra_90_percent+1
		if index!=-1 and unified_spectra_search_info[ind][0]!=-1:
			dif_quant_identified_spectra=dif_quant_identified_spectra+1;
			dif_quant_peptides.append(unified_spectra_search_info[ind][0])
		if index!=-1 and unified_spectra_search_info[ind][0]!=-1 and sum_spectra_to_experiments>=(0.9)*(len(filenames)/replicates):
			dif_quant_peptides_90_percent.append(unified_spectra_search_info[ind][0])
	print("Final report file was successfully created!")
	result_metrics_fid = open(output_folder+"result_metrics_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	common_identified_spectra_unique=len(list(set(common_peptides)))
	common_identified_spectra_unique_50=len(list(set(common_peptides_50)))
	common_identified_spectra_unique_90=len(list(set(common_peptides_90)))
	result_metrics_fid
	print("Number of Common Spectra in more than 50% of experiments:"+ str(common_spectra_50_percent) )
	print("Number of Common Spectra in more than 90% of experiments:"+ str(common_spectra_90_percent) )
	print("Number of Common Peptides in more than 50% of experiments:"+ str(len(common_peptides_50)) )
	print("Number of Common Peptides in more than 90% of experiments:"+ str(len(common_peptides_90)) )
	#print("Number of Common Peptides in more than 90% of experiments:"+ str(common_spectra_90_percent) )
	print("Number of common unique identified spectra in more than 50% of the experiments:"+str(common_identified_spectra_unique_50))
	print("Number of common unique identified spectra in more than 90% of the experiments:"+str(common_identified_spectra_unique_90))
	
	print("Number of Common Identified Spectra:"+ str(common_identified_spectra) )
	print("Number of Unique Common Identified Spectra:"+ str(common_identified_spectra_unique) )
	print("Number of Differentially Quantified Spectra in more than 50% of experiments:"+ str(dif_quant_spectra_50_percent) )
	print("Number of Differentially Quantified Spectra in more than 90% of experiments:"+ str(dif_quant_spectra_90_percent) )
	print("Number of Differentially Quantified Identified Spectra:"+ str(dif_quant_identified_spectra) )
	print("Number of Unique Differentially Quantified Identified Spectra:"+ str(len(list(set(dif_quant_peptides)))) )
	print("Number of Unique Differentially Quantified Identified Spectra in more than 90% of experiments:"+ str(len(list(set(dif_quant_peptides_90_percent)))) )
	result_metrics_fid.write("Number of Common Spectra in more than 50% of experiments:"+ str(common_spectra_50_percent) )
	result_metrics_fid.write("\n")
	result_metrics_fid.write("Number of Common Spectra in more than 90% of experiments:"+ str(common_spectra_90_percent) )
	result_metrics_fid.write("\n")
	result_metrics_fid.write("Number of Common Identified Spectra:"+ str(common_identified_spectra) )
	result_metrics_fid.write("\n")
	result_metrics_fid.write("Number of Differentially Quantified Spectra in more than 50% of experiments:"+ str(dif_quant_spectra_50_percent) )
	result_metrics_fid.write("\n")
	result_metrics_fid.write("Number of Differentially Quantified Spectra in more than 90% of experiments:"+ str(dif_quant_spectra_90_percent) )
	result_metrics_fid.write("\n")
	result_metrics_fid.write("Number of Differentially Quantified Identified Spectra:"+ str(dif_quant_identified_spectra) )
	result_metrics_fid.write("Number of Unique Differentially Quantified Identified Spectra:"+ str(len(list(set(dif_quant_peptides)))) )
	result_metrics_fid.write("\n")
	if len(filenames)/replicates==3:
		print("Number of spectra in experiment 1:"+str(spectra_exp1))
		print("Number of spectra in experiment 2:"+str(spectra_exp2))
		print("Number of spectra in experiment 3:"+str(spectra_exp3))
		print("Number of common spectra in experiments 1 and 2:"+str(common_spectra_experiments_1_2))
		print("Number of common spectra in experiments 1 and 3:"+str(common_spectra_experiments_1_3))
		print("Number of common spectra in experiments 2 and 3:"+str(common_spectra_experiments_2_3))
		print("Number of common spectra in experiments 1,2 and 3:"+str(common_spectra_90_percent))
		print("Number of peptides in experiment 1:"+str(peptides_exp1))
		print("Number of peptides in experiment 2:"+str(peptides_exp2))
		print("Number of peptides in experiment 3:"+str(peptides_exp3))
		print("Number of common peptides in experiments 1 and 2:"+str(common_peptides_experiments_1_2))
		print("Number of common peptides in experiments 1 and 3:"+str(common_peptides_experiments_1_3))
		print("Number of common peptides in experiments 2 and 3:"+str(common_peptides_experiments_2_3))
		print("Number of common peptides in experiments 1,2 and 3:"+str(common_peptides_experiments_1_2_3))
		print("Number of unique peptides in experiment 1:"+str(len(list(set(peptides1)))))
		print("Number of unique peptides in experiment 2:"+str(len(list(set(peptides2)))))
		print("Number of unique peptides in experiment 3:"+str(len(list(set(peptides3)))))
		print("Number of unique common peptides in experiments 1 and 2:"+str(len(list(set(peptides1_2)))))
		print("Number of unique common peptides in experiments 1 and 3:"+str(len(list(set(peptides1_3)))))
		print("Number of unique common peptides in experiments 2 and 3:"+str(len(list(set(peptides2_3)))))
		print("Number of unique common peptides in experiments 1,2 and 3:"+str(len(list(set(peptides1_2_3)))))
		
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 1:"+str(spectra_exp1))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 2:"+str(spectra_exp2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 3:"+str(spectra_exp3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1 and 2:"+str(common_spectra_experiments_1_2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1 and 3:"+str(common_spectra_experiments_1_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 2 and 3:"+str(common_spectra_experiments_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1,2 and 3:"+str(common_spectra_90_percent))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 1:"+str(peptides_exp1))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 2:"+str(peptides_exp2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 3:"+str(peptides_exp3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1 and 2:"+str(common_peptides_experiments_1_2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1 and 3:"+str(common_peptides_experiments_1_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 2 and 3:"+str(common_peptides_experiments_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1,2 and 3:"+str(common_peptides_experiments_1_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 1:"+str(len(list(set(peptides1)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 2:"+str(len(list(set(peptides2)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 3:"+str(len(list(set(peptides3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1 and 2:"+str(len(list(set(peptides1_2)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1 and 3:"+str(len(list(set(peptides1_3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 2 and 3:"+str(len(list(set(peptides2_3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1,2 and 3:"+str(len(list(set(peptides1_2_3)))))
		result_metrics_fid.write("\n")
	if len(filenames)/replicates==4:	
		print("Number of spectra in experiment 1:"+str(spectra_exp1))
		print("Number of spectra in experiment 2:"+str(spectra_exp2))
		print("Number of spectra in experiment 3:"+str(spectra_exp3))
		print("Number of spectra in experiment 3:"+str(spectra_exp4))
		print("Number of common spectra in experiments 1 and 2:"+str(common_spectra_experiments_1_2))
		print("Number of common spectra in experiments 1 and 3:"+str(common_spectra_experiments_1_3))
		print("Number of common spectra in experiments 2 and 3:"+str(common_spectra_experiments_2_3))
		print("Number of common spectra in experiments 1 and 4:"+str(common_spectra_experiments_1_4))
		print("Number of common spectra in experiments 2 and 4:"+str(common_spectra_experiments_2_4))
		print("Number of common spectra in experiments 3 and 4:"+str(common_spectra_experiments_3_4))
		print("Number of common spectra in experiments 1,2 and 4:"+str(common_spectra_experiments_1_2_4))
		print("Number of common spectra in experiments 1,2 and 3:"+str(common_spectra_experiments_1_2_3))
		print("Number of common spectra in experiments 1,3 and 4:"+str(common_spectra_experiments_1_3_4))
		print("Number of common spectra in experiments 2,3 and 4:"+str(common_spectra_experiments_2_3_4))
		print("Number of common spectra in experiments 1,2,3 and 4:"+str(common_spectra_90_percent))
		print("Number of peptides in experiment 1:"+str(peptides_exp1))
		print("Number of peptides in experiment 2:"+str(peptides_exp2))
		print("Number of peptides in experiment 3:"+str(peptides_exp3))
		print("Number of peptides in experiment 3:"+str(peptides_exp4))
		print("Number of common peptides in experiments 1 and 2:"+str(common_peptides_experiments_1_2))
		print("Number of common peptides in experiments 1 and 3:"+str(common_peptides_experiments_1_3))
		print("Number of common peptides in experiments 2 and 3:"+str(common_peptides_experiments_2_3))
		print("Number of common peptides in experiments 1 and 4:"+str(common_peptides_experiments_1_4))
		print("Number of common peptides in experiments 2 and 4:"+str(common_peptides_experiments_2_4))
		print("Number of common peptides in experiments 3 and 4:"+str(common_peptides_experiments_3_4))
		print("Number of common peptides in experiments 1,2 and 3:"+str(common_peptides_experiments_1_2_3))
		print("Number of common peptides in experiments 1,2 and 4:"+str(common_peptides_experiments_1_2_4))
		print("Number of common peptides in experiments 1,3 and 4:"+str(common_peptides_experiments_1_3_4))
		print("Number of common peptides in experiments 2,3 and 4:"+str(common_peptides_experiments_2_3_4))
		print("Number of common peptides in experiments 1,2,3 and 4:"+str(common_peptides_experiments_1_2_3_4))
		print("Number of unique peptides in experiment 1:"+str(len(list(set(peptides1)))))
		print("Number of unique peptides in experiment 2:"+str(len(list(set(peptides2)))))
		print("Number of unique peptides in experiment 3:"+str(len(list(set(peptides3)))))
		print("Number of unique peptides in experiment 4:"+str(len(list(set(peptides4)))))
		print("Number of unique common peptides in experiments 1 and 2:"+str(len(list(set(peptides1_2)))))
		print("Number of unique common peptides in experiments 1 and 3:"+str(len(list(set(peptides1_3)))))
		print("Number of unique common peptides in experiments 2 and 3:"+str(len(list(set(peptides2_3)))))
		print("Number of unique common peptides in experiments 1 and 4:"+str(len(list(set(peptides1_4)))))
		print("Number of unique common peptides in experiments 2 and 4:"+str(len(list(set(peptides2_4)))))
		print("Number of unique common peptides in experiments 3 and 4:"+str(len(list(set(peptides3_4)))))
		print("Number of unique common peptides in experiments 1,2 and 3:"+str(len(list(set(peptides1_2_3)))))
		print("Number of unique common peptides in experiments 1,2 and 4:"+str(len(list(set(peptides1_2_4)))))
		print("Number of unique common peptides in experiments 1,3 and 4:"+str(len(list(set(peptides1_3_4)))))
		print("Number of unique common peptides in experiments 2,3 and 4:"+str(len(list(set(peptides2_3_4)))))
		print("Number of unique common peptides in experiments 1,2,3 and 4:"+str(len(list(set(peptides1_2_3_4)))))
		
		result_metrics_fid.write("Number of spectra in experiment 1:"+str(spectra_exp1))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 2:"+str(spectra_exp2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 3:"+str(spectra_exp3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of spectra in experiment 3:"+str(spectra_exp4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1 and 2:"+str(common_spectra_experiments_1_2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1 and 3:"+str(common_spectra_experiments_1_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 2 and 3:"+str(common_spectra_experiments_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1 and 4:"+str(common_spectra_experiments_1_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 2 and 4:"+str(common_spectra_experiments_2_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 3 and 4:"+str(common_spectra_experiments_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1,2 and 4:"+str(common_spectra_experiments_1_2_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1,2 and 3:"+str(common_spectra_experiments_1_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1,3 and 4:"+str(common_spectra_experiments_1_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 2,3 and 4:"+str(common_spectra_experiments_2_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common spectra in experiments 1,2,3 and 4:"+str(common_spectra_90_percent))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 1:"+str(peptides_exp1))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 2:"+str(peptides_exp2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 3:"+str(peptides_exp3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of peptides in experiment 3:"+str(peptides_exp4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1 and 2:"+str(common_peptides_experiments_1_2))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1 and 3:"+str(common_peptides_experiments_1_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 2 and 3:"+str(common_peptides_experiments_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1 and 4:"+str(common_peptides_experiments_1_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 2 and 4:"+str(common_peptides_experiments_2_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 3 and 4:"+str(common_peptides_experiments_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1,2 and 3:"+str(common_peptides_experiments_1_2_3))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1,2 and 4:"+str(common_peptides_experiments_1_2_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1,3 and 4:"+str(common_peptides_experiments_1_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 2,3 and 4:"+str(common_peptides_experiments_2_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of common peptides in experiments 1,2,3 and 4:"+str(common_peptides_experiments_1_2_3_4))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 1:"+str(len(list(set(peptides1)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 2:"+str(len(list(set(peptides2)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 3:"+str(len(list(set(peptides3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique peptides in experiment 4:"+str(len(list(set(peptides4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1 and 2:"+str(len(list(set(peptides1_2)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1 and 3:"+str(len(list(set(peptides1_3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 2 and 3:"+str(len(list(set(peptides2_3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1 and 4:"+str(len(list(set(peptides1_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 2 and 4:"+str(len(list(set(peptides2_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 3 and 4:"+str(len(list(set(peptides3_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1,2 and 3:"+str(len(list(set(peptides1_2_3)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1,2 and 4:"+str(len(list(set(peptides1_2_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1,3 and 4:"+str(len(list(set(peptides1_3_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 2,3 and 4:"+str(len(list(set(peptides2_3_4)))))
		result_metrics_fid.write("\n")
		result_metrics_fid.write("Number of unique common peptides in experiments 1,2,3 and 4:"+str(len(list(set(peptides1_2_3_4)))))
		
	result_metrics_fid.close()	
	final_report_fid.close();
	selected_spectra_fid.close();
	similarities_fid.close();
	quantitative_fid.close();
	mapping_fid.close();
	unified_spectra_search_fid.close();
	#Missing values imputation method
	print(mode)
	if mode=="1":
		
		queries_so_far_unified_list=queries_so_far;
		to_change=[-1]*queries_so_far_unified_list
		spectra_to_be_added=MSExperiment();
		num_of_spectra_to_be_added=0
		
		if os.path.exists(quantified_search_folder):
			number_of_parsed_files=0;
			search_filenames=list();
			for filename in os.listdir(quantified_search_folder):
				search_filenames.append(filename)
		if len(search_filenames)>0:		
			print("Quantified spectra search files were succesfully identified")
		
		if os.path.exists(filtered_files_folder):
			number_of_parsed_files=0;
			filtered_filenames=list();
			for filename in os.listdir(filtered_files_folder):
				filtered_filenames.append(filename)
		if len(filtered_filenames)>0:		
			print("Filtered files were succesfully identified")
		print("Missing values imputation procedure has started")
		quantified_filenames=list()
			
		for filename in search_filenames:
			print(filename+" is being parsed and examined!")
			quantified_spectra_fid=open(quantified_search_folder+"/"+filename,"r")
			for line1 in quantified_spectra_fid:
				try:
					index_file = line1.index("FILE=")
				except ValueError:
					index_file = -1
				if index_file!=-1:
					splitted_line=line1.split("FILE=")
					splitted_line_phase2=splitted_line[1].split("\\")
					print("It is the Mascot searched file for the quantified file file:"+str(splitted_line_phase2[-1].strip()))
					quantified_filenames.append(splitted_line_phase2[-1].strip())
					break
			quantified_spectra_search_info=list()
			quantified_spectra_search_fid=open(quantified_search_folder+"/"+filename,"r")
			
			number_of_lines=0
			queries_so_far=0
			queries_to_spectra_quant=list()
			order=0;
			#Associate searched files with filenames
			quant_to_filename=-1
			for f in filenames:
				names_parts=f.split(".")
				origin=names_parts[0][9:]
				if quantified_filenames[-1].find(origin)!=-1:
					quant_to_filename=order
				order=order+1
			print("Quant to filename="+str(quant_to_filename))
			flag=0;

			for line1 in quantified_spectra_search_fid:
				if "decoy_peptides" in line1:
					flag=1;
				string_to_search="q"+str(queries_so_far+1)+"_p1="
				if line1.startswith("scans="):
					data1=line1.split("=")
					queries_to_spectra_quant.append(int(data1[1].strip()))
				if flag==0 and line1.startswith(string_to_search):
					two_parts=line1.split("=")
					four_parts=two_parts[1].split(";")
					if two_parts[1].strip()=="-1":
						quantified_spectra_search_info.append([]);
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
					else:	
					
						words=four_parts[0].split(",")
						#print(words)	
						quantified_spectra_search_info.append([]);
						
						quantified_spectra_search_info[queries_so_far].append(words[4])
						quantified_spectra_search_info[queries_so_far].append(words[7].replace(".",","))
						four_parts[1].strip()
						different_proteins_full=four_parts[1].strip().split(",")
						different_proteins="";
						peptide_positions="";
						number_of_proteins_so_far=0;
						for protein in different_proteins_full:
							protein_info=protein.split(":")
							pept_pos=""
							for pos in xrange(1,len(protein_info)):
								if pos==len(protein_info)-1:
									pept_pos=pept_pos+protein_info[pos]
								else:
									pept_pos=pept_pos+protein_info[pos]+":"
							if number_of_proteins_so_far==0:
								different_proteins=different_proteins+protein_info[0].replace("\"","")
								peptide_positions=peptide_positions+pept_pos
							else:
								different_proteins=different_proteins+","+protein_info[0].replace("\"","")
								peptide_positions=peptide_positions+","+pept_pos
							number_of_proteins_so_far=number_of_proteins_so_far+1;
						quantified_spectra_search_info[queries_so_far].append(different_proteins)
						quantified_spectra_search_info[queries_so_far].append(peptide_positions)
					queries_so_far=queries_so_far+1;
			print("Queries so far"+str(queries_so_far))			
			file = pyopenms.FileHandler()
			quantified_file=MSExperiment();
			file.loadExperiment(quantified_folder+"/quantified_"+filenames[quant_to_filename], quantified_file)
			#store the mapping data for this filename in mapping_list and store similarities for this filename
			mapping_list=list();
			similarities_list=list()
			for i in range(30000):
				mapping_list.append(mapping_spectra[i][quant_to_filename])
				similarities_list.append(float(similarities[i][quant_to_filename]))
			print("quant_to_filename="+str(quant_to_filename))
			print("queries_so_far_unified_list="+str(queries_so_far_unified_list))			
			for i in range(queries_so_far_unified_list):
				try:		
					i_new=queries_to_spectra.index(i)
				except ValueError:
					i_new=-1;
				if i_new!=-1 and to_change[i]==-1:
					if unified_spectra_search_info[i_new][0]==-1:
						#check if this spectrum is mapped on the unified list from the examined quantified_file
						try:
							index = mapping_list.index(str(i))
						except ValueError:
							index = -1
						test=1;
						if index!=-1 and similarities_list[index]>float(similarity_threshold):
							test=0;
							#find the order of the spectrum in the quantified file						
							index_in_quantified_file=0;						
							for j in range(28000):
								if j==index:
									break;
								if mapping_list[j]!="-1":
									index_in_quantified_file=index_in_quantified_file+1;
							try:		
								index_in_queries=queries_to_spectra_quant.index(index_in_quantified_file)
							except ValueError:
								index_in_queries=-1;
							if index_in_queries!=-1:
								if quantified_spectra_search_info[index_in_queries][0]!=-1:	
									#save the spectrum from the quantified_file							
									to_change[i]=num_of_spectra_to_be_added
									num_of_spectra_to_be_added=num_of_spectra_to_be_added+1
									spectra_to_be_added.addSpectrum(quantified_file[index_in_quantified_file])
									unified_spectra_search_info[i_new][0]=quantified_spectra_search_info[index_in_queries][0]
		new_unified_spectra_file=MSExperiment()
		count1=0		
		count2=0
		for i in range(queries_so_far_unified_list):
			if to_change[i]==-1:
				count1=count1+1;
				new_unified_spectra_file.addSpectrum(unified_spectra[i])
			else:
				count2=count2+1;
				new_unified_spectra_file.addSpectrum(spectra_to_be_added[to_change[i]])
		print("Spectra kept:"+str(count1))
		print("Spectra changed:"+str(count2))				
		file.storeExperiment(output_folder+"unified_spectrum_list_missing_values_imputed.mzML",new_unified_spectra_file)				
		corrected_unified_spectra_list_fid = open(output_folder+"corrected_unified_spectrum_list_missing_values_imputed.mzML","w")
		unified_spectra_list_fid = open(output_folder+"unified_spectrum_list_missing_values_imputed.mzML","r")
		num_of_lines=0;	
		for line1 in unified_spectra_list_fid:
			num_of_lines=num_of_lines+1;
			corrected_unified_spectra_list_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
		corrected_unified_spectra_list_fid.close();	
		#Update quantified files		
		quantified_spectra_fid.close()
		print("Missing values imputation process in quantified files has started")
		quantified_filenames=list()
		for filename in search_filenames:
			spectra_to_be_added=MSExperiment();
			
			num_of_spectra_to_be_added=0
			print(filename+" is being parsed and examined!")
			quantified_spectra_fid=open(quantified_search_folder+"/"+filename,"r")
			for line1 in quantified_spectra_fid:
				try:
					index_file = line1.index("FILE=")
				except ValueError:
					index_file = -1
				if index_file!=-1:
					splitted_line=line1.split("FILE=")
					splitted_line_phase2=splitted_line[1].split("\\")
					print("It is the Mascot searched file for the quantified file file:"+str(splitted_line_phase2[-1].strip()))
					quantified_filenames.append(splitted_line_phase2[-1].strip())
					break
			quantified_spectra_search_info=list()
			quantified_spectra_search_fid=open(quantified_search_folder+"/"+filename,"r")
			number_of_lines=0
			queries_so_far=0
			queries_to_spectra_quant=list()
			order=0;
			quant_files_to_change=list()
			quant_files_to_change=[-1]*queries_so_far
			#Associate searched files with filenames
			quant_to_filename=-1
			for f in filenames:
				names_parts=f.split(".")
				origin=names_parts[0][9:]
				if quantified_filenames[-1].find(origin)!=-1:
					quant_to_filename=order
				order=order+1
			print("Quant to filename="+str(quant_to_filename))
			flag=0;

			for line1 in quantified_spectra_search_fid:
				if "decoy_peptides" in line1:
					flag=1;
				string_to_search="q"+str(queries_so_far+1)+"_p1="
				if line1.startswith("scans="):
					data1=line1.split("=")
					queries_to_spectra_quant.append(int(data1[1].strip()))
				if flag==0 and line1.startswith(string_to_search):
					two_parts=line1.split("=")
					four_parts=two_parts[1].split(";")
					if two_parts[1].strip()=="-1":
						quantified_spectra_search_info.append([]);
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
						quantified_spectra_search_info[queries_so_far].append(-1)
					else:	
					
						words=four_parts[0].split(",")
						#print(words)	
						quantified_spectra_search_info.append([]);
						
						quantified_spectra_search_info[queries_so_far].append(words[4])
						quantified_spectra_search_info[queries_so_far].append(words[7].replace(".",","))
						four_parts[1].strip()
						different_proteins_full=four_parts[1].strip().split(",")
						different_proteins="";
						peptide_positions="";
						number_of_proteins_so_far=0;
						for protein in different_proteins_full:
							protein_info=protein.split(":")
							pept_pos=""
							for pos in xrange(1,len(protein_info)):
								if pos==len(protein_info)-1:
									pept_pos=pept_pos+protein_info[pos]
								else:
									pept_pos=pept_pos+protein_info[pos]+":"
							if number_of_proteins_so_far==0:
								different_proteins=different_proteins+protein_info[0].replace("\"","")
								peptide_positions=peptide_positions+pept_pos
							else:
								different_proteins=different_proteins+","+protein_info[0].replace("\"","")
								peptide_positions=peptide_positions+","+pept_pos
							number_of_proteins_so_far=number_of_proteins_so_far+1;
						quantified_spectra_search_info[queries_so_far].append(different_proteins)
						quantified_spectra_search_info[queries_so_far].append(peptide_positions)
					queries_so_far=queries_so_far+1;
			quant_files_to_change=[-1]*queries_so_far			
			file = pyopenms.FileHandler()
			quantified_file=MSExperiment();
			file.loadExperiment(quantified_folder+"/quantified_"+filenames[quant_to_filename], quantified_file)
			#store the mapping data for this filename in mapping_list and store similarities for this filename
			mapping_list=list();
			similarities_list=list()
			for i in range(28000):
				mapping_list.append(int(mapping_spectra[i][quant_to_filename]))
				similarities_list.append(float(similarities[i][quant_to_filename]))			
			for i in range(queries_so_far):
				try:		
					index_in_searched_file=queries_to_spectra_quant.index(i)
				except ValueError:
					index_in_searched_file=-1;
				if index_in_searched_file==-1 or quantified_spectra_search_info[index_in_searched_file][0]==-1:
					index_in_filtered_file=0;						
					for j in range(28000):
						if mapping_list[j]!=-1:
							if index_in_filtered_file==i:
								break
							index_in_filtered_file=index_in_filtered_file+1
					#find position in the unified list
					try:
						index_in_unified_list = mapping_list[index_in_filtered_file]
					except ValueError:
						index_in_unified_list = -1					
					if index_in_unified_list !=-1 and similarities_list[index_in_filtered_file]>float(similarity_threshold):
						try:
							index_in_unified_search_file=queries_to_spectra.index(index_in_unified_list)
						except ValueError:
							index_in_unified_search_file=-1;
						if index_in_unified_search_file!=-1:
							if unified_spectra_search_info[index_in_unified_search_file][0]!=-1:
								#get peak list from unified_spectra_list
								unified_peaks=new_unified_spectra_file[index_in_unified_list].get_peaks();
								#get peak list from quantified_files
								quantified_peaks=quantified_file[i].get_peaks();
								
								if len(quantified_peaks)>0 and len(unified_peaks)>0:
									mz_values_old_unified, intensities_old_unified = zip(*unified_peaks)
									mz_values_unified = [float(v) for v in mz_values_old_unified]
									intensities_unified = [float(v) for v in intensities_old_unified]
									mz_values_old_quant, intensities_old_quant = zip(*quantified_peaks)
									mz_values_quant = [float(v) for v in mz_values_old_quant]
									intensities_quant = [float(v) for v in intensities_old_quant]
									other_mz_values_quant=[mz for mz in mz_values_quant if (mz<(125.5) or mz>(131.5))]
									p=0
									while (p<len(other_mz_values_quant)):
										position=mz_values_quant.index(other_mz_values_quant[p])
										del mz_values_quant[position]
										del intensities_quant[position]
										p=p+1
									reporter_mz_values_unified=[mz for mz in mz_values_unified if 125.5<mz<131.5]
									start_position_to_insert=mz_values_unified.index(reporter_mz_values_unified[0])
									p=0
									while (p<len(reporter_mz_values_unified)):
										position=mz_values_unified.index(reporter_mz_values_unified[p])
										del mz_values_unified[position]
										del intensities_unified[position]
										p=p+1
									for p in range(len(mz_values_quant)):
										mz_values_unified.insert(start_position_to_insert+p,mz_values_quant[p])
										intensities_unified.insert(start_position_to_insert+p,intensities_quant[p])
									new_peaks=zip(mz_values_unified,intensities_unified)
									new_peaks2=[list(a) for a in zip(mz_values_unified, intensities_unified)]
									new_unified_spectra_file[index_in_unified_list].set_peaks(numpy.asarray(new_peaks2,dtype=numpy.float32))
									quant_files_to_change[i]=num_of_spectra_to_be_added
									num_of_spectra_to_be_added=num_of_spectra_to_be_added+1
									spectra_to_be_added.addSpectrum(new_unified_spectra_file[index_in_unified_list])
			new_quantified_file=MSExperiment()
			count1=0;
			count2=0;
			for i in range(queries_so_far):
				if quant_files_to_change[i]==-1:
					count1=count1+1;
					new_quantified_file.addSpectrum(quantified_file[i])
				else:
					count2=count2+1;
					new_quantified_file.addSpectrum(spectra_to_be_added[quant_files_to_change[i]])
			print("Spectra kept:"+str(count1))
			print("Spectra changed:"+str(count2))	
			if not os.path.exists(output_folder+"new_quantified_files/"):
				os.makedirs(output_folder+"new_quantified_files/")			
			file.storeExperiment(output_folder+"new_quantified_files/quantified_"+filenames[quant_to_filename],new_quantified_file)				
			corrected_quantified_file_fid = open(output_folder+"new_quantified_files/corrected_quantified_"+filenames[quant_to_filename],"w")
			quantified_file_fid = open(output_folder+"new_quantified_files"+"/quantified_"+filenames[quant_to_filename],"r")
			num_of_lines=0;	
			for line1 in quantified_file_fid:
				num_of_lines=num_of_lines+1;
				corrected_quantified_file_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
			corrected_quantified_file_fid.close();
			quantified_file_fid.close();
	

if __name__ == "__main__":

	
	unified_spectrum_list_filename=sys.argv[1]
	mapping_spectra_filename=sys.argv[2]
	quantitative_values_filename=sys.argv[3]
	classes_filename=sys.argv[4]
	similarities_mapping_spectra_filename=sys.argv[5]
	filtered_files_folder=sys.argv[6]
	quantified_search_folder=sys.argv[7]
	unified_spectra_search_filename=sys.argv[8]
	similarity_threshold=float(sys.argv[9])
	selected_spectra_filename=sys.argv[10]
	scaffold_peptide_report_filename=sys.argv[11]
	quantified_folder=sys.argv[12]
	mode=sys.argv[13]
	replicates=int(sys.argv[14])
	missing_values_filename=sys.argv[15]
	use_inputs_vector_filename=sys.argv[16]
	tstamp = time.strftime('%Y_%m_%d')
	output_folder=str(tstamp)+'_'+str(time.time())+'/'
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	generate_final_report(unified_spectrum_list_filename,mapping_spectra_filename,quantitative_values_filename,classes_filename, similarities_mapping_spectra_filename,filtered_files_folder,quantified_search_folder,unified_spectra_search_filename,similarity_threshold,selected_spectra_filename,scaffold_peptide_report_filename,quantified_folder, mode,output_folder,replicates,missing_values_filename,use_inputs_vector_filename)
	
