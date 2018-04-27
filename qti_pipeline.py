####################################
#This code is the outcome of the project Innovative Processing of TMT Proteomics. This is a common project of InSyBio LTD and Nestle Institute of Health Sciences. 
#This program implements the step 2 code
####################################

#example on calling the code:
#sudo python qti_pipeline.py input_data_example reference_proteins_peptides_peaks_new.txt 10 2 0.8 0.45 0.45 0.05 0.1 initial_good_solutions.txt goal_significances_step1.txt Swissprot_Human_BLA_2014_dec_08.fasta missing_values_example.txt 0.02 900 0.5 50 0.2 0.02 classes_example.txt 10 2 0.8 0.45 0.45 0.05 0.1 goal_significances_step1.txt 5 3 use_inputs_vector_example.txt 1 0.7 0 > output.txt & 


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
import multiprocessing as mp
import glob
from collections import defaultdict
from bisect import bisect_left
from subprocess import Popen
from svmutil import *
import warnings
import logging
from heapq import nlargest
import re
import csv
import pdb
from qti_utils import insybio_automatic_pipeline_construction, apply_best_solution, quantification_and_spectra_alignment, step3_supervised_analysis, generate_final_report
#from qti-utils_test import count_intervals
if __name__ == '__main__':
	
	min_values_step1 = [0.0, 0.0, 0.1, 1.0,  5.0,  3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.01, 0.0, 0.01,  50,    0.0, 0.0, 0.0, 1.0,   3.0,  1.0,  0.0, 0.0, 0.0]
	max_values_step1 = [1.0, 3.0, 2.0, 20.0, 20.0, 8.0, 1.0, 2.0, 5.0, 1.0, 5.0, 0.99, 3.0, 100.0, 300.0, 1.0, 3.0, 3.0, 400.0, 50.0, 20.0, 3.0, 1.0, 2.0]
	randomness=int(sys.argv[34])
	if randomness==0:
		random.seed(0)
	else:
		random.seed()
	
	
	input_data_folder					  = sys.argv[1]
	reference_peptides_filename			= sys.argv[2]
	population_step1					   = int(sys.argv[3])
	generations_step1					  = int(sys.argv[4])
	semi_random_initialization_probability_step1   = float(sys.argv[5])
	two_points_crossover_probability_step1		 = float(sys.argv[6])
	arithmetic_crossover_probability_step1		 = float(sys.argv[7])
	mutation_probability_step1				  = float(sys.argv[8])
	gaussian_mutation_variance_proportion_step1  = float(sys.argv[9])
	initial_good_solutions_step1				 = sys.argv[10]
	goal_significances_filename_step1	  = sys.argv[11]
	database_name						  = sys.argv[12]
	missing_values_filename				= sys.argv[13]
	goal_significances=[]
	with open(goal_significances_filename_step1) as tsv:
		for line in csv.reader(tsv, delimiter='\t'):
			for val in line:
				goal_significances.append(float(val))
	tstamp = time.strftime('%Y_%m_%d')
	output_folder=str(tstamp)+'_'+str(time.time())+'/'
	output_folder_step1=output_folder+'step1/'
	if not os.path.exists(output_folder_step1):
		os.makedirs(output_folder_step1)
	[best_solution_path]=insybio_automatic_pipeline_construction(
			input_data_folder,
			reference_peptides_filename,
			population_step1,
			generations_step1,
			semi_random_initialization_probability_step1,
			two_points_crossover_probability_step1,
			arithmetic_crossover_probability_step1,
			mutation_probability_step1,
			gaussian_mutation_variance_proportion_step1,
			min_values_step1,
			max_values_step1,
			goal_significances,
			initial_good_solutions_step1,
			database_name,
			missing_values_filename,
			output_folder_step1)
	print('Step 1 Successfully Finished!')		
	#step2
	output_folder_step2=output_folder+'step2/'
	#Initialize random generator seed with the current local time in miliseconds
	if not os.path.exists(output_folder_step2):
		os.makedirs(output_folder_step2)
	
	
	minimum_similarity=float(sys.argv[14])
	retention_time_tolerance=float(sys.argv[15])
	quantified_spectra_thres=float(sys.argv[16])
	max_allowed_variance_in_peak_number=int(sys.argv[17])
	abs_sim_thres=float(sys.argv[18])
	min_first_to_second_score_distance=float(sys.argv[19])
	#missing_values_filename = sys.argv[9]
	apply_best_solution(input_data_folder,best_solution_path,quantified_spectra_thres, ubuntu_flag, missing_values_filename, output_folder_step2)
	start_time = time.time()
	[unified_list_filename, mapping_filename, quantitative_filename,similarities_filename]=quantification_and_spectra_alignment(input_data_folder,max_allowed_variance_in_peak_number,minimum_similarity,retention_time_tolerance,quantified_spectra_thres,abs_sim_thres,min_first_to_second_score_distance, missing_values_filename, output_folder_step2)
	print("--- %s seconds ---" % (time.time() - start_time))
	
	print('Step 2 Successfully Finished!')
	
	
	#step3
	min_values=[0.0,0.001,0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	max_values=[2.0,1000.0,1000.0,0.3,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]
	
	output_folder_step3=output_folder+'step3/'
	if not os.path.exists(output_folder_step3):
		os.makedirs(output_folder_step3)
	labels_filename=sys.argv[20]
	population_step3=int(sys.argv[21])
	generations_step3=int(sys.argv[22])
	semi_random_initialization_probability_step3=float(sys.argv[23])
	two_points_crossover_probability_step3=float(sys.argv[24])
	arithmetic_crossover_probability_step3=float(sys.argv[25])
	mutation_probability_step3=float(sys.argv[26])
	gaussian_mutation_variance_proportion_step3=float(sys.argv[27])
	goal_significances_filename_step3=sys.argv[28]
	num_of_folds=int(sys.argv[29])
	step2_output_folder=output_folder
	replicates=int(sys.argv[30])
	#missing_values_filename=sys.argv[30]
	use_inputs_vector_filename=sys.argv[31]
	supervised_mode=int(sys.argv[32])
	goal_significances_fid = open(goal_significances_filename_step3,"r")
	tstamp = time.strftime('%Y_%m_%d')
	
	number_of_lines=0
	for line in goal_significances_fid:
		number_of_lines=number_of_lines+1
	goal_significances_fid.close()
	goal_significances=list()
	goal_significances_fid = open(goal_significances_filename,"r")
	quantitative_fid=open(quantitative_filename,"r")
	number_of_lines=0
	for line in quantitative_fid:
		number_of_lines=number_of_lines+1
		min_values.append(0.0)
		max_values.append(1.0)
	number_of_lines=0
	for line in goal_significances_fid:
		words=line.split("\t")
		for i in range (len(words)):
			goal_significances.append(float(words[i].strip()))
		number_of_lines=number_of_lines+1
	
	selected_spectra_filename=step3_supervised_analysis(unified_list_filename,mapping_filename,quantitative_filename,labels_filename,population_step3,generations_step3,semi_random_initialization_probability_step3,two_points_crossover_probability_step3,arithmetic_crossover_probability_step3,mutation_probability_step3,gaussian_mutation_variance_proportion_step3,min_values,max_values,goal_significances, database_name,num_of_folds,output_folder,output_folder_step3,replicates, missing_values_filename,use_inputs_vector_filename,supervised_mode)
	print('Step 3 Successfully Finished!')
	#step4
	
	
	output_folder_step4=output_folder+'step4/'
	if not os.path.exists(output_folder_step4):
		os.makedirs(output_folder_step4)
	unified_spectrum_list_filename=output_folder_step3+'/corrected+'+unified_list_filename
	mapping_spectra_filename=mapping_filename
	quantitative_values_filename=quantitative_filename
	classes_filename=labels_filename
	similarities_mapping_spectra_filename=similarities_filename
	filtered_files_folder=output_folder_step2+'filtered_files/'
	quantified_search_folder=raw_input("Search the quantified mzml files file with Mascot, upload the results on a new folder on the server and write here the absolute path of the folder where they have been upploaded.")
	
	unified_spectra_search_filename = raw_input("Search the unified spectral list file with Mascot, upload it in the server and write here the absolute path of the .dat mascot searched unified spectral list file and press enter.")
	
	similarity_threshold=float(sys.argv[33])
	scaffold_peptide_report_filename = raw_input("Analyze the searched unified spectral list file with Scaffold export the quantified peptide report, upload it in the server and write here the absolute path of the .txt  quantified peptide report exported from Scaffold using the unified spectral list file and press enter.")
	
	
	quantified_folder=output_folder_step2+'quantified_files/'
	mode=1
	#replicates=int(sys.argv[14])
	#missing_values_filename=sys.argv[15]
	#use_inputs_vector_filename=sys.argv[16]
	#tstamp = time.strftime('%Y_%m_%d')
	#output_folder=str(tstamp)+'_'+str(time.time())+'/'
	#if not os.path.exists(output_folder):
	#	os.makedirs(output_folder)
	[new_unified_spectrum_list,new_quantified_files_folder]=generate_final_report(unified_spectrum_list_filename,mapping_spectra_filename,quantitative_values_filename,classes_filename, similarities_mapping_spectra_filename,filtered_files_folder,quantified_search_folder,unified_spectra_search_filename,similarity_threshold,selected_spectra_filename,scaffold_peptide_report_filename,quantified_folder, mode,output_folder_step4,replicates,missing_values_filename,use_inputs_vector_filename)
	unified_spectra_search_filename = raw_input("Search the corrected unified spectral list file with Mascot, upload it in the server and write here the absolute path of the .dat mascot searched corrected unified spectral list file and press enter.")
	quantified_search_folder=raw_input("Search the corrected quantified mzml files file (they can be found in Step4 folder) with Mascot, upload the results on a new folder on the server and write here the absolute path of the folder where they have been upploaded.")
	
	
	scaffold_peptide_report_filename = raw_input("Analyze the searched corrected unified spectral list file with Scaffold export the quantified peptide report, upload it in the server and write here  the absolute path of the .txt  quantified peptide report exported from Scaffold using the corrected unified spectral list file and press enter.")
	
	
	mode=0
	
	[new_unified_spectrum_list,new_quantified_files_folder]=generate_final_report(new_unified_spectrum_list,mapping_spectra_filename,quantitative_values_filename,classes_filename, similarities_mapping_spectra_filename,filtered_files_folder,quantified_search_folder,unified_spectra_search_filename,similarity_threshold,selected_spectra_filename,scaffold_peptide_report_filename,quantified_folder, mode,output_folder_step4,replicates,missing_values_filename,use_inputs_vector_filename)
	print('Step 4 Successfully Finished!')