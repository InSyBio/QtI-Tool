#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################
#This code is the outcome of the project Innovative Processing of TMT Proteomics. This is a common project of InSyBio LTD and Nestle Institute of Health Sciences. 
#This program implements the step 2 code
####################################

import multiprocessing as mp
from pyopenms import *
import os, shutil
import random
import math
import copy
import sys
import time
import subprocess
from heapq import nlargest
import numpy as np
from collections import defaultdict
from bisect import bisect_left
from subprocess import Popen

def apply_best_solution(input_folder,best_solution_file,quantified_spectra_thres, missing_values_filename, output_folder):
	best_solution_fid = open(best_solution_file,"r")
	number_of_lines=0
	best_solution=list()
	missing_values_fid      = open(missing_values_filename,'r')
	#parse the missing values filename
    missing_values_list=list()
    missing_values_raws=0;
    print('missing values filename is being parsed')
    for line1 in missing_values_fid:
        words=line1.split('\t')
        missing_values_list.append([])
        for i in range(len(words)):
            missing_values_list[missing_values_raws].append(words[i].strip())
        missing_values_raws=missing_values_raws+1;

	for line in best_solution_fid:
		number_of_lines=number_of_lines+1
		words=line.split("\t")
		print(words)
		print("Elements found="+str(len(words)))
		for i in range(len(words)):
			
			best_solution.append(float(words[i]))
		break
	print(best_solution)
	print_best_solution_in_detail(best_solution)
	write_best_solution_in_detail(best_solution,output_folder)
	filenames=list()
	#load mzML files for evaluation
	print(input_folder)
	if os.path.exists(input_folder):
		number_of_parsed_files=0
		filenames2=list()
		fully_quantified_average=0
		quantified_average=0
		
		for filename in os.listdir(input_folder):
			missing_values_raw=list()
			for i in range(len(missing_values_list)):
				if filename==missing_values_list[i][0]:
					for j in xrange(1, len(missing_values_list[0])):
						missing_values_raw.append(missing_values_list[i][j])
			print(number_of_parsed_files)
			exp = MSExperiment()
			file = pyopenms.FileHandler()
			refiltered_file2=MSExperiment()
			file.loadExperiment(input_folder+filename, exp)
			parsed_files2=exp
			filtered_file=filter_experiments(best_solution,parsed_files2)
			if best_solution[16]<1:
				peaks_picked=MSExperiment()
				refiltered_file=peaks_picked
				pp=PeakPickerHiRes()
				param=pyopenms.Param()
				pp.setParameters(param)
				pp.pickExperiment(filtered_file,refiltered_file)
				filename_out = output_folder+"filtered_files/filtered_"+str(filename)
				file.storeExperiment(filename_out, refiltered_file)
			else:			
				filename_out = output_folder+"filtered_files/filtered_"+str(filename)
				file.storeExperiment(filename_out, filtered_file)
			filenames.append(filename)
			
			if best_solution[16]<1:
				peaks_picked=MSExperiment()
				refiltered_file=peaks_picked
				pp=PeakPickerHiRes()
				param=pyopenms.Param()
				pp.setParameters(param)
				pp.pickExperiment(filtered_file,refiltered_file)
			
				reporter_ions_variation=0
				fully_quantified=0
				identified_xtandem=0
				identified_xtandem2=0
				proteins_scaffold=0
				quantified=0
				exact_reporter_ions=0
				ms_ms2=0
				to_delete=list()
				print("File size="+str(refiltered_file.size()))
				for k in range(refiltered_file.size()):
					if refiltered_file[k].getMSLevel() == 1:	
						to_delete.append(k)	
						continue
					peaks = refiltered_file[k].get_peaks()
					#print("peaks")
					#print(peaks)
					if len(peaks)==0:
						mz_values=list()
						intensities=list()
					else:
						mz_values_old, intensities_old = zip(*peaks)
						#mz_values_old=list(peaks[0])
						mz_values = [float(v) for v in mz_values_old]
						#intensities_old=list(peaks[1])
						intensities = [float(v) for v in intensities_old]
					exact_reporter_ions=0
					dictionary_result=count_intervals(mz_values, [125.5, 126.5,127.5, 128.5, 129.5, 130.5, 131.5])
					number_of_non_missing_channels=0
					for chan in range(len(missing_values_raw)):
						if(missing_values_raw[chan]=='0'):
							number_of_non_missing_channels+=1;
					if dictionary_result[126.5]>=1 and missing_values_raw[0]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if dictionary_result[127.5]>=1 and missing_values_raw[1]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[128.5]>=1 and missing_values_raw[2]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[129.5]>=1 and missing_values_raw[3]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[130.5]>=1 and missing_values_raw[4]=='0':
						exact_reporter_ions=exact_reporter_ions+1				
					if dictionary_result[131.5]>=1 and missing_values_raw[5]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if exact_reporter_ions>=quantified_spectra_thres*number_of_non_missing_channels:
						refiltered_file2.addSpectrum(refiltered_file[k])
						quantified=quantified+1
					else:
						to_delete.append(k)
							
					if exact_reporter_ions==number_of_non_missing_channels:
						fully_quantified=fully_quantified+1
					
					ms_ms2=ms_ms2+1
				print("ms_ms2="+str(ms_ms2))
				
			else:
				to_delete=list()
				refiltered_file=filtered_file
				reporter_ions_variation=0
				fully_quantified=0
				identified_xtandem=0
				identified_xtandem2=0
				proteins_scaffold=0
				quantified=0
				exact_reporter_ions=0
				ms_ms2=0
				print("File size="+str(refiltered_file.size()))
				for k in range(refiltered_file.size()):
					if refiltered_file[k].getMSLevel() == 1:
						to_delete.append(k)	
						continue
					peaks = refiltered_file[k].get_peaks()
					if len(peaks)==0:
						mz_values=list()
						intensities=list()
					else:					
						mz_values_old, intensities_old = zip(*peaks)
						#mz_values_old=list(peaks[0])
						mz_values = [float(v) for v in mz_values_old]
						#intensities_old=list(peaks[1])
						intensities = [float(v) for v in intensities_old]
						
					exact_reporter_ions=0
					dictionary_result=count_intervals(mz_values, [125.5, 126.5,127.5, 128.5, 129.5, 130.5, 131.5])
					number_of_non_missing_channels=0
					for chan in range(len(missing_values_raw)):
						if(missing_values_raw[chan]=='0'):
							number_of_non_missing_channels+=1;
					if dictionary_result[126.5]>=1 and missing_values_raw[0]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if dictionary_result[127.5]>=1 and missing_values_raw[1]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[128.5]>=1 and missing_values_raw[2]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[129.5]>=1 and missing_values_raw[3]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[130.5]>=1 and missing_values_raw[4]=='0':
						exact_reporter_ions=exact_reporter_ions+1				
					if dictionary_result[131.5]>=1 and missing_values_raw[5]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if exact_reporter_ions>=quantified_spectra_thres*number_of_non_missing_channels:
						refiltered_file2.addSpectrum(refiltered_file[k])
						quantified=quantified+1
					else:
						to_delete.append(k)
							
					if exact_reporter_ions==number_of_non_missing_channels:
						fully_quantified=fully_quantified+1
			
					ms_ms2=ms_ms2+1
				print("ms_ms2="+str(ms_ms2))
				
			print(quantified)
			print(fully_quantified)	
			fully_quantified_average=fully_quantified_average+fully_quantified
			quantified_average=quantified_average+quantified
			filename_out = output_folder+"quantified_spectra_files/quantified_filtered_"+str(filename)
			file.storeExperiment(filename_out, refiltered_file2)
			del refiltered_file2
			number_of_parsed_files=number_of_parsed_files+1
		print(str(number_of_parsed_files)+" mzML files were successfully filtered and stored")
		fully_quantified_average=fully_quantified_average/float(number_of_parsed_files)
		quantified_average=quantified_average/float(number_of_parsed_files)
		del refiltered_file2
		del refiltered_file
		del parsed_files2
		del exp
		del file
		print("Average number of quantified MSMS spectra per experiment is:"+str(quantified_average))
		print("Average number of fully quantified MSMS spectra per experiment is:"+str(fully_quantified_average))

def count_intervals(sequence, intervals):
    count = defaultdict(int)
    intervals.sort()
    for item in sequence:
        pos = bisect_left(intervals, item)
        if pos == len(intervals):
            count[None] += 1
        else:
            count[intervals[pos]] += 1
    return count

def write_best_solution_in_detail(individual,output_folder):
	#Individual (list): Representing an individual solution
	detailed_best_solution_fid=open(output_folder+"detailed_best_solution"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	if individual[0]>individual[6] and individual[6]>individual[9] and individual[9]>individual[15]:
		detailed_best_solution_fid.write("Denoising->Precursor Removal ->Normalization ->PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")		
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")	
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("PPM tolerance="+str(individual[3]))
			detailed_best_solution_fid.write("\n")	
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")	
	
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")	
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")	

		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")	
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")	
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")	
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")	
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")	
			
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")	
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")	
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")	

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")	
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")	 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")	
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")	
			
	elif individual[0]>=individual[6] and individual[6]>=individual[15] and individual[15]>=individual[9]:
		detailed_best_solution_fid.write("#Denoising->Precursor Removal ->PeakPicking ->Normalization ->PickFiltering")
				
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")	
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")	
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")	
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")	
	
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")	
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[0]>=individual[9] and individual[9]>=individual[6] and individual[6]>=individual[15]:
		detailed_best_solution_fid.write("#Denoising->Normalization ->Precursor Removal->PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width=",+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")

	elif individual[0]>=individual[9] and individual[9]>=individual[15] and individual[15]>=individual[6]:
		detailed_best_solution_fid.write("#Denoising->Normalization ->PeakPicking ->Precursor Removal->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")		
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")		
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")

	elif individual[0]>=individual[15] and individual[15]>=individual[9] and individual[9]>=individual[6]:
		detailed_best_solution_fid.write("#Denoising->PeakPicking ->Normalization ->Precursor Removal->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")	
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peak Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[0]>=individual[15] and individual[15]>=individual[6] and individual[6]>=individual[9]:
		detailed_best_solution_fid.write("#Denoising->PeakPicking ->Precursor Removal-> Normalization ->PickFiltering")
		detailed_best_solution_fid.write("\n")
		
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")		
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[6]>=individual[0] and individual[0]>=individual[9] and individual[9]>=individual[15]:
		detailed_best_solution_fid.write("#Precursor Removal-> Denoising->Normalization ->PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")			
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")

		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))	
		detailed_best_solution_fid.write("\n")

		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[6]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[9]:
		detailed_best_solution_fid.write("#Precursor Removal-> Denoising->PeakPicking ->Normalization ->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")

		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
		detailed_best_solution_fid.write("\n")

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold=",str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[9]>=individual[0] and individual[0]>=individual[6] and individual[6]>=individual[15]:
		detailed_best_solution_fid.write("#Normalization-> Denoising->Precursor Removal -> PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")	

		#Normalization:		
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[9]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[6]:
		detailed_best_solution_fid.write("#Normalization-> Denoising-> PeakPicking ->Precursor Removal -> PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	elif individual[6]>=individual[9] and individual[9]>=individual[0]:
		detailed_best_solution_fid.write("#Precursor Removal-> Normalization ->Denoising ->PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")
	
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")

		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")

	elif individual[9]>=individual[6] and individual[6]>=individual[0]:
		detailed_best_solution_fid.write("#Normalization-> Precursor Removal->Denoising ->PeakPicking ->PickFiltering")
		detailed_best_solution_fid.write("\n")
			
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
		detailed_best_solution_fid.write("\n")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")

		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")		

		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied") 
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
	
	else:
		detailed_best_solution_fid.write("#default")
		detailed_best_solution_fid.write("\n")
		
		#Denoising		
		if individual[1]<1:
			detailed_best_solution_fid.write("No denoising is applied")
			detailed_best_solution_fid.write("\n")
			nothing=0
		elif individual[1]<2:
			detailed_best_solution_fid.write("Gaussian denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Gaussian width="+str(individual[2]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("PPM tolerance"+str(individual[3]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Savitzky Golay denoising")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("frame_length="+str(int(round(individual[4]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("polynomial_order="+str(int(round(individual[5]))))
			detailed_best_solution_fid.write("\n")
		print("Denoising finished")
		#Precursor Removal
		if individual[7]<1:
			detailed_best_solution_fid.write("no precursor removal")
			detailed_best_solution_fid.write("\n")
			nothing=0
		else:
			detailed_best_solution_fid.write("ParentPeakMower Method")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("ParentPeakMower_window_size="+str(float(individual[8])))
			detailed_best_solution_fid.write("\n")
		print("Precursor removal finished")
		#Normalization:
		if individual[10]<1:
			detailed_best_solution_fid.write("No normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<2:
			detailed_best_solution_fid.write("Direct Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<3:
			detailed_best_solution_fid.write("Scaler Normalization")
			detailed_best_solution_fid.write("\n")
		elif individual[10]<4:
			detailed_best_solution_fid.write("Logarithmic Normalization")
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("Bern Normalization")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("Bern Normalization Threshold="+str(individual[11]))
			detailed_best_solution_fid.write("\n")			
		print("normalization finished")
		#PeakPicking
		if individual[16]<1:
			detailed_best_solution_fid.write("Default PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			
		elif individual[16]<2:
			detailed_best_solution_fid.write("PeakPickerHiRes is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[17]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:win_len="+str(individual[18]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("PeakPickerCWT is applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("signal_to_noise="+str(individual[21]))
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("peak_width="+str(individual[22]))
			detailed_best_solution_fid.write("\n")
		print("Peak picking finished")
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			detailed_best_solution_fid.write("No Peak Filtering")
			detailed_best_solution_fid.write("\n")
		elif individual[12]<2:
			detailed_best_solution_fid.write("Threshold Mower Peakk Filtering Method is Applied")
			detailed_best_solution_fid.write("\n") 
			detailed_best_solution_fid.write("Threshold="+str(individual[13]))
			detailed_best_solution_fid.write("\n")
		else:
			detailed_best_solution_fid.write("NLargest Peak Filtering Method is Applied")
			detailed_best_solution_fid.write("\n")
			detailed_best_solution_fid.write("N="+str(int(round(individual[14]))))
			detailed_best_solution_fid.write("\n")
		print("Peak filtering finished")
	detailed_best_solution_fid.close()	
		
def print_best_solution_in_detail(individual):
	#Individual (list): Representing an individual solution
	#parsed_data (list): a list of MSExperiment storing parsed mzML files
	
	#print(len(parsed_data))
	if individual[0]>individual[6] and individual[6]>individual[9] and individual[9]>individual[15]:
		print("Denoising->Precursor Removal ->Normalization ->PeakPicking ->PickFiltering")
				
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
	
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))

		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
			
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width=",+str(individual[22]))

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
			
	elif individual[0]>=individual[6] and individual[6]>=individual[15] and individual[15]>=individual[9]:
		print("#Denoising->Precursor Removal ->PeakPicking ->Normalization ->PickFiltering")
				
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
	
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[0]>=individual[9] and individual[9]>=individual[6] and individual[6]>=individual[15]:
		print("#Denoising->Normalization ->Precursor Removal->PeakPicking ->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))

	elif individual[0]>=individual[9] and individual[9]>=individual[15] and individual[15]>=individual[6]:
		print("#Denoising->Normalization ->PeakPicking ->Precursor Removal->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))		
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))		
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))

	elif individual[0]>=individual[15] and individual[15]>=individual[9] and individual[9]>=individual[6]:
		print("#Denoising->PeakPicking ->Normalization ->Precursor Removal->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))	
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[0]>=individual[15] and individual[15]>=individual[6] and individual[6]>=individual[9]:
		print("#Denoising->PeakPicking ->Precursor Removal-> Normalization ->PickFiltering")
		
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))		
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width=",+str(individual[22]))
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[6]>=individual[0] and individual[0]>=individual[9] and individual[9]>=individual[15]:
		print("#Precursor Removal-> Denoising->Normalization ->PeakPicking ->PickFiltering")
			
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))			
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))

		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))	

		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[6]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[9]:
		print("#Precursor Removal-> Denoising->PeakPicking ->Normalization ->PickFiltering")
			
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))

		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[9]>=individual[0] and individual[0]>=individual[6] and individual[6]>=individual[15]:
		print("#Normalization-> Denoising->Precursor Removal -> PeakPicking ->PickFiltering")	

		#Normalization:		
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[9]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[6]:
		print("#Normalization-> Denoising-> PeakPicking ->Precursor Removal -> PickFiltering")
			
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	elif individual[6]>=individual[9] and individual[9]>=individual[0]:
		print("#Precursor Removal-> Normalization ->Denoising ->PeakPicking ->PickFiltering")
	
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))

		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
		
		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))

	elif individual[9]>=individual[6] and individual[6]>=individual[0]:
		print("#Normalization-> Precursor Removal->Denoising ->PeakPicking ->PickFiltering")
			
		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))

		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))		

		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
	
	else:
		print("#default")
		
		#Denoising		
		if individual[1]<1:
			print("No denoising is applied")
			nothing=0
		elif individual[1]<2:
			print("Gaussian denoising")
			print("Gaussian width="+str(individual[2]))
			print("PPM tolerance"+str(individual[3]))
		else:
			print("Savitzky Golay denoising")
			print("frame_length="+str(int(round(individual[4]))))
			print("frame_length="+str(int(round(individual[5]))))
	
		#Precursor Removal
		if individual[7]<1:
			print("no precursor removal")
			nothing=0
		else:
			print("ParentPeakMower Method")
			print("ParentPeakMower_window_size="+str(float(individual[8])))

		#Normalization:
		if individual[10]<1:
			print("No normalization")
		elif individual[10]<2:
			print("Direct Normalization")
		elif individual[10]<3:
			print("Scaler Normalization")
		elif individual[10]<4:
			print("Logarithmic Normalization")
		else:
			print("Bern Normalization")
			print("Bern Normalization Threshold="+str(individual[11]))			

		#PeakPicking
		if individual[16]<1:
			print("Default PeakPickerHiRes is applied")
			
		elif individual[16]<2:
			print("PeakPickerHiRes is applied")
			print("signal_to_noise="+str(individual[17]))
			print("SignalToNoise:win_len="+str(individual[18]))
			print("SignalToNoise:bin_count="+str(int(round(individual[19]))))
			print("SignalToNoise:min_required_elements="+str(int(round(individual[20]))))
		else:
			print("PeakPickerCWT is applied")
			print("signal_to_noise="+str(individual[21]))
			print("peak_width="+str(individual[22]))

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			print("No Peak Filtering")
		elif individual[12]<2:
			print("Threshold Mower Peakk Filtering Method is Applied") 
			print("Threshold="+str(individual[13]))
		else:
			print("NLargest Peak Filtering Method is Applied")
			print("N="+str(int(round(individual[14]))))
		
def filter_experiments(individual,parsed_data):
	#Individual (list): Representing an individual solution
	#parsed_data (list): a list of MSExperiment storing parsed mzML files
	exp_picked=MSExperiment()
	parsed_data_picked=exp_picked
	#print(len(parsed_data))
	if individual[0]>individual[6] and individual[6]>individual[9] and individual[9]>individual[15]:
		print("#Denoising->Precursor Removal ->Normalization ->PeakPicking ->PickFiltering")
		
		
				
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1				
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum

		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[0]>=individual[6] and individual[6]>=individual[15] and individual[15]>=individual[9]:
		print("#Denoising->Precursor Removal ->PeakPicking ->Normalization ->PickFiltering")
				
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)


		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists			
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[0]>=individual[9] and individual[9]>=individual[6] and individual[6]>=individual[15]:
		print("#Denoising->Normalization ->Precursor Removal->PeakPicking ->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)			
		

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked

	elif individual[0]>=individual[9] and individual[9]>=individual[15] and individual[15]>=individual[6]:
		print("#Denoising->Normalization ->PeakPicking ->Precursor Removal->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)			
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked


	elif individual[0]>=individual[15] and individual[15]>=individual[9] and individual[9]>=individual[6]:
		print("#Denoising->PeakPicking ->Normalization ->Precursor Removal->PickFiltering")
			
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		
		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)	

		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[0]>=individual[15] and individual[15]>=individual[6] and individual[6]>=individual[9]:
		print("#Denoising->PeakPicking ->Precursor Removal-> Normalization ->PickFiltering")
		
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		
		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					#filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists		

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[6]>=individual[0] and individual[0]>=individual[9] and individual[9]>=individual[15]:
		print("#Precursor Removal-> Denoising->Normalization ->PeakPicking ->PickFiltering")
			
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum	
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		

		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum		

		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[6]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[9]:
		print("#Precursor Removal-> Denoising->PeakPicking ->Normalization ->PickFiltering")
			
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)


		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists			

		

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked

	elif individual[9]>=individual[0] and individual[0]>=individual[6] and individual[6]>=individual[15]:
		print("#Normalization-> Denoising->Precursor Removal -> PeakPicking ->PickFiltering")
				

		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum			

		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[9]>=individual[0] and individual[0]>=individual[15] and individual[15]>=individual[6]:
		print("#Normalization-> Denoising-> PeakPicking ->Precursor Removal -> PickFiltering")
		
		
			
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum		

		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			if individual[16]<1:
				for j in range(0,parsed_data.size()):
					if parsed_data[j].getMSLevel() == 1:		
						continue
					filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
			else:
				for j in range(0,parsed_data_picked.size()):
					if parsed_data_picked[j].getMSLevel() == 1:		
						continue
					filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Listsering peaklists
		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	elif individual[6]>=individual[9] and individual[9]>=individual[0]:
		print("#Precursor Removal-> Normalization ->Denoising ->PeakPicking ->PickFiltering")
					
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum		
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked

	elif individual[9]>=individual[6] and individual[6]>=individual[0]:
		print("#Normalization-> Precursor Removal->Denoising ->PeakPicking ->PickFiltering")
			
		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)			

		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked
	
	else:
		print("#default")
		
		#Denoising		
		if individual[1]<1:
			#no normalization
			nothing=0
		elif individual[1]<2:
			#Gaussian Normalization
			filter = GaussFilter()
			param=pyopenms.Param()
			param.setValue("gaussian_width",individual[2],"Gaussian width")
			param.setValue("ppm_tolerance",individual[3],"Gaussian width, depending on the m/z position. The higher the value, the wider the peak and therefore the wider the gaussian.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
		else:
			#Savitzky Golay Normalization
			filter=SavitzkyGolayFilter()
			param=pyopenms.Param()
			param.setValue("frame_length",int(round(individual[4])),"The number of subsequent data points used for smoothing. This number has to be uneven. If it is not, 1 will be added.")
			if (int(round(individual[5]))>=int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1	
			param.setValue("polynomial_order",int(round(individual[5])), "Order or the polynomial that is fitted.")
			filter.setParameters(param)
			filter.filterExperiment(parsed_data)
	
		#Precursor Removal
		if individual[7]<1:
			#no precursor removal
			nothing=0
		else:
			#ParentPeakMower Method
			param=pyopenms.Param()
			param.setValue("window_size",float(individual[8]),"The size of the m/z window where the peaks are removed, +/- window_size.")
			filter=ParentPeakMower()
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum

		#Normalization:
		if individual[10]<1:
			#No normalization
			nothing=0
		elif individual[10]<2:
			#Normalizer
			filter=Normalizer()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<3:
			#Scaler
			filter=Scaler()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		elif individual[10]<4:
			#SqrtMower
			filter=SqrtMower()
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum
		else:
			#BernNorm
			filter=BernNorm()
			param=pyopenms.Param()
			param.setValue("threshold",individual[11],"Threshold of the Bern et al. normalization")
			filter.setParameters(param)
			for j in range(0,parsed_data.size()):
				if parsed_data[j].getMSLevel() == 1:		
					continue
				filter.filterSpectrum(parsed_data[j]) #for filtering Spectrum		
		#PeakPicking
		if individual[16]<1:
			#Default PeakPickerHiRes is applied
			nothing=0
		elif individual[16]<2:
			#PeakPickerHiRes
			pp=PeakPickerHiRes()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[17],"Minimal signal-to-noise ratio for a peak to be picked (0.0 disables SNT estimation!)")
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)
		else:
			#PeakPickerCWT
			pp=PeakPickerCWT()
			param=pyopenms.Param()
			param.setValue("signal_to_noise",individual[21],"Minimal signal to noise ratio for a peak to be picked.")
			param.setValue("peak_width",float(individual[22]),"Approximate fwhm of the peaks.")
			if individual[23]>0.5:
				param.setValue("estimate_peak_width","true","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")
			else:
				param.setValue("estimate_peak_width","false","Flag if the average peak width shall be estimated. Attention: when this flag is set, the peak_width is ignored.")			
			pp.setParameters(param)
			pp.pickExperiment(parsed_data,parsed_data_picked)

		#Peak Filtering
		if individual[12]<1 or individual[16]<1:
			#No Peak Filtering
			nothing=0
		elif individual[12]<2:
			#Threshold Mower 
			filter=ThresholdMower()
			param=pyopenms.Param()
			param.setValue("threshold",individual[13],"Intensity threshold, peaks below this threshold are discarded")
			filter.setParameters(param)
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		else:
			#Nlargest
			filter=NLargest()
			param=pyopenms.Param()
			param.setValue("n",int(round(individual[14])),"The number of peaks to keep")
			filter.setParameters(param)
			nr_ms2_spectra=0
			for j in range(0,parsed_data_picked.size()):
				if parsed_data_picked[j].getMSLevel() == 1:		
					continue
				filter.filterPeakSpectrum(parsed_data_picked[j]) #for filtering Peak Lists
		if individual[16]<1:
			return parsed_data
		else:		
			return parsed_data_picked

def quantification_and_spectra_alignment(input_folder,max_allowed_variance_in_peak_number,minimum_similarity_score,retention_time_tolerance,quantified_spectra_thres,abs_sim_thres,min_first_to_second_score_distance, missing_values_filename, output_folder):
	missing_values_fid      = open(missing_values_filename,'r')
	#parse the missing values filename
	missing_values_list=list()
	missing_values_raws=0;
	print('missing values filename is being parsed')
	for line1 in missing_values_fid:
		words=line1.split('\t')
		missing_values_list.append([])
		for i in range(len(words)):
			missing_values_list[missing_values_raws].append(words[i].strip())
		missing_values_raws=missing_values_raws+1;	
	spectra_peaks=list()
	exp_spectra_list = MSExperiment()
	"""#Table list with the quantitative information of the unified spectra list among all experiments last row contains quantitative classes"""
	#The dot product of peaks in the peak list with themeselves
	dot_product=list()
	intensities_sum=list()
	spectra_quantitative_dataset=list()
	maximum_number_of_spectra=30000
	#Table-list with the mapping of individual experiment's spectra on the unified spectra list	
	spectra_experiment_mapping=[[] for j in range(maximum_number_of_spectra)]
	#similarities mapping	
	spectra_experiment_similarities_mapping=[[] for j in range(maximum_number_of_spectra)]
	#Number of spectra being added in the unified lists so far
	num_of_spectra=0
	#parse all mzML peak-list files to construct the output tables
	if os.path.exists(output_folder+"filtered_files/"):
		number_of_parsed_files=0
		filenames=list()
		#For every file in the filtered_files folder
		for filename in os.listdir(output_folder+"filtered_files/"):
			#store the filename to the filesnames list
			filenames.append(filename)
		filenames.sort()
		for f in range(len(filenames)):
			print(filenames[f])
			missing_values_raw=list()
			for i in range(len(missing_values_list)):
				if filename[f][9:]==missing_values_list[i][0]:
					for j in xrange(1, len(missing_values_list[0])):
						missing_values_raw.append(missing_values_list[i][j])
			exp = MSExperiment()
			#create pyopenms File object and load a peak list in the format of a mzML file
			file = pyopenms.FileHandler()
			file.loadExperiment(output_folder+"filtered_files/"+filenames[f], exp)
			parsed_file=exp
			#store the filename to the filesnames list
			
			#If this is the first file being parsed then store all its quantified spectra in the unified spectra peaks list, in the unified spectra quantification list and update the mapping list.
			if number_of_parsed_files==0:
				#This variable measures the number of quantified spectra included so far in the spectra list.
				quantified=0
				number_of_intensities=0
				#mean_intensities=0.0
				for k in range(parsed_file.size()):
					if k%1000==0:
						print(k)
					#For MS spectra just state in the mapping list that they are not included in the spectra peak list
					if parsed_file[k].getMSLevel() == 1:
						#spectra_experiment_mapping.append([])
						spectra_experiment_mapping[k].append(-1)
						spectra_experiment_similarities_mapping[k].append(-1)		
						continue
					peaks = parsed_file[k].get_peaks()
					if len(peaks)==0:
						mz_values_old=list()
						intensities_old=list()
					else:					
						mz_values_old, intensities_old = zip(*peaks)
						#mz_values_old=list(peaks[0])
						#intensities_old=list(peaks[1])
						mz_values = [float(v) for v in mz_values_old]
						intensities = [float(v) for v in intensities_old]					
					
					#get quantification information					
					exact_reporter_ions=0
					dictionary_result=count_intervals(mz_values, [125.5, 126.5,127.5, 128.5, 129.5, 130.5, 131.5])
					number_of_non_missing_channels=0
					for chan in range(len(missing_values_raw)):
						if(missing_values_raw[chan]=='0'):
							number_of_non_missing_channels+=1;
					if dictionary_result[126.5]>=1 and missing_values_raw[0]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if dictionary_result[127.5]>=1 and missing_values_raw[1]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[128.5]>=1 and missing_values_raw[2]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[129.5]>=1 and missing_values_raw[3]=='0':
						exact_reporter_ions=exact_reporter_ions+1					
					if dictionary_result[130.5]>=1 and missing_values_raw[4]=='0':
						exact_reporter_ions=exact_reporter_ions+1				
					if dictionary_result[131.5]>=1 and missing_values_raw[5]=='0':
						exact_reporter_ions=exact_reporter_ions+1
					if exact_reporter_ions>=quantified_spectra_thres*number_of_non_missing_channels:
						quantified=quantified+1
						#add a new row in unified spectra peak lists and add the peak list of this spectrum
						spectra_peaks.append([])
						exp_spectra_list.addSpectrum(parsed_file[k])
						
						spectra_peaks[num_of_spectra].append(peaks)
						
						d_product=0.0
						int_sum=0.0
						#precalculate dot products and intensities sum for spectra comparison step
						for m in range(len(mz_values)):
							if mz_values[m]<126 or mz_values[m]>132:
								d_product = d_product + intensities[m]*intensities[m]
								int_sum=int_sum+intensities[m]
						dot_product.append(d_product)
						intensities_sum.append(int_sum)
						#Update Spectra Experiment Mapping
						spectra_experiment_mapping[k].append(num_of_spectra)
						spectra_experiment_similarities_mapping[k].append(1)
						#Add a new row at spectra quantitative dataset
						spectra_quantitative_dataset.append([])
						#find maximum peak in all channels and store it to the unified spectra quantification list
						for reporter in range(6):
							if dictionary_result[126.5+reporter]>=1  and (missing_values_raw[reporter] == 0):
								reporter_mz_values=[mz for mz in mz_values if (125.5+reporter)<=mz<=(126.5+reporter) ]
								quantitative_values=[[] for q in range(len(reporter_mz_values))]
								for mz in range(len(reporter_mz_values)):
									quantitative_values[mz]=intensities[mz_values.index(reporter_mz_values[mz])]
								spectra_quantitative_dataset[num_of_spectra].append(max(quantitative_values))
							else:
								#If no peak is found for this channel add -1
								spectra_quantitative_dataset[num_of_spectra].append(-1)
								
								
						num_of_spectra=num_of_spectra+1
					else:
						#Update the mapping list withe -1 in this position as this spectrum is not added to spectra list.						
						spectra_experiment_mapping[k].append(-1)
						spectra_experiment_similarities_mapping[k].append(-1)
				for k in range((30000-parsed_file.size())):
					spectra_experiment_mapping[parsed_file.size()+k].append(-1)
					spectra_experiment_similarities_mapping[parsed_file.size()+k].append(-1)
			#If this is not the first file being parsed
			else:
				#This variable measures the number of quantified spectra included so far in the spectra list.
				quantified=0
				spectra_found=0
				spectra_found_twice=0
				not_found=0
				#mean_intensities=mean_intensities/number_of_intensities
				c=(minimum_similarity_score*2)/(1000.0)
				association=[-1]*parsed_file.size()
				similarity_list=[0]*parsed_file.size()
				newly_added=0
				arg_fabs=math.fabs
				spectra_found_list=[0]*num_of_spectra
				new_peaks=list();
				new_dot_products=list();
				new_intensities_sum=list();
				new_dictionary_results=list();
				new_intensities=list();
				new_mz_values=list();
				new_positions=list();
				#create a list of Searchparams objects
				obj_searchparams = []
				#print parsed_file.size()
				start_time2 = time.time()
				exp_spectra_list_rt = list()
				for sp in xrange(0,num_of_spectra-1):
					exp_spectra_list_rt.append(exp_spectra_list[sp].getRT())				
				for k in range(parsed_file.size()):
					parsed_file_k_MSLevel = parsed_file[k].getMSLevel()
					parsed_file_k_peaks = parsed_file[k].get_peaks()
					parsed_file_k_RT = parsed_file[k].getRT()
					obj_searchparams.append(Searchparams(k,
						parsed_file_k_MSLevel,
						parsed_file_k_peaks,
						parsed_file_k_RT,
						num_of_spectra,
						retention_time_tolerance, 
						max_allowed_variance_in_peak_number,
						intensities_sum,
						c,
						dot_product,
						minimum_similarity_score,
						spectra_found_list,
						arg_fabs,
						exp_spectra_list_rt,
						spectra_peaks,
						missing_values_raw,
						quantified_spectra_thres,abs_sim_thres,
						min_first_to_second_score_distance))
				#pool = mp.Pool(processes=mp.cpu_count())
				print "Fill objects array at: {} secs".format(time.time() - start_time2)
				start_time2 = time.time()
				processes = 5
				#calculate chunksize
				poolChunksize = 1000
				print "Start parallel jobs in {} processorsm with chunksize: {}".format(processes,poolChunksize)
				pool = mp.Pool(processes)				
				results = pool.map(searchparams_spectrum_in_unified_list, obj_searchparams, chunksize=poolChunksize)
				pool.close()
				pool.join()
				print "Parallel job finished in: {} secs".format(time.time() - start_time2)
				#results = [[] for i in range(len(parsed_file.size()))]
				#for k in range(parsed_file.size()):
				#	results[k]=search_spectrum_in_unified_list(k,parsed_file,num_of_spectra,retention_time_tolerance, max_allowed_variance_in_peak_number,intensities_sum,c,dot_product,minimum_similarity_score,spectra_found_list,arg_fabs,exp_spectra_list,spectra_peaks)
				#	results[k]=searchparams_spectrum_in_unified_list(obj_searchparams)
				for k in range(parsed_file.size()):	
					if results[k]['mapping']==-2:
						spectra_experiment_mapping[k].append(-1)
						spectra_experiment_similarities_mapping[k].append(-1)
					elif results[k]['spectra_found']==1:
						spectra_found=spectra_found+1
						association[k]=results[k]['association']
						spectra_found_list[results[k]['association']]=spectra_found_list[results[k]['association']]+1
						similarity_list[k]=results[k]['maximum']
					else:
						not_found=not_found+1
						#spectra_experiment_mapping[k].append(num_of_spectra+newly_added)
						exp_spectra_list.addSpectrum(parsed_file[k])
						newly_added=newly_added+1;
						new_positions.append(k)
						spectra_experiment_similarities_mapping[k].append(1)
						new_peaks.append(results[k]['peaks'])
						new_dot_products.append(results[k]['d_product2'])
						new_intensities_sum.append(results[k]['int_sum2'])
						new_dictionary_results.append(results[k]['dictionary_result'])
						new_intensities.append(results[k]['intensities'])
						new_mz_values.append(results[k]['mz_values'])
				for n in range(len(new_peaks)):
					spectra_experiment_mapping[new_positions[n]].append(num_of_spectra+n)
					#Insert new spectra in spectra peaks list
					spectra_peaks.append([])
					spectra_peaks[num_of_spectra+n].append(new_peaks[n])
					dot_product.append(new_dot_products[n])
					intensities_sum.append(new_intensities_sum[n])
					#Add a new row at spectra quantitative dataset
					spectra_quantitative_dataset.append([])
					for rep in range(number_of_parsed_files):
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)
						spectra_quantitative_dataset[num_of_spectra+n].append(-1)

					#find maximum peak in all channels and store it to the unified spectra quantification list
					for reporter in range(6):
						if new_dictionary_results[n][126.5+reporter]>=1  and (missing_values_raw[reporter] == 0):
							reporter_mz_values=[mz for mz in new_mz_values[n] if (125.5+reporter)<=mz<(126.5+reporter)]
							quantitative_values=[[] for q in range(len(reporter_mz_values))]
							for mz in range(len(reporter_mz_values)):
								quantitative_values[mz]=new_intensities[n][new_mz_values[n].index(reporter_mz_values[mz])]
							spectra_quantitative_dataset[num_of_spectra+n].append(max(quantitative_values))
						else:
							#If no peak is found for this channel add -1
							spectra_quantitative_dataset[num_of_spectra+n].append(-1)

													
		
				#print("Spectra found="+str(spectra_found))
				#print("Spectra found twice="+str(spectra_found_twice))
				print("not found:"+str(not_found))
				print("found:"+str(spectra_found))
				n_spec=num_of_spectra				
				for k in range(n_spec):
					positions=[val for val, x in enumerate(association) if x == k]
					if len(positions)>=1:
						max1=0
						position_max=0
						for p in range(len(positions)):
							if similarity_list[positions[p]]>max1:
								max1=similarity_list[positions[p]]
								position_max=p
								#Update Spectra Experiment Mappings
							spectra_experiment_mapping[positions[p]].append(k)
							spectra_experiment_similarities_mapping[positions[p]].append(max1)
						peaks = parsed_file[positions[position_max]].get_peaks()
						if len(peaks)==0:
							mz_values_old=list()
							intensities_old=list()
						else:					
							mz_values_old, intensities_old = zip(*peaks)
							#mz_values_old=list(peaks[0])
							#intensities_old=list(peaks[1])
							mz_values = [float(v) for v in mz_values_old]
							intensities = [float(v) for v in intensities_old]				
						#get quantification information					
						exact_reporter_ions=0
						dictionary_result=count_intervals(mz_values, [125.5, 126.5,127.5, 128.5, 129.5, 130.5, 131.5])
						number_of_non_missing_channels=0
						for chan in range(len(missing_values_raw)):
							if(missing_values_raw[chan]=='0'):
								number_of_non_missing_channels+=1;
						if dictionary_result[126.5]>=1 and missing_values_raw[0]=='0':
							exact_reporter_ions=exact_reporter_ions+1
						if dictionary_result[127.5]>=1 and missing_values_raw[1]=='0':
							exact_reporter_ions=exact_reporter_ions+1					
						if dictionary_result[128.5]>=1 and missing_values_raw[2]=='0':
							exact_reporter_ions=exact_reporter_ions+1					
						if dictionary_result[129.5]>=1 and missing_values_raw[3]=='0':
							exact_reporter_ions=exact_reporter_ions+1					
						if dictionary_result[130.5]>=1 and missing_values_raw[4]=='0':
							exact_reporter_ions=exact_reporter_ions+1				
						if dictionary_result[131.5]>=1 and missing_values_raw[5]=='0':
							exact_reporter_ions=exact_reporter_ions+1
						
						#find maximum peak in all channels and store it to the unified spectra quantification list
						for reporter in range(6):
							if dictionary_result[126.5+reporter]>=1  and (missing_values_raw[reporter] == 0):
								reporter_mz_values=[mz for mz in mz_values if (125.5+reporter)<=mz<=(126.5+reporter)]
								quantitative_values=list()
								quantitative_values=[[] for q in range(len(reporter_mz_values))]
								for mz in range(len(reporter_mz_values)):
									quantitative_values[mz]=intensities[mz_values.index(reporter_mz_values[mz])]
								if len(quantitative_values)==0:
									spectra_quantitative_dataset[k].append(-1)
								else:
									#This could be changed to average!!! (TODO)
									spectra_quantitative_dataset[k].append(max(quantitative_values))
							else:
								#If no peak is found for this channel add -1
								spectra_quantitative_dataset[k].append(-1)

					else:
						#This spectrum was not associated with no other spectrum
						#update the quantitiative dataset table
						for reporter in range(6):
							spectra_quantitative_dataset[k].append(-1)
				num_of_spectra=num_of_spectra+newly_added
				for k in range((30000-parsed_file.size())):
					spectra_experiment_mapping[parsed_file.size()+k].append(-1)
					spectra_experiment_similarities_mapping[parsed_file.size()+k].append(-1)
			number_of_parsed_files=number_of_parsed_files+1
		#write mapping
		mapping_spectra_file_fid = open("mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
		for i in range(len(filenames)):
			mapping_spectra_file_fid.write(filenames[i])
			if i<len(filenames)-1:
				mapping_spectra_file_fid.write("\t")
		mapping_spectra_file_fid.write("\n")
		for i in range(len(spectra_experiment_mapping)):
			for j in range(len(spectra_experiment_mapping[i])): 
				if j<len(spectra_experiment_mapping[i])-1:
					
					mapping_spectra_file_fid.write(str(spectra_experiment_mapping[i][j]))
					mapping_spectra_file_fid.write("\t")
				else:
					mapping_spectra_file_fid.write(str(spectra_experiment_mapping[i][j]))
			if i<len(spectra_experiment_mapping)-1:
				mapping_spectra_file_fid.write("\n")	
		mapping_spectra_file_fid.close()
		#write similarities mapping
		similarities_mapping_spectra_file_fid = open("similarities_mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
		for i in range(len(filenames)):
			similarities_mapping_spectra_file_fid.write(filenames[i])
			if i<len(filenames)-1:
				similarities_mapping_spectra_file_fid.write("\t")
		similarities_mapping_spectra_file_fid.write("\n")
		for i in range(len(spectra_experiment_similarities_mapping)):
			for j in range(len(spectra_experiment_similarities_mapping[i])): 
				if j<len(spectra_experiment_similarities_mapping[i])-1:
					
					similarities_mapping_spectra_file_fid.write(str(spectra_experiment_similarities_mapping[i][j]))
					similarities_mapping_spectra_file_fid.write("\t")
				else:
					similarities_mapping_spectra_file_fid.write(str(spectra_experiment_similarities_mapping[i][j]))
			if i<len(spectra_experiment_mapping)-1:
				similarities_mapping_spectra_file_fid.write("\n")	
		similarities_mapping_spectra_file_fid.close()
		#write quantitative values
		quantitative_values_fid = open("quantitative_values"+str(time.strftime("%Y_%m_%d"))+".txt","w")
		for i in range(len(spectra_quantitative_dataset)):
			#print("i="+str(i))
			for j in range(len(spectra_quantitative_dataset[i])): 
				if j<len(spectra_quantitative_dataset[i])-1:
					
					#print("j="+str(j))
					quantitative_values_fid.write(str(spectra_quantitative_dataset[i][j]))
					quantitative_values_fid.write("\t")
				else:
					quantitative_values_fid.write(str(spectra_quantitative_dataset[i][j]))
			if i<len(spectra_quantitative_dataset)-1:
				quantitative_values_fid.write("\n")	
		
		quantitative_values_fid.close()
		file = MzMLFile()	
		file.store("unified_spectrum_list_"+str(time.strftime("%Y_%m_%d"))+".mzML",exp_spectra_list)		
		print(str(number_of_parsed_files)+" mzML files were successfully filtered and stored")

def searchparams_spectrum_in_unified_list(searchparams):
    #print searchparams.parsed_file_k_MSLevel
    return searchparams.search_spectrum_in_unified_list()

class Searchparams():

	def __init__(self, k, parsed_file_k_MSLevel,parsed_file_k_peaks,parsed_file_k_RT, num_of_spectra, retention_time_tolerance, max_allowed_variance_in_peak_number, intensities_sum, c, dot_product, minimum_similarity_score, spectra_found_list, arg_fabs, exp_spectra_list_rt, spectra_peaks, missing_values_raw,quantified_spectra_thres,abs_sim_thres, min_first_to_second_score_distance):
		self.k = k
		#self.parsed_file = parsed_file
		self.num_of_spectra = num_of_spectra
		self.retention_time_tolerance = retention_time_tolerance
		self.max_allowed_variance_in_peak_number = max_allowed_variance_in_peak_number
		self.intensities_sum = intensities_sum
		self.c = c
		self.dot_product = dot_product
		self.minimum_similarity_score = minimum_similarity_score
		self.spectra_found_list = spectra_found_list
		self.arg_fabs = arg_fabs
		self.exp_spectra_list_rt = exp_spectra_list_rt
		self.spectra_peaks = spectra_peaks
		self.parsed_file_k_MSLevel = parsed_file_k_MSLevel
		self.parsed_file_k_peaks = parsed_file_k_peaks
		self.parsed_file_k_RT = parsed_file_k_RT
		self.missing_values_raw=missing_values_raw
		self.quantified_spectra_thres=quantified_spectra_thres
		self.abs_sim_thres=abs_sim_thres
		self.min_first_to_second_score_distance=min_first_to_second_score_distance

	def search_spectrum_in_unified_list(self):
		try:
			
			#For MS spectra just state in the mapping list that they are not included in the spectra peak list
			if self.parsed_file_k_MSLevel == 1:
				#print "Job {} has finished at 1".format(self.k)
				return {'mapping':-2, 'similarity':-1, 'k':self.k,'spectra_found':0,'association':0,'peaks':0,'d_product2':0,'int_sum2':0,'mz_values':0,'dictionary_result':0,'intensities':0,'maximum':0}
			peaks = self.parsed_file_k_peaks
			if len(peaks)==0:
				mz_values_old=list()
				intensities_old=list()
			else:					
				mz_values_old, intensities_old = zip(*peaks)
				#mz_values_old=list(peaks[0])
				#intensities_old=list(peaks[1])
			mz_values = [float(v) for v in mz_values_old]
			intensities = [float(v) for v in intensities_old]			
			#get quantification information					
			exact_reporter_ions=0
			dictionary_result=count_intervals(mz_values, [125.5, 126.5,127.5, 128.5, 129.5, 130.5, 131.5])
			number_of_non_missing_channels=0
			for chan in range(len(self.missing_values_raw)):
				if(self.missing_values_raw[chan]=='0'):
					number_of_non_missing_channels+=1;
			if dictionary_result[126.5]>=1 and self.missing_values_raw[0]=='0':
				exact_reporter_ions=exact_reporter_ions+1
			if dictionary_result[127.5]>=1 and self.missing_values_raw[1]=='0':
				exact_reporter_ions=exact_reporter_ions+1					
			if dictionary_result[128.5]>=1 and self.missing_values_raw[2]=='0':
				exact_reporter_ions=exact_reporter_ions+1					
			if dictionary_result[129.5]>=1 and self.missing_values_raw[3]=='0':
				exact_reporter_ions=exact_reporter_ions+1					
			if dictionary_result[130.5]>=1 and self.missing_values_raw[4]=='0':
				exact_reporter_ions=exact_reporter_ions+1				
			if dictionary_result[131.5]>=1 and self.missing_values_raw[5]=='0':
				exact_reporter_ions=exact_reporter_ions+1
			#for spectra with more than three channels being quantified update all three files					
			#print(exact_reporter_ions)
			if exact_reporter_ions>=self.quantified_spectra_thres*number_of_non_missing_channels:
				found=-1
				minimum_score=1
				flag_once=0
				d_product2=0.0
				int_sum2=0.0
				#precalculate dot products and intentisites sum
				for m in range(len(mz_values)):
					if mz_values[m]<126 or mz_values[m]>132:
						d_product2 = d_product2 + intensities[m]*intensities[m]
						int_sum2=int_sum2+intensities[m]
				similarity_scores=[0]*self.num_of_spectra
				disimilar=0	
				len1=len(peaks)
				#print(len1)
			
				RT_time_current_spectra=self.parsed_file_k_RT
				if len1>8:
					for sp in xrange(0,self.num_of_spectra-1):
						#if sp>100 and sp<150 and k%500==0:
						#	print("Number of peaks:"+str(len(spectra_peaks[sp][0])))
						#	print("Peaks"+str(spectra_peaks[sp][0]))
						if self.arg_fabs(self.exp_spectra_list_rt[sp]-RT_time_current_spectra)<self.retention_time_tolerance:
							if self.spectra_found_list[sp]<1:
								len2=len(self.spectra_peaks[sp][0])
								if len2>8 and self.arg_fabs(len1-len2)<self.max_allowed_variance_in_peak_number:
									expected_product=self.c*self.intensities_sum[sp]*int_sum2
									sim_score=similarity_measurement(mz_values,intensities, self.spectra_peaks[sp][0],self.dot_product[sp],d_product2,expected_product, self.minimum_similarity_score)
									similarity_scores[sp]=sim_score
									if sim_score>self.abs_sim_thres:
										break
				max_list=nlargest(2, similarity_scores)
				maximum=max_list[0]
				second_maximum=max_list[1]
				if self.k%100==0:
					print("Maximum Similarity Score="+str(maximum))
					print("Second Maximum Similarity Score="+str(second_maximum))
				#print("Maximum Similarity Score="+str(maximum))
				if ((maximum>0.7) or (maximum>0.2 and (maximum-second_maximum)>self.min_first_to_second_score_distance)):
					return {'mapping':0, 'similarity':0, 'k':self.k,'spectra_found':1,'association':similarity_scores.index(maximum),'peaks':peaks,'d_product2':d_product2,'int_sum2':int_sum2,'mz_values':mz_values,'dictionary_result':dictionary_result,'intensities':intensities,'maximum':maximum}		
					#spectra_found=spectra_found+1
					#association[k]=similarity_scores.index(maximum)
					#spectra_found_list[similarity_scores.index(maximum)]=spectra_found_list[similarity_scores.index(maximum)]+1
					#similarity_list[k]=maximum
				else:
					#print "Job {} has finished at 2".format(self.k)
					if self.num_of_spectra>60000:
						return {'mapping':-2, 'similarity':-1, 'k':self.k,'spectra_found':0,'association':0,'peaks':0,'d_product2':0,'int_sum2':0,'mz_values':0,'dictionary_result':0,'intensities':0,'maximum':0}
					else:
						return {'mapping':-1, 'similarity':-1, 'k':self.k,'spectra_found':0,'association':0,'peaks':peaks,'d_product2':d_product2,'int_sum2':int_sum2,'mz_values':mz_values,'dictionary_result':dictionary_result,'intensities':intensities,'maximum':0}
			else:
				#print "Job {} has finished at 3".format(self.k)
				return {'mapping':-2, 'similarity':-1, 'k':self.k,'spectra_found':0,'association':0,'peaks':0,'d_product2':0,'int_sum2':0,'mz_values':0,'dictionary_result':0,'intensities':0,'maximum':0}
		except:
			#print e
			return {}
def similarity_measurement(mz_values1,intensities1, peaks2,dot_product11,dot_product22, expected_product,mz_tolerance, arg1=math.sqrt,arg2=math.fabs):
	score=0
	mz_values_old2, intensities_old2 = zip(*peaks2)
	#mz_values_old2=list(peaks2[0])
	#intensities_old2=list(peaks2[1])
	mz_values2 = [float(v) for v in mz_values_old2]
	intensities2 = [float(v) for v in intensities_old2]
	dot_product12=0.0
	arg3=mz_values2.index
	for i in xrange(0,len(mz_values1)-1):
		if mz_values1[i]<126 or mz_values1[i]>132:		
			closest_value=takeClosest(mz_values2,mz_values1[i])
			if arg2(closest_value-mz_values1[i])<(mz_tolerance+mz_tolerance):
				dot_product12 = dot_product12 + intensities2[arg3(closest_value)]*intensities1[i]	
	score=(dot_product12-expected_product)/arg1(dot_product11*dot_product22)	
	return score

#a function that counts the elements of a list within an interval
def count_intervals(sequence, intervals):
	count = defaultdict(int)
	#intervals.sort()
	for item in sequence:
		pos = bisect_left(intervals, item)
		if pos == len(intervals):
			count[None] += 1
		else:
			count[intervals[pos]] += 1
	return count

#Assumes myList is sorted. Returns closest value to myNumber. If two numbers are equally close, return the smallest number.
def takeClosest(myList, myNumber):
	pos = bisect_left(myList, myNumber)
	if pos == 0:
		return myList[0]
	if pos == len(myList):
		return myList[-1]
	before = myList[pos - 1]
	after = myList[pos]
	if after - myNumber < myNumber - before:
		return after
	else:
		return before

if __name__ == "__main__":
	
	
	output_folder='step2/'+str(tstamp)+'_'+str(time.time())+'/'
	#Initialize random generator seed with the current local time in miliseconds
	random.seed()
	input_data_folder=sys.argv[1]
	best_solution_path=sys.argv[2]
	minimum_similarity=float(sys.argv[3])
	retention_time_tolerance=float(sys.argv[4])
	quantified_spectra_thres=float(sys.argv[5])
	max_allowed_variance_in_peak_number=int(sys.argv[6])
	abs_sim_thres=float(sys.argv[7])
	min_first_to_second_score_distance=float(sys.argv[8])
	missing_values_filename = sys.argv[9]
	apply_best_solution(input_data_folder,best_solution_path,quantified_spectra_thres, ubuntu_flag, missing_values_filename, output_folder)
	start_time = time.time()
	quantification_and_spectra_alignment(input_data_folder,max_allowed_variance_in_peak_number,minimum_similarity,retention_time_tolerance,quantified_spectra_thres,abs_sim_thres,min_first_to_second_score_distance, missing_values_filename, output_folder)
	print("--- %s seconds ---" % (time.time() - start_time))
