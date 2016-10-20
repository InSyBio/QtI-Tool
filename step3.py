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
import multiprocessing as mp
import glob
from shutil import copyfile
import platform
import numpy
import scipy.stats as st
import scipy
from collections import defaultdict
from bisect import bisect_left
from subprocess import Popen
from svmutil import *
import warnings
import logging



###
import re
import csv
import pdb

def step3_supervised_analysis(unified_list_filename,mapping_filename,quantitative_filename,labels_filename,population,generations,semi_random_initialization_probability,two_points_crossover_probability,arithmetic_crossover_probability,mutation_probability,gaussian_mutation_variance_proportion,min_values,max_values,goal_significances, database_name,spectraST_databases_filename,num_of_folds,step2_output_folder,output_folder,replicates, missing_values_filename,use_inputs_vector_filename,supervised_mode):
	#unified_list_filename:(String) it is the filename of the unified spectra list
	#mapping_filename:(String) it is the filename of the spectra mapping txt file
	#quantitative_filename: (String) it is the filename of the quantitative spectra txt file
	#labels_filename: (String) it is the filename of the labels txt file
	#population:(integer) (default: 10) it is the number of individual solutions which are evolved on parallel
	#generations: (integer) (default: 10) it is the maximum nulber of generations which we allow for the population to get evolved
	#semi_random_initialization_probability: (float) (default: 0.8) it is the proportion (multiplied with 100) of the individuals which will be initialized with the semi automatic method. The rest will be initialized using random method
	#two_points_crossover_probability: (float) (default: 0.45) it is the probability (multiplied with 100) for using two points crossover operator
	#arithmetic_crossover_probability: (float) (default: 0.45) it is the probability (multiplied with 100) for using  arithmetic crossover operator
	#mutation_probability: (float) (default: 0.05) it is the probability for applying gaussian mutation operator
	#gaussian_mutation_variance_proportion: (float) (default:0.1) it is the proportion of the abs(max_val-min_val) which is used as variance of the gaussian mutation operator
	#min_values: (list) it is a list with the minimum allowed values for the variables which should be optimized
	#max_values: (list) it is a list with the maximum allowed values for the variables which should be optimized
	#goal_significances: (list of floats) Includes the significances (values from 0 to 1) of the individual goals
	#database_name: (String) the string of the name of the fasta database to be used for scaffold
	#num_of_folds: (integer) the number of k folds in k-fold cross validation
	#replicates:(integer) the number of replicates per experiment
	#missing_values_filename: (String) the name of the missing values file
	#use_inputs_vector_filename: (String) the filename of a file were user states the number of TMT channels which are used in this dataset
	#supervised_mode: (integer) 1 for supervised learning, 0 for unsupervised learning
	#parse input files
	unified_spectra = MSExperiment()
	file = MzMLFile()
	file.load(step2_output_folder+unified_list_filename, unified_spectra)
	print("Unified Spectra file was successfully parsed!")
	mapping_spectra=list()
	mapping_fid=open(step2_output_folder+mapping_filename,"r")
	number_of_lines=0
	for line1 in mapping_fid:
		mapping_spectra.append([])
		words=line1.split("\t")
		for i in range(len(words)):
			mapping_spectra[number_of_lines].append(words[i])
		number_of_lines=number_of_lines+1
	print("Mapping file was successfully parsed!")	
	#parse use_inputs_vector_filename
	use_inputs_vector_fid=open(use_inputs_vector_filename,'r')
	use_inputs_vector=list()
	for line1 in use_inputs_vector_fid:
		words=line1.split('\t')
		for i in range(len(words)):
			use_inputs_vector.append(words[i].strip())
	number_of_tmt_channels_used=use_inputs_vector.count("1")	
	quantitative_spectra=list()
	quantitative_fid=open(step2_output_folder+quantitative_filename,"r")
	number_of_lines=0
	for line1 in quantitative_fid:
		quantitative_spectra.append([])
		words=line1.split("\t")
		for i in range(len(words)):
			if use_inputs_vector[i%6]=="1":
				quantitative_spectra[number_of_lines].append(float(words[i]))
		number_of_lines=number_of_lines+1
	print("Quantitative file was successfully parsed!")	
	labels=list()
	labels_fid=open(labels_filename,"r")
	number_of_lines=0
	filenames_in_labels=list()
	for line1 in labels_fid:
		labels.append([])
		words=line1.split("\t")
		filenames_in_labels.append(words[0])
		for i in range(1,len(words)):
			labels[number_of_lines].append(words[i].strip())
		number_of_lines=number_of_lines+1
	#parse the missing values file
	missing_values_fid2      = open(missing_values_filename,'r')
	missing_values_list=list()
	missing_values_raws=0
	for line1 in missing_values_fid2:
		words=line1.split('\t')
		missing_values_list.append([])
		for i in range(len(words)):
			missing_values_list[missing_values_raws].append(words[i].strip())
		missing_values_raws+=1	
	print("Labels file were successfully parsed!")
	number_of_samples=number_of_lines*number_of_tmt_channels_used
	number_of_experiments=number_of_lines
	print("number_of_experiments"+str(number_of_experiments))
	l=len(quantitative_spectra)
	quantitative_spectra.append([])
	#Merge labels with quantitative spectra
	for i in range(len(mapping_spectra[0])):
		ind=filenames_in_labels.index(mapping_spectra[0][i].strip())
		for j in range(6):
			if use_inputs_vector[j]=="1":
				quantitative_spectra[l].append(labels[ind][j])
	print("Quantitative values were successfully merged with Spectra labels")
	#Transpose list of lists samples rows and number of spectra +1 columns
	quantitative_spectra_dataset=map(list, zip(*quantitative_spectra))	
	used_spectra=[1]*(l)
	num_of_values_per_experiment=replicates*use_inputs_vector.count("1")	
	for i in range(l):

		not_found_experiments=0
		for j in range(((len(quantitative_spectra[0]))/(num_of_values_per_experiment))):
			found=0			
			for k in range(num_of_values_per_experiment):
				if quantitative_spectra[i][j*num_of_values_per_experiment+k]!=-1:
					found=1
			if found==0:
				not_found_experiments=not_found_experiments+1		
		if supervised_mode==1:
			if not_found_experiments>0.5*(number_of_experiments/replicates):
				used_spectra[i]=0
		else:
			if not_found_experiments>0*(number_of_experiments/replicates):
				used_spectra[i]=0	
	if supervised_mode==1:	
		print("Commonly Quantified Spectra (Present in more than 50% of the experiments)="+str(used_spectra.count(1)))
	else:
		print("Commonly Quantified Spectra(Present in 100% of the experiments)="+str(used_spectra.count(1)))
	print("Dataset was successfully created")
	number_of_used_tmt_channels=use_inputs_vector.count("1")
	#Normalization
	for i in range(l):
		for j in range(number_of_experiments):
			maximum=-1
			position=-1
			for k in range(number_of_used_tmt_channels):
				if quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]>maximum:
					maximum=quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]
					position=k
			if position!=-1:
				if float(quantitative_spectra_dataset[j*number_of_used_tmt_channels][i])==-1:
					for k in range(number_of_used_tmt_channels):
						quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]=-1
				else:				
					reference=float(quantitative_spectra_dataset[j*number_of_used_tmt_channels][i])				
					for k in range(number_of_used_tmt_channels):
						if float(quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i])!=-1:
							quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]=quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]/reference	
	quantitative_spectra_dataset_normalized=copy.deepcopy(quantitative_spectra_dataset)	
	quantitative_spectra_dataset_normalized=map(list, zip(*quantitative_spectra_dataset_normalized))
	normalized_fid = open(output_folder+"normalized_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	for i in range(len(quantitative_spectra_dataset_normalized)):
		for j in range(len(quantitative_spectra_dataset_normalized[0])):
			normalized_fid.write(str(quantitative_spectra_dataset_normalized[i][j])+"\t")
		normalized_fid.write("\n")
	normalized_fid.close()
	print("Dataset was successfully normalized!")
	#Missing Values Calculation
	averages=[0]*l
	number_of_missing_values=[0]*l
	for i in range(l):
		num_of_elements=0
		for j in range(number_of_samples):
			if(float(quantitative_spectra_dataset[j][i])!=-1):
				averages[i]=averages[i]+float(quantitative_spectra_dataset[j][i])
				num_of_elements=num_of_elements+1
			else:
				number_of_missing_values[i]=number_of_missing_values[i]+1
		if(num_of_elements==0):
			averages[i]=0
		else:
			averages[i]=averages[i]/float(num_of_elements)
	global_average=scipy.mean(averages)
	for i in range(l):
		for j in range(number_of_experiments):
			missing=1
			average=0.0
			num_of_elements=0
			for k in range(number_of_used_tmt_channels):
				if float(quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i])!=-1:
					average=average+float(quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i])
					num_of_elements=num_of_elements+1
					missing=0			
			if missing==1:
				for k in range(number_of_used_tmt_channels):
					quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]=global_average
			else:
				average=average/float(num_of_elements)				
				for k in range(number_of_used_tmt_channels):
					if quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]==-1:
						quantitative_spectra_dataset[j*number_of_used_tmt_channels+k][i]=average
	quantitative_spectra_dataset_normalized_mv_completed=copy.deepcopy(quantitative_spectra_dataset)	
	quantitative_spectra_dataset_normalized_mv_completed=map(list, zip(*quantitative_spectra_dataset_normalized_mv_completed))	
	missing_values_fid = open(output_folder+"normalized_missing_values_completed_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	for i in range(len(quantitative_spectra_dataset_normalized_mv_completed)):
		for j in range(len(quantitative_spectra_dataset_normalized_mv_completed[0])):
			missing_values_fid.write(str(quantitative_spectra_dataset_normalized_mv_completed[i][j])+"\t")
		missing_values_fid.write("\n")
	missing_values_fid.close()		
	print("Missing values were successfully imputed!")
	#Create Different Classification Problems
	categories=[quantitative_spectra[-1][i] for i in range(len(quantitative_spectra[-1]))]
	unique_categories=list(set(categories))
	classification_problems=[[] for i in range(((len(unique_categories)-2)*(len(unique_categories)-1))/2)]
	number_of_classification_problems=0
	print("Katigoria 0="+str(unique_categories[0]))	
	for i in xrange(1,len(unique_categories)-1):
		for k in xrange(i+1,len(unique_categories)):
			for j in range(len(categories)):
				if categories[j]==unique_categories[i]:
					classification_problems[number_of_classification_problems].append(1)
				elif categories[j]==unique_categories[k]:
					classification_problems[number_of_classification_problems].append(-1)
				else:
					classification_problems[number_of_classification_problems].append(0)
			number_of_classification_problems=number_of_classification_problems+1
	for i in range(len(classification_problems)):	
		print(len(classification_problems[i]))	
	print("Classification problems were created successfully!")
	#Initialize Population of Individual Solutions
	individuals=initialize(population,semi_random_initialization_probability,min_values, max_values)
	print("Initial population has been formulated.")	
	max_eval_per_generation=[0]*generations
	average_eval_per_generation=[0]*generations
	sum_ranked_eval_per_generation=[0]*generations
	average_ranked_eval_per_generation=[0]*generations
	selected_individuals=[[0 for x in range(len(min_values))] for x in range(population)] 
	temp_individuals=[[0 for x in range(len(min_values))] for x in range(population)] 
	best_solutions_fid = open(output_folder+"best_solutions_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	best_solution_fid = open(output_folder+"best_solution_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	average_performance_fid = open(output_folder+"average_performance_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	best_performance_fid=open(output_folder+"best_performance_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	final_solutions_fid = open(output_folder+"final_solutions_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	#Apply Evolutionary Process to Optimize Individual Solutions
	print("Running Optimization!")	
	premature_termination=0
	for rep in range(generations):
		print("Generation:"+str(rep+1))
		#evaluate population of solutions
		evaluation_values=evaluate_individuals(quantitative_spectra_dataset,quantitative_spectra_dataset_normalized,unified_spectra,individuals,goal_significances,database_name,spectraST_databases_filename,num_of_folds,used_spectra,classification_problems, number_of_used_tmt_channels, supervised_mode,output_folder)		
		average_performance=0
		for i in range(len(evaluation_values)):
			average_performance=average_performance+evaluation_values[i][-1]
		average_performance=average_performance/float(len(evaluation_values))
		print("Best Performance="+str(evaluation_values[0][-1]))
		print("Average Performance="+str(average_performance))	
		print("Convergence Percentage="+str(math.fabs(evaluation_values[0][-1]-average_performance)/average_performance))	
		if math.fabs(evaluation_values[0][-1]-average_performance)<0.001*average_performance:
			premature_termination=1			
			break
		#Estimate non dominated fronts
		assigned=0
		fronts=[0]*population
		front=1
		chromosomes_temp=copy.deepcopy(individuals)
		eval_temp=copy.deepcopy(evaluation_values)
		evaluation_values_temp=	copy.deepcopy(evaluation_values)
		print("evaluation_values")
		print(evaluation_values_temp)
		while assigned<population:
			number_of_solutions=len(chromosomes_temp)
			non_dominated_solutions=[0]*number_of_solutions
			index=[i for i in range(len(eval_temp[0]))]
			eval_temp_index=zip(eval_temp[0],index)
			ordered_list= sorted(range(len(eval_temp_index)), key=lambda k: eval_temp_index[k],reverse=True) 
			non_dominated_solutions[0]=ordered_list[0]			
			number_of_non_dominated_solutions=1		
			for i in range(1,number_of_solutions):
				n=0
				condition=0
				condition2=1
				while n<number_of_non_dominated_solutions and condition==0:
					solution1=[0]*(len(eval_temp)-1)
					solution2=[0]*(len(eval_temp)-1)
					for j in range(len(eval_temp)-1):
						solution1[j]=eval_temp[j][ordered_list[i]]
						solution2[j]=eval_temp[j][ordered_list[n]]					
					check=dominate(solution1,solution2)
					if check==3:
						condition=1
						condition2=0
					elif check==1:
						if number_of_non_dominated_solutions==1:
							condition=1
							non_dominated_solutions[0]=ordered_list[i]
						else:
							number_of_non_dominated_solutions=number_of_non_dominated_solutions-1	
							del non_dominated_solutions[n]
					n=n+1
				if condition2==1:
					non_dominated_solutions[number_of_non_dominated_solutions]=ordered_list[i]
					number_of_non_dominated_solutions=number_of_non_dominated_solutions+1
			sorted_non_dominated_solutions=sorted(non_dominated_solutions, reverse=True)
			for i in range(number_of_non_dominated_solutions):
				assigned=assigned+1
				fronts[sorted_non_dominated_solutions[i]]=front
				for j in range(len(eval_temp)):
					eval_temp[j][sorted_non_dominated_solutions[i]]	=0		
			front=front+1					
		print("Calculated Pareto Frontiers:")		
		print(fronts)
		evaluation_values=copy.deepcopy(evaluation_values_temp)		
		#apply selection operator: Ranked base selection is used
		#find and write to file maximum and average performances		
		max_eval=0
		max_position2=0
		for i in range(population):
			if evaluation_values[-1][i]>max_eval:
				max_eval=evaluation_values[-1][i]
				max_position2=i
		max_eval_per_generation[rep]=max_eval
		best_performance_fid.write(str(max_eval_per_generation[rep])+"\n")
		sum_eval=0
		for i in range(population):
			sum_eval=sum_eval+evaluation_values[-1][i]
		average_eval=sum_eval/population
		average_eval_per_generation[rep]=average_eval
		average_performance_fid.write(str(average_eval_per_generation[rep])+"\n")
		#Tune fitness values by locating and using solution niches
		sigma_share=0.5/(float(len(individuals[0]))**(0.1))
		for i in range(1,max(fronts)+1):
			ind=[y for y,x in enumerate(fronts) if x == i]
			#Calculate max values per goal per pareto frontier			
			max_significances=[0]*(len(evaluation_values)-1)
			for j in range(len(ind)):
				for goal in range(len(evaluation_values)-1):
					if evaluation_values[goal][ind[j]]>=max_significances[goal]:
						max_significances[goal]=evaluation_values[goal][ind[j]]			
			for j in range(len(ind)):
				m=0
				for k in range(len(ind)):
					d=0
					for gene in range(len(individuals[0])):
						d=d+((individuals[ind[j]][gene]-individuals[ind[k]][gene])/float(max_values[gene]-min_values[gene]))**2
					d=math.sqrt(d/(len(individuals[0])))
					if d<=sigma_share:
						m=m+(1-((d/float(sigma_share))**2))
				if m==0:
					m=1
				print("m="+str(m))				
				for goal in range(len(evaluation_values)-1):
					evaluation_values[goal][ind[j]]=float(max_significances[goal])/m
		print("evaluation_values=")		
		print(evaluation_values)
		for i in range(len(evaluation_values[0])):
			evaluation_values[-1][i]=0
			for j in range(len(evaluation_values)-1):
				evaluation_values[-1][i]=evaluation_values[-1][i]+evaluation_values[j][i]
			evaluation_values[-1][i]=evaluation_values[-1][i]/float(len(evaluation_values)-1)
		max_eval=0
		max_position=0
		for i in range(population):
			if evaluation_values[-1][i]>max_eval:
				max_eval=evaluation_values[-1][i]
				max_position=i
		sum_ranked=0
		for i in range(population):
			sum_ranked=sum_ranked+evaluation_values[-1][i]
		
		sum_ranked_eval_per_generation[rep]=sum_ranked
		average_ranked_eval_per_generation[rep]=sum_ranked/population
	
		sum_prop=[0]*(population+1)
		for i in range(1, population+1):
			sum_prop[i]=sum_prop[i-1]+evaluation_values[-1][i-1]/float(sum_ranked_eval_per_generation[rep]) ##Check this out
		for i in range(1, population):
			random_number=random.uniform(0,1)
			for j in range(0,population):
				if random_number>=sum_prop[j] and random_number<sum_prop[j+1]:
					for k in range(0,len(individuals[0])):					
						selected_individuals[i][k]=individuals[j][k]	
			
		for k in range(len(individuals[0])):		
			selected_individuals[0][k]=individuals[max_position2][k]
		for gen in range(len(individuals[max_position2])):
			best_solutions_fid.write(str(individuals[max_position2][gen])+"\t")
		best_solutions_fid.write("\n")
		#apply crossover operator
		for i in range(1,population-1,2):
			if random_number<two_points_crossover_probability:
				#print("Two Point Crossover")
				cross_point1=0
				cross_point2=0
				while cross_point1==cross_point2:
					cross_point1=math.ceil((len(individuals[0])-4)*random.uniform(0,1))
					if cross_point1<math.floor((2*len(individuals[0])-1)/3):
						width=math.ceil(random.uniform(0,1)*(math.floor(len(individuals[0])-1)/3 -2))
						cross_point2=cross_point1+width
					else:
						width=math.ceil(random.uniform(0,1))*(math.floor(len(individuals[0])/3 -1)-2-(cross_point1-math.floor(2*len(individuals[0])/3)))
						cross_point2=cross_point1+width
				if cross_point1>cross_point2:
					temp_cross_point=cross_point1
					cross_point1=cross_point2
					cross_point2=temp_cross_point
				width=int(width)
				cross_point1=int(cross_point1)
				cross_point2=int(cross_point2)
				for j in range(cross_point1,cross_point2+1):
					temp_individuals[i][j]=selected_individuals[i+1][j]
				for j in range(cross_point1,cross_point2+1):
					selected_individuals[i+1][j]=selected_individuals[i][j]
				for j in range(cross_point1,cross_point2+1):
					selected_individuals[i][j]=temp_individuals[i][j]

				
			elif random_number>=two_points_crossover_probability and random_number<(two_points_crossover_probability+arithmetic_crossover_probability):
				alpha=random.uniform(0,1)			
				for j in range(0,len(individuals[0])):
					temp_individuals[i][j]=alpha*selected_individuals[i][j]+(1-alpha)*selected_individuals[i+1][j]
					temp_individuals[i+1][j]=(1-alpha)*selected_individuals[i][j]+(alpha)*selected_individuals[i+1][j]
				for k in range(0,len(individuals[0])):				
					selected_individuals[i][k]=temp_individuals[i][k]
					selected_individuals[i+1][k]=temp_individuals[i+1][k]			
		#apply mutation operator
		for i in range(1,population):
			for j in range(0,len(individuals[0])):
				random_number=random.uniform(0,1)
				if random_number<mutation_probability:
					selected_individuals[i][j]=selected_individuals[i][j]+random.gauss(0,0.1*(max_values[j]-min_values[j]))
			#Correct values out of boundaries
			for j in range(0, len(min_values)):
				if selected_individuals[i][j]<min_values[j]:
					selected_individuals[i][j]=min_values[j]
				if selected_individuals[i][j]>max_values[j]:
					selected_individuals[i][j]=max_values[j]
		#update the population with the offspings
		individuals=copy.deepcopy(selected_individuals)
	for gen in range(len(individuals[0])):
		if gen<	len(individuals[0])-1:
			best_solution_fid.write(str(individuals[0][gen])+"\t")
		else:
			best_solution_fid.write(str(individuals[0][gen]))
	best_solutions_fid.write("\n")
	#Write all final solutions in a file
	for memb in range(population):
		for gen in range(len(individuals[memb])):
			if gen<len(individuals[memb])-1:
				final_solutions_fid.write(str(individuals[memb][gen])+"\t")
			else:
				final_solutions_fid.write(str(individuals[memb][gen]))
		if memb<population-1:
			final_solutions_fid.write("\n")	
	l=len(individuals[0])-16
	selected_first=[0]*l
	for spectra in range(l):
		if used_spectra[spectra]==0:
			continue
		
		data1=list()
		for j in range(len(quantitative_spectra_dataset)):
			if (float(quantitative_spectra_dataset_normalized[spectra][j])!=-1):						
					data1.append(quantitative_spectra_dataset[j][spectra])				
		if len(data1)>1:
			if numpy.var(data1)>0.1:
				selected_first[spectra]=1
			if spectra==1000:
				print(numpy.var(data1))					
	print("Selected according to variance="+str(selected_first.count(1)))
	#Perform wilcoxon rank sam
	selected=[0]*l
	for spectra in range(l):
		if used_spectra[spectra]==0 or selected_first[spectra]==0:
			continue
		for i in range(len(classification_problems)):
			if selected[spectra]==1:
				break
			category1_mean=0.0
			category1_samples=0
			category2_mean=0.0
			category2_samples=0
			data1=list()
			data2=list()
			for j in range(len(quantitative_spectra_dataset)):
				if classification_problems[i][j]==1:
					if (float(quantitative_spectra_dataset_normalized[spectra][j])!=-1):						
						data1.append(quantitative_spectra_dataset[j][spectra])
				elif classification_problems[i][j]==-1:
					if (float(quantitative_spectra_dataset_normalized[spectra][j])!=-1):
						data2.append(quantitative_spectra_dataset[j][spectra])
			if len(data1)>1 and len(data2)>1:			
				[z,pvalue]=st.ranksums(data1,data2)
				if pvalue<individuals[0][3]:
					selected[spectra]=selected[spectra]+1
	if supervised_mode==1:
		for spectra in range(l):
			if selected[spectra]==1:
				selected[spectra]=1
			else:
				selected[spectra]=0
	else:
		for spectra in range(l):		
			if used_spectra[spectra]==0:
				selected[spectra]=0
			else:
				selected[spectra]=1
	print(len(selected))
	print(len(quantitative_spectra_dataset_normalized))
	normalized_fid = open(output_folder+"selected_normalized_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	for i in range(len(quantitative_spectra_dataset_normalized)-2):
		if selected[i]==1:
			normalized_fid.write(str(i)+"\t")
			for j in range(len(quantitative_spectra_dataset_normalized[0])):
				normalized_fid.write(str(quantitative_spectra_dataset_normalized[i][j])+"\t")
			normalized_fid.write("\n")
	normalized_fid.close()	
	missing_values_fid = open(output_folder+"selected_normalized_missing_values_completed_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	for i in range(len(quantitative_spectra_dataset_normalized_mv_completed)-2):
		if selected[i]==1:
			missing_values_fid.write(str(i)+"\t")
			for j in range(len(quantitative_spectra_dataset_normalized_mv_completed[0])):
				missing_values_fid.write(str(quantitative_spectra_dataset_normalized_mv_completed[i][j])+"\t")
			missing_values_fid.write("\n")
	missing_values_fid.close()
	print("Initial Feature Selection Completed Succesfully!")
	print("Selected Spectra:")
	print(str(selected.count(1)))
	print("Length of selected="+str(len(selected)))		
	number_of_selected_spectra=0
	number_of_selected_ptms=0
	refiltered=MSExperiment()
	metadata=list()
	for i in range(len(individuals[0])):
		if i>=4 and i<=15:
			if individuals[0][i]>1:
				number_of_selected_ptms=number_of_selected_ptms+1
		if i>=16:
			if selected[i-16]==1:
				number_of_selected_spectra=number_of_selected_spectra+1
				str_metadata=str(unified_spectra[i-16].getNativeID())+"_1"
				unified_spectra[i-16].setNativeID(str_metadata)
				str_metadata2="controllerType=0 controllerNumber=1 scan="+str(number_of_selected_spectra)
				metadata.append(str_metadata2)
				refiltered.addSpectrum(unified_spectra[i-16])
				refiltered[number_of_selected_spectra-1].setNativeID(str_metadata2)
	filename_out = output_folder+"differentially_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML"
	file = MzMLFile()		
	file.store(filename_out, refiltered)
	
	number_of_selected_spectra=0
	number_of_selected_ptms=0
	refiltered=MSExperiment()
	metadata=list()
	for i in range(len(individuals[0])):
		if i>=4 and i<=15:
			if individuals[0][i]>1:
				number_of_selected_ptms=number_of_selected_ptms+1
		if i>=16:
			if selected[i-16]==1 and individuals[0][i]>0.5:
				number_of_selected_spectra=number_of_selected_spectra+1
				unified_spectra[i-16].setNativeID(unified_spectra[i-16].getNativeID()+"_1")
				metadata.append(unified_spectra[i-16].getNativeID())
				refiltered.addSpectrum(unified_spectra[i-16])

	filename_out = output_folder+"differentially_quantified_spectra_list_used_in_classification_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML"
	file = MzMLFile()		
	file.store(filename_out, refiltered)	
	refiltered=MSExperiment()
	metadata=list()
	for i in range(len(individuals[0])):
		if i>=16:
			if used_spectra[i-16]==1:
				number_of_selected_spectra=number_of_selected_spectra+1
				refiltered.addSpectrum(unified_spectra[i-16])
	filename_out = output_folder+"commonly_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML"
	file = MzMLFile()		
	file.store(filename_out, refiltered)
	print("Correction of mzML files started!")	
	differentially_quantified_spectra_list_fid = open(output_folder+"differentially_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+".mzML","w")
	differentially_quantified_spectra_list_intermediate_fid = open(output_folder+"differentially_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML","r")
	num_of_lines=0
	
	for line1 in differentially_quantified_spectra_list_intermediate_fid:
		num_of_lines=num_of_lines+1
		differentially_quantified_spectra_list_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
	differentially_quantified_spectra_list_intermediate_fid.close()
	differentially_quantified_spectra_list_fid.close()
	differentially_quantified_spectra_list_fid = open(output_folder+"differentially_quantified_spectra_list_used_in_classification_"+str(time.strftime("%Y_%m_%d"))+".mzML","w")
	differentially_quantified_spectra_list_intermediate_fid = open(output_folder+"differentially_quantified_spectra_list_used_in_classification_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML","r")
	num_of_lines=0
	for line1 in differentially_quantified_spectra_list_intermediate_fid:
		num_of_lines=num_of_lines+1
		differentially_quantified_spectra_list_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
	differentially_quantified_spectra_list_intermediate_fid.close()
	differentially_quantified_spectra_list_fid.close()
	commonly_quantified_spectra_list_fid = open(output_folder+"commonly_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+".mzML","w")
	commonly_quantified_spectra_list_intermediate_fid = open(output_folder+"commonly_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML","r")
	num_of_lines=0	
	for line1 in commonly_quantified_spectra_list_intermediate_fid:
		num_of_lines=num_of_lines+1
		commonly_quantified_spectra_list_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
	commonly_quantified_spectra_list_intermediate_fid.close()
	commonly_quantified_spectra_list_fid.close()
	corrected_unified_spectra_list_fid = open(output_folder+"corrected_"+str(unified_list_filename),"w")
	unified_spectra_list_fid = open(step2_output_folder+unified_list_filename,"r")
	num_of_lines=0	
	for line1 in unified_spectra_list_fid:
		num_of_lines=num_of_lines+1
		corrected_unified_spectra_list_fid.write(line1.replace("controllerType=0 controllerNumber=1","controllerType=0 controllerNumber=1_"+str(num_of_lines))) 
	corrected_unified_spectra_list_fid.close()
	unified_spectra_list_fid.close()
	os.remove(output_folder+"differentially_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML")
	os.remove(output_folder+"differentially_quantified_spectra_list_used_in_classification_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML")	
	os.remove(output_folder+"commonly_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+"_intermediate.mzML")	
	print("Correction of mzML files finished!")	
	#Generate final solution report
	final_detailed_solutions_fid = open(output_folder+"final_detailed_solution_2nd_optimization_problem_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
	if individuals[0][0]<1:
		final_detailed_solutions_fid.write("Linear Kernel is selected.\n")	
	else:
		final_detailed_solutions_fid.write("Radial Basis Function kernel is selected.\n")
		final_detailed_solutions_fid.write("Parameter gamma of Radial Basis Function="+str(individuals[0][2])+"\n")
	final_detailed_solutions_fid.write("Regularization Parameter C of SVM="+str(individuals[0][1])+"\n")
	final_detailed_solutions_fid.write("Wilcoxon Rank Sum Feature Selection Threshold="+str(individuals[0][3])+"\n")
	final_detailed_solutions_fid.write("The following additional PTMs were selected to be used:\n")
	if individuals[0][4]>1:
		final_detailed_solutions_fid.write("Phosphorylation (ST)\n")
	if individuals[0][5]>1:
		final_detailed_solutions_fid.write("Pyro-glutamate (N-term Q)\n")
	if individuals[0][6]>1:
		final_detailed_solutions_fid.write("Carbamylation (K)\n")
	if individuals[0][7]>1:
		final_detailed_solutions_fid.write("N-terminal methionine\n")
	if individuals[0][8]>1:
		final_detailed_solutions_fid.write("N-terminal acetylation\n")
	if individuals[0][9]>1:
		final_detailed_solutions_fid.write("Acetylation (K) f\n")
	if individuals[0][10]>1:
		final_detailed_solutions_fid.write("Dihydroxy tryptophan (W)\n")
	if individuals[0][11]>1:
		final_detailed_solutions_fid.write("Methylation (K)\n")
	if individuals[0][12]>1:
		final_detailed_solutions_fid.write("Formylation (ST)\n")
	if individuals[0][13]>1:
		final_detailed_solutions_fid.write("Iron (ED)\n")
	if individuals[0][14]>1:
		final_detailed_solutions_fid.write("Iodoacetmaide (M)\n")
	if individuals[0][15]>1:
		final_detailed_solutions_fid.write("Iodination (Y)\n")
	#Write down differential quantified peptides report
	input_xml_fid = open("xtandem_input_step3/input.xml","w")
	input_xml_fid.write("<?xml version=\"1.0\"?>\n")
	input_xml_fid.write("<bioml>\n")
	input_xml_fid.write("\t<note>\n")
	input_xml_fid.write("\tEach one of the parameters for x! tandem is entered as a labeled note node.\n")
	input_xml_fid.write("\tAny of the entries in the default_input.xml file can be over-ridden by\n")
	input_xml_fid.write("\tadding a corresponding entry to this file. This file represents a minimum\n")
	input_xml_fid.write("\tinput file, with only entries for the default settings, the output file\n")
	input_xml_fid.write("\tand the input spectra file name.\n")
	input_xml_fid.write("\tSee the taxonomy.xml file for a description of how FASTA sequence list\n")
	input_xml_fid.write("\tfiles are linked to a taxon name.\n")
	input_xml_fid.write("\t</note>\n")
	message="\t<note>model refinement parameters</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine\">yes</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, modification mass\"></note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, sequence path\"></note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, tic percent\">20</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, spectrum synthesis\">yes</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, maximum valid expectation value\">0.1</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, potential N-terminus modifications\">+229.16@["
	if individuals[0][5]>=1:
		message=message+",+39.99@["
	if individuals[0][8]>=1:
		message=message+",+42.01@["
	message=message+"</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, potential C-terminus modifications\"></note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, unanticipated cleavage\">yes</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, potential modification mass\">+0.98@A, +0.98@Q, +15.99@M"
	if individuals[0][4]>=1:
		message=message+",+79.97@S,+79.97@T"
	if individuals[0][6]>=1:
		message=message+",+43.01@K"
	if individuals[0][9]>=1:
		message=message+",+42.01@K"
	if individuals[0][10]>=1:
		message=message+",+19.99@W"
	if individuals[0][11]>=1:
		message=message+",+14.02@K"
	if individuals[0][12]>=1:
		message=message+",+27.79@S,+27.79@T"
	if individuals[0][13]>=1:
		message=message+",+53.92@D,+53.92@E"
	if individuals[0][15]>=1:
		message=message+",+125.90@Y"
	message=message+"</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, point mutations\">no</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, use potential modifications for full refinement\">yes</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, point mutations\">no</note>\n"
	message=message+"\t<note type=\"input\" label=\"refine, potential modification motif\"></note>\n"
	message=message+"\t<note>The format of this parameter is similar to residue, modification mass, with the addition of a modified PROSITE notation sequence motif specification. For example, a value of 80@[ST!]PX[KR] indicates a modification of either S or T when followed by P, and residue and the a K or an R. A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it is NOT followed by a P, then either an S or a T, NOT followed by a P. Positive and negative values are allowed. </note>\n"
	input_xml_fid.write(message)
	input_xml_fid.write("\t<note type=\"input\" label=\"list path, default parameters\">default_input.xml</note>\n")
	input_xml_fid.write("\n")
	input_xml_fid.write("\t<note type=\"input\" label=\"list path, taxonomy information\">taxonomy.xml</note>\n")
	input_xml_fid.write("\n")
	input_xml_fid.write("\t<note type=\"input\" label=\"protein, taxon\">human</note>\n")
	input_xml_fid.write("\n")
	input_xml_fid.write("\t<note type=\"input\" label=\"protein, taxon\">Escherichia coli</note>\n")
	input_xml_fid.write("\n")
	input_xml_fid.write("\t<note type=\"input\" label=\"spectrum, path\">"+output_folder+"differentially_quantified_spectra_list_"+str(time.strftime("%Y_%m_%d"))+".mzML</note>\n")
	input_xml_fid.write("\n")
	input_xml_fid.write("\t<note type=\"input\" label=\"output, path\">"+output_folder+"quantified_peptides_proteins_"+str(time.strftime("%Y_%m_%d"))+".xml</note>\n")
	input_xml_fid.write("</bioml>\n")
	input_xml_fid.close()
	command = ["./tandem.exe","xtandem_input_step3/input.xml"]
	result_p=Popen(command)
	print("Result xTandem="+str(result_p))	
	#reading xTandem results
	result_p.communicate()	
	#write down SVM classification models
	print("Classification models training and testing started!")
	inputs=list()
	outputs=list()
	num_of_samples=0
	for sample in range(len(quantitative_spectra_dataset)):
		inputs.append([])
		for feature in range(len(individuals[0])-16):
			if selected[feature]==1 and individuals[0][feature+16]>0.5:				
				inputs[num_of_samples].append(quantitative_spectra_dataset[sample][feature])
		outputs.append(int(quantitative_spectra_dataset[sample][-1]))
		num_of_samples=num_of_samples+1
	print("Datasets Created Successfully")
	fold_size=math.floor(float(len(quantitative_spectra_dataset))/num_of_folds)
	accuracy=0
	if not os.path.exists(output_folder+"classification_models"):
    		os.makedirs(output_folder+"classification_models/")
	for k in range(num_of_folds):
		training_inputs=list()
		training_outputs=list()
		testing_inputs=list()
		testing_outputs=list()
		testing_samples=0
		training_samples=0
		for j in range(len(quantitative_spectra_dataset)):
			if j>=k*fold_size and j<(k+1)*fold_size:			
				testing_inputs.append([])
				for z in range(len(inputs[j])):
					testing_inputs[testing_samples].append(inputs[j][z])
				testing_outputs.append(outputs[j])
				testing_samples=testing_samples+1
			else:
				training_inputs.append([])
				for z in range(len(inputs[j])):
					training_inputs[training_samples].append(inputs[j][z])
				training_outputs.append(outputs[j])
				training_samples=training_samples+1	
		prob = svm_problem(training_outputs, training_inputs)
		if individuals[0][0]<1:
			#linear svm
			param=svm_parameter("-t 0 -q -c "+str(individuals[0][1]))
		else:
			#RBF svm
			param=svm_parameter("-t 2 -q -c "+str(individuals[0][1])+" -g "+str(individuals[0][2]))
		model = svm_train(prob, param)
		p_labs, p_acc, p_vals = svm_predict(testing_outputs, testing_inputs, model)
		accuracy=accuracy+p_acc[0]
		svm_save_model(output_folder+"classification_models/model_"+str(k)+"_"+str(time.strftime("%Y_%m_%d")), model)
	accuracy=accuracy/(float(num_of_folds)*100)
	#parse filtered files and create differentially_quantified mzml files
	if not os.path.exists(output_folder+"differentially_quantified_files"):
    		os.makedirs(output_folder+"differentially_quantified_files/")
	print("Differentially quantified spectra individual experiment mzML files are being created!")
	for i in range(len(mapping_spectra[0])):
		print(i)
		filtered_file = MSExperiment()
		file = MzMLFile()
		file.load(step2_output_folder+"filtered_files/"+str(mapping_spectra[0][i].strip()), filtered_file)
		refiltered=MSExperiment()
		number_of_selected_spectra=0		
		for j in range(filtered_file.size()):
			if int(mapping_spectra[j+1][i])!=-1:
				if int(mapping_spectra[j+1][i])<len(selected):				
					if selected[int(mapping_spectra[j+1][i])]==1:
						number_of_selected_spectra=number_of_selected_spectra+1
						refiltered.addSpectrum(filtered_file[j])
		file = MzMLFile()
		file.store(output_folder+"differentially_quantified_files/differentially_quantified_"+str(mapping_spectra[0][i].strip()),refiltered)
	print("successfully finished!")

def initialize(population,semi_random_initialization_probability,min_values, max_values):
	#population:(integer) (default: 10) it is the number of individual solutions which are evolved on parallel
	#semi_random_initialization_probability: (float) (default: 0.8) it is the proportion (multiplied with 100) of the individuals which will be initialized with the semi automatic method. The rest will be initialized using random method
	#min_values: (list) it is a list with the minimum allowed values for the variables which should be optimized
	#max_values: (list) it is a list with the maximum allowed values for the variables which should be optimized
	
	individuals=[[0 for x in range(len(min_values))] for x in range(population)] 
			
	for i in range(len(individuals)):
		if random.uniform(0,1)<semi_random_initialization_probability:
			individuals[i][0]=random.uniform(min_values[0],max_values[0])
			individuals[i][1]=random.gauss(1.0,0.1*(max_values[1]-min_values[1]))	
			individuals[i][2]=random.gauss(1.0,0.1*(max_values[2]-min_values[2]))		
			individuals[i][3]=random.gauss(0.25,0.1*(max_values[3]-min_values[3]))		
			individuals[i][4]=random.gauss(1.5,0.1*(max_values[4]-min_values[4]))	
			individuals[i][5]=random.gauss(1.5,0.1*(max_values[5]-min_values[5]))
			individuals[i][6]=random.gauss(1.5,0.1*(max_values[6]-min_values[6]))	
			individuals[i][7]=random.gauss(1.5,0.1*(max_values[7]-min_values[7]))	
			individuals[i][8]=random.gauss(1.5,0.1*(max_values[8]-min_values[8]))	
			individuals[i][9]=random.gauss(1.5,0.1*(max_values[9]-min_values[9]))	
			individuals[i][10]=random.gauss(1.5,0.1*(max_values[10]-min_values[10]))	
			individuals[i][11]=random.gauss(1.5,0.1*(max_values[11]-min_values[11]))	
			individuals[i][12]=random.gauss(0.5,0.1*(max_values[12]-min_values[12]))
			individuals[i][13]=random.gauss(0.5,0.1*(max_values[13]-min_values[13]))
			individuals[i][14]=random.gauss(0.5,0.1*(max_values[14]-min_values[14]))
			individuals[i][15]=random.gauss(0.5,0.1*(max_values[15]-min_values[15]))
			#individuals[i][16]=random.gauss(0.5,0.1*(max_values[16]-min_values[16]))	
			for j in range(16,len(min_values)):
				r=random.uniform(0.0,1.0)
				if r<=0.8:
					individuals[i][j]=random.uniform(0.5,1.0)
				else:
					individuals[i][j]=random.uniform(0.0,0.5)
			for j in range(0, len(min_values)):
				if individuals[i][j]<min_values[j]:
					individuals[i][j]=min_values[j]
				if individuals[i][j]>max_values[j]:
					individuals[i][j]=max_values[j]
					
		else:
			for j in range(0, len(min_values)):
				individuals[i][j]=random.uniform(min_values[j],max_values[j])
	return individuals
#a function that counts the elements of a list within an interval
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

def dominate(solution1, solution2):
	check=2
	ffs=len(solution1)
	dominate1=1
	equal1=1
	f=0
	while f<ffs and dominate1==1:
		if solution1[f]>solution2[f]:
			equal1=0
		elif solution1[f]==solution2[f]:
			do_nothing=1
		else:
			dominate1=0
		f=f+1
	if dominate1==1 and equal1==0:
		check=1
	elif dominate1==1 and equal1==1:
		check=2
	else:
		dominate2=1
		equal2=1
		f=0
		while f<ffs and dominate2==1:
			if solution2[f]>solution1[f]:
				equal2=0
			elif solution2[f]==solution1[f]:
				do_nothing=1
			else:
				dominate2=0
			f=f+1
		if dominate2==1 and equal2==0:
			check=3
	return check

def evaluate_individuals_test(dataset,dataset_with_missing_values,unified_spectra,individuals,goal_significances,database_name,spectraST_databases_filename,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder):
	evaluation_values=[[] for i in range(8)]
	
	for i in range(len(individuals)):
		r1=random.uniform(0,1)
		r2=random.uniform(0,1)
		r3=random.uniform(0,1)
		r4=random.uniform(0,1)
		r5=random.uniform(0,1)
		r6=random.uniform(0,1)
		r7=random.uniform(0,1)
		
		evaluation_values[0].append(r1)
		evaluation_values[1].append(r2)
		evaluation_values[2].append(r3)
		evaluation_values[3].append(r4)
		evaluation_values[4].append(r5)
		evaluation_values[5].append(r6)
		evaluation_values[6].append(r7)
		evaluation_values[7].append(r1+r2+r3+r4+r5+r6+r7)
	return evaluation_values
def evaluate_individuals(dataset,dataset_with_missing_values,unified_spectra,individuals,goal_significances,database_name,spectraST_databases_filename,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder):
	#dataset:(List) a two dimensional list including the quantitative values and the classes information and formulated as a dataset. Each line corresponds to one sample, the first n columns are the quantitative values of all spectra in the unified spectra list and last line is the classes information. This dataset has been normalized and missing values have been imputed
	#dataset_with_missing_values:(List) a two dimensional list including the quantitative values and the classes information and formulated as a dataset. Each line corresponds to one sample, the first n columns are the quantitative values of all spectra in the unified spectra list and last line is the classes information. This dataset has been normalized but missing values have not been imputed
	#individuals: (two dimensional list) including the population of individual solutions
	#goal_significances: (list of floats) Includs the significances (values from 0 to 1) of the individual goals
	#database_name: (String) the string of the name of the fasta database to be used for scaffold
	#SpectraST_databases_filename: (String) the filename of the file which includes the names of the libraries which will be used by spectraST software
	#num_of_folds: (integer) the number of k folds in k-fold cross validation
	#used_spectra: (list) list with info for selected or not selected spectra
	#classification_problems: (list) list of different multiclassification problems definitions
	#number_of_used_tmt_channels: (int) number of used tmt channels as described in the relevant input file
	#supervised_mode: (int) 1 for supervised_learning and 0 for unsupervised learning
	#output_folder: (string) step3 folder
	evaluation_values=[[] for j in range(8)]
	for ind in range(len(individuals)):
		
		
		#evaluation_values=list()
		l=len(individuals[ind])-16
		selected_first=[0]*l
		for spectra in range(l):
			if used_spectra[spectra]==0:
				continue
			

			data1=list()
			for j in range(len(dataset)):
				if (float(dataset_with_missing_values[spectra][j])!=-1):						
					data1.append(dataset[j][spectra])
			if len(data1)>1:
				if numpy.var(data1)>0.1:
					selected_first[spectra]=1

		#Perform wilcoxon rank sam
		selected=[0]*l
		for spectra in range(l):
			if used_spectra[spectra]==0 or selected_first[spectra]==0:
				continue
			for i in range(len(classification_problems)):
				if selected[spectra]==1:
					break
				category1_mean=0.0
				category1_samples=0
				category2_mean=0.0
				category2_samples=0
				data1=list()
				data2=list()
				for j in range(len(dataset)):
					if classification_problems[i][j]==1:
						if (float(dataset_with_missing_values[spectra][j])!=-1):						
							data1.append(dataset[j][spectra])
					elif classification_problems[i][j]==-1:
						if (float(dataset_with_missing_values[spectra][j])!=-1):
							data2.append(dataset[j][spectra])
				if len(data1)>1 and len(data2)>1:			
					[z,pvalue]=st.ranksums(data1,data2)
					if pvalue<individuals[ind][3]:
						selected[spectra]=selected[spectra]+1
		if supervised_mode==1:		
			for spectra in range(l):
				if selected[spectra]==1:
					selected[spectra]=1
				else:
					selected[spectra]=0
		else:
			for spectra in range(l):
				if used_spectra[spectra]==0:
					selected[spectra]=0
				else:
					selected[spectra]=1
		print("Initial Feature Selection Completed Succesfully!")
		refiltered=MSExperiment()
		number_of_selected_spectra=0
		number_of_selected_ptms=0
		for i in range(len(individuals[ind])):
			if i>=4 and i<=15:
				if individuals[ind][i]>1:
					number_of_selected_ptms=number_of_selected_ptms+1
			if i>=16:
				if supervised_mode==1:
					if selected[i-16]==1 and individuals[ind][i]>0.5:
						number_of_selected_spectra=number_of_selected_spectra+1
						refiltered.addSpectrum(unified_spectra[i-16])
				else:
					if selected[i-16]==1:
						number_of_selected_spectra=number_of_selected_spectra+1
						refiltered.addSpectrum(unified_spectra[i-16])
		if number_of_selected_spectra==0:
			goal1=0
		goal1=1/(0.01+number_of_selected_spectra)
		goal2=1/(0.01+number_of_selected_ptms)
		goal3=selected.count(1)/float(l)
		#xTandem Search
		#create or empty input/output folders
		if not os.path.exists(output_folder+'xtandem_input_step3_'+str(ind)):
				os.makedirs(output_folder+'xtandem_input_step3_'+str(ind))
		if not os.path.exists(output_folder+'xtandem_output_step3_'+str(ind)):
				os.makedirs(output_folder+'xtandem_output_step3_'+str(ind))
		xtandem_input_folder = output_folder+'xtandem_input_step3_'+str(ind)
		for the_file in os.listdir(xtandem_input_folder):
			file_path = os.path.join(xtandem_input_folder, the_file)
			try:
				if os.path.isfile(file_path):
					os.unlink(file_path)
			except Exception, e:
				print e
		xtandem_output_folder = output_folder+'xtandem_output_step3_'+str(ind)
		for the_file in os.listdir(xtandem_output_folder):
			file_path = os.path.join(xtandem_output_folder, the_file)
			try:
				if os.path.isfile(file_path):
					os.unlink(file_path)
			except Exception, e:
				print e
		filename_out = output_folder+'differentially_quantified_spectra_list_'+str(ind)+'.mzML'
		file = MzMLFile()		
		file.store(filename_out, refiltered)
		input_xml_fid = open(output_folder+"xtandem_input_step3_"+str(ind)+"/input.xml","w")
		input_xml_fid.write("<?xml version=\"1.0\"?>\n")
		input_xml_fid.write("<bioml>\n")
		input_xml_fid.write("\t<note>\n")
		input_xml_fid.write("\tEach one of the parameters for x! tandem is entered as a labeled note node.\n")
		input_xml_fid.write("\tAny of the entries in the default_input.xml file can be over-ridden by\n")
		input_xml_fid.write("\tadding a corresponding entry to this file. This file represents a minimum\n")
		input_xml_fid.write("\tinput file, with only entries for the default settings, the output file\n")
		input_xml_fid.write("\tand the input spectra file name.\n")
		input_xml_fid.write("\tSee the taxonomy.xml file for a description of how FASTA sequence list\n")
		input_xml_fid.write("\tfiles are linked to a taxon name.\n")
		input_xml_fid.write("\t</note>\n")
		message="\t<note>model refinement parameters</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine\">yes</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, modification mass\"></note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, sequence path\"></note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, tic percent\">20</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, spectrum synthesis\">yes</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, maximum valid expectation value\">0.1</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, potential N-terminus modifications\">+229.16@["
		if individuals[ind][5]>=1:
			message=message+",+39.99@["
		if individuals[ind][8]>=1:
			message=message+",+42.01@["
		message=message+"</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, potential C-terminus modifications\"></note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, unanticipated cleavage\">yes</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, potential modification mass\">+0.98@A, +0.98@Q, +15.99@M"
		if individuals[ind][4]>=1:
			message=message+",+79.97@S,+79.97@T"
		if individuals[ind][6]>=1:
			message=message+",+43.01@K"
		if individuals[ind][9]>=1:
			message=message+",+42.01@K"
		if individuals[ind][10]>=1:
			message=message+",+19.99@W"
		if individuals[ind][11]>=1:
			message=message+",+14.02@K"
		if individuals[ind][12]>=1:
			message=message+",+27.79@S,+27.79@T"
		if individuals[ind][13]>=1:
			message=message+",+53.92@D,+53.92@E"
		if individuals[ind][15]>=1:
			message=message+",+125.90@Y"
		message=message+"</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, point mutations\">no</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, use potential modifications for full refinement\">yes</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, point mutations\">no</note>\n"
		message=message+"\t<note type=\"input\" label=\"refine, potential modification motif\"></note>\n"
		message=message+"\t<note>The format of this parameter is similar to residue, modification mass, with the addition of a modified PROSITE notation sequence motif specification. For example, a value of 80@[ST!]PX[KR] indicates a modification of either S or T when followed by P, and residue and the a K or an R. A value of 204@N!{P}[ST]{P} indicates a modification of N by 204, if it is NOT followed by a P, then either an S or a T, NOT followed by a P. Positive and negative values are allowed. </note>\n"
		input_xml_fid.write(message)
		input_xml_fid.write("\t<note type=\"input\" label=\"list path, default parameters\">default_input.xml</note>\n")
		input_xml_fid.write("\n")
		input_xml_fid.write("\t<note type=\"input\" label=\"list path, taxonomy information\">taxonomy.xml</note>\n")
		input_xml_fid.write("\n")
		input_xml_fid.write("\t<note type=\"input\" label=\"protein, taxon\">human</note>\n")
		input_xml_fid.write("\n")
		input_xml_fid.write("\t<note type=\"input\" label=\"protein, taxon\">Escherichia coli</note>\n")
		input_xml_fid.write("\n")
		input_xml_fid.write("\t<note type=\"input\" label=\"spectrum, path\">"+output_folder+"differentially_quantified_spectra_list_"+str(ind)+".mzML</note>\n")
		input_xml_fid.write("\n")
		input_xml_fid.write("\t<note type=\"input\" label=\"output, path\">"+output_folder+"xtandem_output_step3_"+str(ind)+"/output.xml</note>\n")
		input_xml_fid.write("</bioml>\n")
		input_xml_fid.close()
		command = ["./tandem.exe",output_folder+"xtandem_input_step3_"+str(ind)+"/input.xml"]
		result_p=Popen(command)
		#reading xTandem results
		result_p.communicate()			
		xtandem_output_folder = output_folder+'xtandem_output_step3_'+str(ind)
		identified_xtandem=1
		identified_xtandem2=0
		proteins=list()
		for the_file in os.listdir(xtandem_output_folder):
			file_path = os.path.join(xtandem_output_folder, the_file)
			xtandem_result_file=file_path
			output_xml_fid = open(file_path,"r")
			for line in output_xml_fid:
				if line.find("<protein expect=")!=-1:
					if line.find("label=\"")!=-1:
						words=line.split("label=\"")
						words2=words[1].split("\"")
						try:
							index_in_proteins=proteins.index(words2[0])
						except ValueError:
							index_in_proteins=-1
						if index_in_proteins==-1:
							proteins.append(words2[0])
				if (line.find("total spectra assigned")!=-1):
					words=line.split(">")
					words2=words[1].split("<")
					identified_xtandem=int(words2[0])
				if (line.find("total unique assigned")!=-1):
					words=line.split(">")
					words2=words[1].split("<")
					identified_xtandem2=int(words2[0])
			output_xml_fid.close()
		if identified_xtandem>0:
			goal5=	identified_xtandem2/float(identified_xtandem)	
		else:
			goal5=0
		if goal5>1:
			goal5=1	
		if float(len(selected))>0.0:
			goal4=identified_xtandem2/float(len(selected))
		else:
			goal4=0
		goal6=len(proteins)/float(500)
		evaluation_values[0].append(goal1)
		evaluation_values[1].append(goal2)
		evaluation_values[2].append(goal3)
		evaluation_values[3].append(goal4)
		evaluation_values[4].append(goal5)
		evaluation_values[5].append(goal6)
	obj_individuals = []
	#create a list of Individuals objects
	for i, individual in enumerate(individuals):
		obj_individuals.append(Individual(individual, i, dataset,dataset_with_missing_values,unified_spectra,goal_significances,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder))
	pool = mp.Pool(processes=mp.cpu_count()-1)
	
	results = pool.map(evaluate_individual, obj_individuals, chunksize=1)
	
	pool.close()
	pool.join()
	for i in range(len(results)):
		evaluation_values[6].append(float(results[i][0]))	
	for i in range(len(evaluation_values[0])):
		evaluation_values[7].append((float(evaluation_values[0][i])*float(goal_significances[0])+float(evaluation_values[1][i])*float(goal_significances[1])+float(evaluation_values[2][i])*float(goal_significances[2])+float(evaluation_values[3][i])*float(goal_significances[3])+float(evaluation_values[5][i])*float(goal_significances[5])+float(evaluation_values[6][i])*float(goal_significances[6]))/6.0)
	for i in range(len(evaluation_values[0])):
		print '-----------------Individual:'+str(i)+'----------------------------'
		print '--------------------------------------------------------------'
		print '|     |  1  |  Minimizing selected by the classifier spetra | %s ' % evaluation_values[0][i]
		print '|     |  2  |  Minimizing PTMs    			   | %s ' % evaluation_values[1][i]
		print '|  G  |  3  |  Differentially Quantified Spectra            | %s ' % evaluation_values[2][i]
		print '|  O  |  4  |  Xtandem Identidied Peptides      		   | %s ' % evaluation_values[3][i]
		print '|  A  |  5  |  Percentage of Validated Peptides Matches	   | %s ' % evaluation_values[4][i]
		print '|  L  |  6  |  Quantified Proteins			   | %s ' % evaluation_values[5][i]
		print '|  S  |  7  |  Classifiers Accuracy			   | %s ' % evaluation_values[6][i]
		print '|     |  8  |  Weighted Sum of Goals 			   | %s ' % evaluation_values[7][i]
		print '--------------------------------------------------------------'	
	return evaluation_values
def evaluate_individual(individual):
	return individual.evaluate()	
class Individual():
	def __init__(self, individual, process_i, dataset,dataset_with_missing_values,unified_spectra,goal_significances,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder):
		self.individual = individual
		self.dataset= dataset
		self.process_i=process_i
		self.dataset_with_missing_values = dataset_with_missing_values
		self.unified_spectra = unified_spectra
		self.goal_significances = goal_significances
		self.num_of_folds = num_of_folds
		self.used_spectra = used_spectra
		self.classification_problems = classification_problems
		self.number_of_used_tmt_channels = number_of_used_tmt_channels
		self.supervised_mode = supervised_mode
		self.output_folder=output_folder
	
	def evaluate(self):
		print '\n\n> Training models for individual ## %s ################################' % (self.process_i)
		evaluation_values=list()
		classification_performances=list()
		l=len(self.individual)-16
		selected_first=[0]*l
		for spectra in range(l):
			if self.used_spectra[spectra]==0:
				continue
			

			data1=list()
			data2=list()
			for j in range(len(self.dataset)):
				if (float(self.dataset_with_missing_values[spectra][j])!=-1):						
					data1.append(self.dataset[j][spectra])
			if len(data1)>1:
				if numpy.var(data1)>0.1:
					selected_first[spectra]=1

		#Perform wilcoxon rank sam
		selected=[0]*l
		for spectra in range(l):
			if self.used_spectra[spectra]==0 or selected_first[spectra]==0:
				continue
			for i in range(len(self.classification_problems)):
				if selected[spectra]==1:
					break
				category1_mean=0.0
				category1_samples=0
				category2_mean=0.0
				category2_samples=0
				data1=list()
				data2=list()
				for j in range(len(self.dataset)):
					if self.classification_problems[i][j]==1:
						if (float(self.dataset_with_missing_values[spectra][j])!=-1):						
							data1.append(self.dataset[j][spectra])
					elif self.classification_problems[i][j]==-1:
						if (float(self.dataset_with_missing_values[spectra][j])!=-1):
							data2.append(self.dataset[j][spectra])
				if len(data1)>1 and len(data2)>1:			
					[z,pvalue]=st.ranksums(data1,data2)
					if pvalue<self.individual[3]:
						selected[spectra]=selected[spectra]+1
		if self.supervised_mode==1:		
			for spectra in range(l):
				if selected[spectra]==1:
					selected[spectra]=1
				else:
					selected[spectra]=0
		else:
			for spectra in range(l):
				if self.used_spectra[spectra]==0:
					selected[spectra]=0
				else:
					selected[spectra]=1
		print("Initial Feature Selection Completed Succesfully!")
		print("Selected Spectra:")
		print(str(selected.count(1)))
		number_of_selected_spectra=0
		number_of_selected_ptms=0
		for i in range(len(self.individual)):
			if i>=4 and i<=15:
				if self.individual[i]>1:
					number_of_selected_ptms=number_of_selected_ptms+1
			if i>=16:
				if self.supervised_mode==1:
					if selected[i-16]==1 and self.individual[i]>0.5:
						number_of_selected_spectra=number_of_selected_spectra+1
						
				else:
					if selected[i-16]==1:
						number_of_selected_spectra=number_of_selected_spectra+1
		
		#Classification models testing
		print("Classification models training and testing started!")
		if float(len(selected))>0.0:
			inputs=list()
			outputs=list()
			num_of_samples=0
			for sample in range(len(self.dataset)):
				inputs.append([])
				for feature in range(len(self.individual)-16):
					if selected[feature]==1 and self.individual[feature+16]>0.5:				
						inputs[num_of_samples].append(self.dataset[sample][feature])
				outputs.append(int(self.dataset[sample][-1]))
				num_of_samples=num_of_samples+1
			print("Datasets Created Successfully")
			fold_size=math.floor(float(len(self.dataset))/self.num_of_folds)
			accuracy=0
			for k in range(self.num_of_folds):
				training_inputs=list()
				training_outputs=list()
				testing_inputs=list()
				testing_outputs=list()
				testing_samples=0
				training_samples=0
				for j in range(len(self.dataset)):
					if j>=k*fold_size and j<(k+1)*fold_size:
						testing_inputs.append([])
						for z in range(len(inputs[j])):
							testing_inputs[testing_samples].append(inputs[j][z])
						testing_outputs.append(outputs[j])
						testing_samples=testing_samples+1
					else:
						training_inputs.append([])
						for z in range(len(inputs[j])):
							training_inputs[training_samples].append(inputs[j][z])
						training_outputs.append(outputs[j])
						training_samples=training_samples+1	
				prob = svm_problem(training_outputs, training_inputs)
				if self.individual[0]<1:
					#linear svm
					param=svm_parameter("-t 0 -q -c "+str(self.individual[1]))
				else:
					#RBF svm
					param=svm_parameter("-t 2 -q -c "+str(self.individual[1])+" -g "+str(self.individual[2]))
				model = svm_train(prob, param)
				p_labs, p_acc, p_vals = svm_predict(testing_outputs, testing_inputs, model)

				accuracy=accuracy+p_acc[0]
			accuracy=accuracy/(float(self.num_of_folds)*100)

			if self.supervised_mode==1:	
				goal6=accuracy	
			else:
				goal6=1
		else:
			goal6=0
		
		
		evaluation_values.append(goal6)
				
		return evaluation_values
if __name__ == "__main__":
	#Initialize random generator seed with the current local time in miliseconds
	random.seed()
	min_values=[0.0,0.001,0.001,0.001,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0]
	max_values=[2.0,1000.0,1000.0,0.3,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0,2.0]
	

	unified_list_filename=sys.argv[1]
	mapping_filename=sys.argv[2]	
	quantitative_filename=sys.argv[3]
	labels_filename=sys.argv[4]
	population=int(sys.argv[5])
	generations=int(sys.argv[6])
	semi_random_initialization_probability=float(sys.argv[7])
	two_points_crossover_probability=float(sys.argv[8])
	arithmetic_crossover_probability=float(sys.argv[9])
	mutation_probability=float(sys.argv[10])
	gaussian_mutation_variance_proportion=float(sys.argv[11])
	goal_significances_filename=sys.argv[12]
	database_name=sys.argv[13]
	num_of_folds=int(sys.argv[14])
	step2_output_folder=sys.argv[15]
	replicates=int(sys.argv[16])
	missing_values_filename=sys.argv[17]
	use_inputs_vector_filename=sys.argv[18]
	supervised_mode=int(sys.argv[19])
	goal_significances_fid = open(goal_significances_filename,"r")
	tstamp = time.strftime('%Y_%m_%d')
	output_folder='step3/'+str(tstamp)+'_'+str(time.time())+'/'
	if not os.path.exists(output_folder):
		os.makedirs(output_folder)
	number_of_lines=0
	for line in goal_significances_fid:
		number_of_lines=number_of_lines+1
	goal_significances_fid.close()
	goal_significances=list()
	goal_significances_fid = open(goal_significances_filename,"r")
	quantitative_fid=open(step2_output_folder+quantitative_filename,"r")
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
	
	step3_supervised_analysis(unified_list_filename,mapping_filename,quantitative_filename,labels_filename,population,generations,semi_random_initialization_probability,two_points_crossover_probability,arithmetic_crossover_probability,mutation_probability,gaussian_mutation_variance_proportion,min_values,max_values,goal_significances, database_name,num_of_folds,step2_output_folder,output_folder,replicates, missing_values_filename,use_inputs_vector_filename,supervised_mode)
