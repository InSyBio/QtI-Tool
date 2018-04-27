####################################
#This code is the outcome of the project Innovative Processing of TMT Proteomics. This is a common project of InSyBio LTD and Nestle Institute of Health Sciences. 
#This program implements the step 2 code
####################################
from pyopenms import *
import os, shutil
import random
import math
import pyopenms
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

#Step 1 functions

class Utils():

	@staticmethod
	def get_file_line_count(fname):
		with open(fname) as f:
			for i, l in enumerate(f):
				pass
		return i + 1

	@staticmethod
	def remove_files_in_dir(path):
		filelist = [ os.path.join(path, f) for f in os.listdir(path) if os.path.isfile(os.path.join(path, f)) ]
		for f in filelist:
			try:
				os.remove(f)
			except Exception, e:
				print e

	@staticmethod
	def generate_random_list(uniform, min_values, max_values):
		# generate non-uniform random values
		if not uniform:
			random_list = [ random.uniform(min_values[0],  max_values[0]),
							random.uniform(min_values[1],  max_values[1]),
							random.gauss  (0.2,			0.1),
							random.gauss  (10.0,		   0.1*(max_values[3]-min_values[3])),
							random.gauss  (11.0,		   0.1*(max_values[4]-min_values[4])),
							random.gauss  (4,			  0.1*(max_values[5]-min_values[5])),
							random.uniform(min_values[6],  max_values[6]),
							random.uniform(min_values[7],  max_values[7]),
							random.gauss  (2,			  0.1*(max_values[5]-min_values[5])),
							random.uniform(min_values[9],  max_values[9]),
							random.uniform(min_values[10], max_values[10]),
							random.gauss  (0.1,			0.1),
							random.uniform(min_values[12], max_values[12]),
							random.gauss  (0.05,		   0.05),
							random.gauss  (200,			0.1*(max_values[14]-min_values[14])),
							random.uniform(min_values[15], max_values[15]),
							random.uniform(min_values[16], max_values[16]),
							random.gauss  (1,			  0.1*(max_values[17]-min_values[17])),
							random.gauss  (200,			0.1*(max_values[18]-min_values[18])),
							random.gauss  (30,			 0.1*(max_values[19]-min_values[19])),
							random.gauss  (10,			 0.1*(max_values[20]-min_values[20])),
							random.gauss  (1,			  0.1*(max_values[21]-min_values[21])),
							random.gauss  (0.15,		   0.15),
							random.gauss  (0.25,		   0.25) ]

			# limit the random vals between the min and max
			random_list = [min(max(min_values[i], val), max_values[i]) for i, val in enumerate(random_list)]

		# generate uniform random values
		else:
			random_list = [random.uniform(min_values[i], max_values[i]) for i, v in enumerate(min_values)]

		return random_list


class Filters():

	@staticmethod
	def _apply(filter, data):
		for d in data:
			if(d.getMSLevel() != 1):
				filter.filterSpectrum(d)

	@staticmethod
	def _apply_peak(filter, data):
		for d in data:
			if(d.getMSLevel() != 1):
				filter.filterPeakSpectrum(d)

	@classmethod
	def apply(cls, methods, *args, **kwargs):
		print '\n>> Applying filters: %s' %(' -> '.join(methods))
		for m in methods:
			getattr(cls, m)(*args, **kwargs)


	@staticmethod
	def denoise(individual, parsed_data=None, parsed_data_picked=None):
		if individual[1] < 1:
			return

		if individual[1] < 2:
			#Gaussian Normalization
			filter = GaussFilter()
			param  = pyopenms.Param()
			param.setValue('gaussian_width', individual[2], 'Gaussian width')
			param.setValue('ppm_tolerance',  individual[3], ('Gaussian width, depending on the m/z position. '
															 'The higher the value, the wider the peak and therefore '
															 'the wider the gaussian.'))
		else:
			#Savitzky Golay Normalization
			filter = SavitzkyGolayFilter()
			param  = pyopenms.Param()
			param.setValue('frame_length', int(round(individual[4])), ('The number of subsequent data points used for smoothing. '
																	   'This number has to be uneven. '
																	   'If it is not, 1 will be added.'))
			if (int(round(individual[5])) >= int(round(individual[4]))):
				individual[5]=int(round(individual[4]))-1

			param.setValue('polynomial_order', int(round(individual[5])), 'Order or the polynomial that is fitted.')

		filter.setParameters(param)
		filter.filterExperiment(parsed_data)

	@staticmethod
	def precursor_removal(individual, parsed_data=None, parsed_data_picked=None):
		if individual[7] < 1:
			return

		#ParentPeakMower Method
		filter = ParentPeakMower()
		param  = pyopenms.Param()
		param.setValue('window_size',float(individual[8]), ('The size of the m/z window '
															'where the peaks are removed, +/- window_size.'))
		filter.setParameters(param)

		if individual[16] < 1:
			Filters._apply(filter, parsed_data)
		elif parsed_data_picked:
			Filters._apply_peak(filter, parsed_data_picked)
		else:
			Filters._apply(filter, parsed_data)

	@staticmethod
	def normalizer(individual, parsed_data=None, parsed_data_picked=None):
		if individual[10] < 1:
			return

		if individual[10] < 2:
			filter = Normalizer()
		elif individual[10] < 3:
			filter = Scaler()
		elif individual[10] < 4:
			filter = SqrtMower()
		else:
			filter = BernNorm()
			param  = pyopenms.Param()
			param.setValue('threshold', individual[11], 'Threshold of the Bern et al. normalization')
			filter.setParameters(param)

		if individual [16] < 1:
			Filters._apply(filter, parsed_data)
		elif parsed_data_picked:
			Filters._apply_peak(filter, parsed_data_picked)
		else:
			Filters._apply(filter, parsed_data)

	@staticmethod
	def peak_picking(individual, parsed_data=None, parsed_data_picked=None):
		if individual[16] < 1:
			return

		if individual[16] < 2:
			pp	= PeakPickerHiRes()
			param = pyopenms.Param()
			param.setValue('signal_to_noise', individual[17], ('Minimal signal-to-noise ratio '
															   'for a peak to be picked (0.0 disables SNT '
															   'estimation!)'))

			#param.setValue('ms_levels',[2],'List of MS levels for which the peak picking is applied. Other scans are copied to the output without changes.')
			#param.setValue('SignalToNoise:win_len',individual[18],'window length in Thomson')
			#param.setValue('SignalToNoise:bin_count',int(round(individual[19])),'number of bins for intensity values')
			#param.setValue('SignalToNoise:min_required_elements',int(round(individual[20])),'minimum number of elements required in a window (otherwise it is considered sparse)')
			#param.setValue('SignalToNoise:write_log_messages','false','Write out log messages in case of sparse windows or median in rightmost histogram bin')
		else:
			#PeakPickerCWT
			pp	= PeakPickerCWT()
			param = pyopenms.Param()
			param.setValue('signal_to_noise', individual[21], ('Minimal signal to noise ratio '
															   'for a peak to be picked.'))
			param.setValue('peak_width', float(individual[22]), 'Approximate fwhm of the peaks.')

		pp.setParameters(param)
		pp.pickExperiment(parsed_data, parsed_data_picked)

	@staticmethod
	def peak_filtering(individual, parsed_data=None, parsed_data_picked=None):
		if individual[12] < 1 or individual[16] < 1:
			return

		if individual[12] < 2:
			filter = ThresholdMower()
			param  = pyopenms.Param()
			param.setValue('threshold', individual[13], ('Intensity threshold, peaks below '
														 'this threshold are discarded'))
			filter.setParameters(param)
		else:
			filter = NLargest()
			param  = pyopenms.Param()
			param.setValue('n', int(round(individual[14])), 'The number of peaks to keep')
			filter.setParameters(param)

		Filters._apply_peak(filter, parsed_data_picked)


def insybio_automatic_pipeline_construction(input_data_folder, reference_peptides_filename,
											population, generations, semi_random_initialization_probability,
											two_points_crossover_probability, arithmetic_crossover_probability,
											mutation_probability, gaussian_mutation_variance_proportion, min_values,
											max_values, goal_significances, initial_good_solutions, database_name, missing_values_filename,output_folder ):
	# input_data_folder: (String) it is the data folder including the mzML which will be used during the evaluation phase
	# reference_peptides_filename: (String) The filename of the reference peptides peaks.
	#							  Its format should be tab delimited with a column being constituted from
	#							  the name of the peptides and the mz values for its reference peaks
	# population:  (integer) (default: 10) it is the number of individual solutions which are evolved on parallel
	# generations: (integer) (default: 10) it is the maximum nulber of generations which we allow for the population to get evolved
	# semi_random_initialization_probability: (float) (default: 0.8) it is the proportion (multiplied with 100) of the individuals
	#										 which will be initialized with the semi automatic method. The rest will be initialized using random method
	# two_points_crossover_probability: (float) (default: 0.45) it is the probability (multiplied with 100) for using two points crossover operator
	# arithmetic_crossover_probability: (float) (default: 0.45) it is the probability (multiplied with 100) for using  arithmetic crossover operator
	# mutation_probability: (float) (default: 0.05) it is the probability for applying gaussian mutation operator
	# gaussian_mutation_variance_proportion: (float) (default:0.1) it is the proportion of the abs(max_val-min_val) which is used as variance of the gaussian mutation operator
	# min_values: (list) it is a list with the minimum allowed values for the variables which should be optimized
	# max_values: (list) it is a list with the maximum allowed values for the variables which should be optimized
	# goal_significances: (list of floats) Includs the significances (values from 0 to 1) of the individual goals
	# initial good solutions: (String): Name of the file which includes encoded good solutions
	# database_name: (String) the string of the name of the fasta database to be used for scaffold
	# missing_values_filename: (String) the string denoting the name of the file describing a priori known missing channels

	# Initialize Population of Individual Solutions
	tstart = time.time()
	tstamp = time.strftime('%Y_%m_%d')
	individuals = init_individuals_step1(population, semi_random_initialization_probability, min_values, max_values, initial_good_solutions)
	
	best_performance_fid	= open(output_folder+'best_performance_'+str(tstamp)+'.txt', 'w')
	average_performance_fid = open(output_folder+'average_performance_'+str(tstamp)+'.txt', 'w')
	best_solutions_fid	  = open(output_folder+'best_solutions_'+str(tstamp)+'.txt', 'w')

	best_solution_fid	   = open(output_folder+'best_solution_'+str(tstamp)+'.txt', 'w')
	final_solutions_fid	 = open(output_folder+'final_solutions_'+str(tstamp)+'.txt', 'w')
	missing_values_fid	  = open(missing_values_filename,'r')
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
		
		
	selected_individuals = [[0 for x in range(len(min_values))] for x in range(population)]

	#Apply Evolutionary Process to Optimize Individual Solutions
	print '\n\nRunning Optimization ...'

	for rep in range(generations):
		print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print '\nGeneration %s' % (rep+1)
		print '------------\n'
		print ('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
		print 'Individuals: '

		for i, v in enumerate(individuals):
			print '\n## %s' % (i+1)
			print ', '.join(map(str, v))
		print '\n'

		evaluation_values = evaluate_individuals_step1(input_data_folder, reference_peptides_filename,
												 individuals, goal_significances, database_name,missing_values_list)


		number_of_solutions = len(individuals)

		fronts = get_fronts(number_of_solutions, population, evaluation_values)

		print 'Pareto Frontiers are %s' % fronts

		#apply selection operator: Ranked base selection is used
		#find and write to file maximum and average performances
		eval_slice	= evaluation_values[-1][:population]
		max_eval	  = max(eval_slice)
		max_position2 = (eval_slice).index(max_eval)
		sum_eval	  = sum(eval_slice)

		best_performance_fid.write('%s\n' % max_eval)
		average_performance_fid.write('%s\n' % (sum_eval/population))

		#Tune fitness values by locating and using solution niches
		sigma_share = 0.5/(float(len(individuals[0]))**(0.1))
		for i in range(1, max(fronts)+1):
			ind = [y for y, x in enumerate(fronts) if x == i]

			#Calculate max values per goal per pareto frontier
			max_significances = [0]*len(goal_significances)

			for i in ind:
				for goal in range(len(goal_significances)):
					if evaluation_values[goal][i]>=max_significances[goal]:
						max_significances[goal]=evaluation_values[goal][i]

			for i in ind:
				m = 0
				for j in ind:
					d = 0
					for gene in range(len(individuals[0])):
						d = d + ((individuals[i][gene] - individuals[j][gene])/float(max_values[gene] - min_values[gene]))**2

					d = math.sqrt(d/(len(individuals[0])))

					if d<= sigma_share:
						m = m+(1-((d/float(sigma_share))**2))

				if m == 0: m = 1

				for goal in range(len(goal_significances)):
					evaluation_values[goal][i] = float(max_significances[goal])/m

		print 'evaluation_values after niche calculations : '
		for i, v in enumerate(evaluation_values):
			print '> %s' % (i+1)
			print ', '.join(map(str, v))
		print '\n'

		# doing some resetting on last evaluation_values
		for i in range(len(evaluation_values[0])):
			evaluation_values[-1][i] = 0

			for j in range(len(evaluation_values)-1):
				evaluation_values[-1][i] = evaluation_values[-1][i] + evaluation_values[j][i]

			evaluation_values[-1][i] = evaluation_values[-1][i]/float(len(evaluation_values)-1)

		eval_slice   = evaluation_values[-1][:population]
		max_eval	 = max(eval_slice)
		max_position = (eval_slice).index(max_eval)
		sum_ranked   = sum(eval_slice)

		sum_prop=[0]*(population+1)
		for i in range(1, population+1):
			sum_prop[i] = sum_prop[i-1] + evaluation_values[-1][i-1]/float(sum_ranked) ##Check this out

		for i in range(1, population):
			random_number = random.uniform(0,1)
			for j in range(0, population):
				if random_number >= sum_prop[j] and random_number < sum_prop[j+1]:
					for k in range(0,len(individuals[0])):
						selected_individuals[i][k] = individuals[j][k]

		for k in range(len(individuals[0])):
			selected_individuals[0][k] = individuals[max_position2][k]

		best_solutions_fid.write('%s\n' % ('\t'.join(map(str, individuals[max_position2]))))

		#apply crossover operator
		for i in range(1, population-1, 2):
			if random_number < two_points_crossover_probability:
				#Two Point Crossover
				cross_point1 = 0
				cross_point2 = 0
				while cross_point1 == cross_point2:
					cross_point1 = math.ceil((len(individuals[0])-4)*random.uniform(0,1))
					if cross_point1 < math.floor((2*len(individuals[0])-1)/3):
						width = math.ceil(random.uniform(0,1)*(math.floor(len(individuals[0])-1)/3 -2))
					else:
						width = math.ceil(random.uniform(0,1))*(math.floor(len(individuals[0])/3 -1)-2-(cross_point1-math.floor(2*len(individuals[0])/3)))
					cross_point2 = cross_point1 + width

				# because math.floor and math.ceil return floats
				cross_point1 = int(cross_point1)
				cross_point2 = int(cross_point2)

				#SWAPPING
				if cross_point1 > cross_point2:
					cross_point1, cross_point2 = cross_point2, cross_point1

				for j in range(cross_point1, cross_point2 + 1):
					selected_individuals[i][j], selected_individuals[i+1][j] = selected_individuals[i+1][j], selected_individuals[i][j]


			elif random_number >= two_points_crossover_probability and random_number < (two_points_crossover_probability+arithmetic_crossover_probability):
				#arithmetic crossover
				alpha = random.uniform(0,1)
				for j in range(0,len(individuals[0])):
					(selected_individuals[i][j], selected_individuals[i+1][j]) = (alpha * selected_individuals[i][j] + (1-alpha) * selected_individuals[i+1][j],
																				  (1-alpha) * selected_individuals[i][j] + alpha * selected_individuals[i+1][j])

		#apply mutation operator
		for i in range(1, population):
			for j in range(0, len(individuals[0])):
				if random.uniform(0,1) < mutation_probability:
					print '~~ Mutated_gene: %s' % j
					selected_individuals[i][j] = selected_individuals[i][j] + random.gauss(0,0.1*(max_values[j]-min_values[j]))

			#Correct values out of boundaries
			for j in range(0, len(min_values)):
				if selected_individuals[i][j] < min_values[j]:
					selected_individuals[i][j] = min_values[j]
				if selected_individuals[i][j] > max_values[j]:
					selected_individuals[i][j] = max_values[j]

		#check termination criteria (more elaborate termination criteria should be added)

		#update the population with the offspings
		individuals = copy.deepcopy(selected_individuals)


	best_solutions_fid.close()
	best_performance_fid.close()
	average_performance_fid.close()

	#TODO i have a feeling that it should be max_position index INSYBIO REPLY: NO THIS IS CORRECT
	#best_solution_fid.write('\t'.join(individuals[max_position]))
	best_solution_fid.write('\t'.join(map(str, individuals[0])))
	best_solution_fid.close()

	#Write all final solutions in a file
	member_line = ['\t'.join(map(str, member)) for member in individuals]
	final_solutions_fid.write('\n'.join(member_line))
	final_solutions_fid.close()


	print '==============='
	print '#   THE END   #'
	print '==============='
	print '\n exec_time in %.02fs' % (time.time() - tstart)
	return [output_folder+'best_solution.txt']


def init_individuals_step1(population, semi_random_initialization_probability, min_values, max_values, initial_good_solutions):
	#population:(integer) (default: 10) it is the number of individual solutions which are evolved on parallel
	#semi_random_initialization_probability: (float) (default: 0.8) it is the proportion (multiplied with 100) of the individuals which will be initialized with the semi automatic method. The rest will be initialized using random method
	#min_values: (list) it is a list with the minimum allowed values for the variables which should be optimized
	#max_values: (list) it is a list with the maximum allowed values for the variables which should be optimized
	#initial_good_solutions: (String) the filename of a file including some encoded good solutions to be used on the initial population of solutions

	individuals = []

	with open(initial_good_solutions) as tsv:
		for line in csv.reader(tsv, delimiter='\t'):
			individuals.append([float(val) for val in line])

	while(len(individuals) < population):
		flip_coin = random.uniform(0, 1) < semi_random_initialization_probability
		individuals.append(Utils.generate_random_list(flip_coin, min_values, max_values))

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

def quantify(mz_values,missing_values_raw):
	intervals = [125.5, 126.5, 127.5, 128.5, 129.5, 130.5, 131.5]
	count = count_intervals(mz_values, intervals)
	intervals.remove(125.5)

	match_ions_count = 0
	for val in range(len(intervals)):
		if count[intervals[val]] >= 1 and missing_values_raw[val]=='0':
			match_ions_count += 1
	num_of_quantified_peaks=0;
	num_of_non_missing_channels=0;
	for i in range(6):
		if missing_values_raw[i]=='0':
			num_of_quantified_peaks=num_of_quantified_peaks+count[126.5+i];
			num_of_non_missing_channels+=1;
	variation = abs(6.0 - num_of_quantified_peaks)
	#print('quantify message 2')
	return (
		match_ions_count >= 0.5*num_of_non_missing_channels, # quantified
		match_ions_count == num_of_non_missing_channels, # fully quantified
		variation
	)


def get_fronts(number_of_solutions, population, evaluation_values):
	#Estimate non dominated fronts
	assigned = 0
	fronts   = [0]*population
	front	= 1

	eval_temp  = copy.deepcopy(evaluation_values)

	def _dominate(solution1, solution2):
		check	 = 2
		ffs	   = len(solution1)
		dominate1 = 1
		equal1	= 1
		f		 = 0
		while f < ffs and dominate1 == 1:
			if solution1[f] > solution2[f]:
				equal1 = 0
			elif solution1[f] != solution2[f]:
				dominate1 = 0
			f += 1

		if dominate1 == 1 and equal1 == 0:
			check = 1
		elif dominate1 == 1 and equal1 == 1:
			check = 2
		else:
			dominate2 = 1
			equal2	= 1
			f		 = 0
			while f < ffs and dominate2 == 1:
				if solution2[f] > solution1[f]:
					equal2=0
				elif solution2[f] != solution1[f]:
					dominate2 = 0
				f += 1
			if dominate2 == 1 and equal2 == 0:
				check = 3
		return check


	while assigned < population:
		non_dominated_solutions = [0]*number_of_solutions

		index = range(len(eval_temp[0]))

		eval_temp_index = zip(eval_temp[0],index)
		ordered_list = sorted(range(len(eval_temp_index)), key = lambda k: eval_temp_index[k],reverse = True)
		non_dominated_solutions[0] = ordered_list[0]
		number_of_non_dominated_solutions = 1

		for i in range(1, number_of_solutions):
			n		  = 0
			condition  = 0
			condition2 = 1
			while n < number_of_non_dominated_solutions and condition == 0:
				solution1 = [0]*(len(eval_temp)-1)
				solution2 = [0]*(len(eval_temp)-1)
				for j in range(len(eval_temp)-1):
					solution1[j] = eval_temp[j][ordered_list[i]]
					solution2[j] = eval_temp[j][ordered_list[n]]

				check = _dominate(solution1,solution2)

				if check == 3:
					condition  = 1
					condition2 = 0
				elif check == 1:
					if number_of_non_dominated_solutions == 1:
						condition = 1
						non_dominated_solutions[0] = ordered_list[i]
					else:
						number_of_non_dominated_solutions = number_of_non_dominated_solutions-1
						del non_dominated_solutions[n]
				n += 1
			if condition2 == 1:
				non_dominated_solutions[number_of_non_dominated_solutions] = ordered_list[i]
				number_of_non_dominated_solutions = number_of_non_dominated_solutions+1

		sorted_non_dominated_solutions = sorted(non_dominated_solutions, reverse = True)

		for i in range(number_of_non_dominated_solutions):
			assigned += 1
			fronts[sorted_non_dominated_solutions[i]] = front
			for j in range(len(eval_temp)):
				eval_temp[j][sorted_non_dominated_solutions[i]] = 0

		front += 1

	return fronts

def evaluate_individuals_step1(input_data_folder, reference_peptides_filename, individuals, goal_significances, database_name, missing_values_list):
	#input_data_folder:(String) it is the data folder including the mzML which will be used during the evaluation phase
	#reference_peptides_filename: (String) The filename of the reference peptides peaks.
	#							 Its format should be tab delimited with a column being constituted from
	#							 the name of the peptides and the mz values for its reference peaks
	#individuals: (two dimensional list) including the population of individual solutions
	#goal_significances: (list of floats) Includes the significances (values from 0 to 1) of the individual goals
	#database_name: (String) the string of the name of the fasta database to be used for scaffold
	#missing_values_list: (String) the list containing information about a priori missing channels

	#load mzML files for evaluation
	if not os.path.exists(input_data_folder):
		raise Exception('input directory [%s] for mzML does not exist'%input_data_folder)

	random_mzML_name = random.choice([name for name
						 in os.listdir(input_data_folder)
						 if name.endswith('.mzML') ])
	missing_values_experiment=list()
	for i in range(len(missing_values_list)):
		if random_mzML_name in missing_values_list[i][0]:
			for j in xrange(1, len(missing_values_list[0])):
				missing_values_experiment.append(missing_values_list[i][j])
		
	random_mzML_path = os.path.join(input_data_folder, random_mzML_name)
	exp = MSExperiment()
	file = pyopenms.FileHandler()
	file.loadExperiment(random_mzML_path, exp)
	print 'Input file was successfully parsed and loaded [%s]'%random_mzML_name

	reference_peptides_peaks = []
	with open(reference_peptides_filename) as tsv:
		for line in csv.reader(tsv, delimiter='\t'):
			reference_peptides_peaks.append([float(val) for val in line[1:]])

	print('Missing values row')
	print(missing_values_experiment)
	#create a list of Individuals objects
	obj_individuals = []
	for i, individual in enumerate(individuals):
		obj_individuals.append(Individual_step1(individual, i,os.path.join(os.getcwd(), 'individuals') , reference_peptides_peaks, goal_significances,random_mzML_path,missing_values_experiment))


	evaluation_values = [[] for i in range(8)]
	pool = mp.Pool(processes=1)
	results = pool.map(evaluate_individual_step1, obj_individuals, chunksize=1)
	pool.close()
	pool.join()

	return [list(val) for val in zip(*results)]


def evaluate_individual_step1(individual):
	return individual.evaluate()


class Individual_step1():

	def __init__(self, individual, i, output_dir, reference_peptides_peaks, goal_significances,random_mzML_path, missing_values_raw):
		self.goal_significances = goal_significances
		self.individual = individual
		self.i = i
		self.output_dir = os.path.join(output_dir, str(i))
		self.filtered_dir = os.path.join(self.output_dir, 'filtered_files')
		self.reference_peptides_peaks = reference_peptides_peaks
		self.random_mzML_path = random_mzML_path
		self.missing_values_raw = missing_values_raw

	@staticmethod
	def pick_experiment(zfile):
		ms = MSExperiment()
		pp = PeakPickerHiRes()
		pp.setParameters(pyopenms.Param())
		pp.pickExperiment(zfile, ms)
		return ms

	@staticmethod
	def individual_mz_matrix(individual, random_mzML_path, path_for_filtered_files):
		filtered_exp = filter_experiments_new(random_mzML_path, individual)
		file = pyopenms.FileHandler()
		print('testi1')
		file.storeExperiment(path_for_filtered_files, filtered_exp)
		print('testi2')
		if individual[16] < 1:
			filtered_exp = Individual.pick_experiment(filtered_exp)

		filtered_msms = [f for f in filtered_exp if f.getMSLevel() != 1]
		print("test1")
		mz_matrix = []
		for f in filtered_msms:

			peaks = f.get_peaks()
			if len(peaks):
				mz_values_old, intensities_old = zip(*peaks)
				mz_matrix.append([float(v) for v in mz_values_old])
			else:
				mz_matrix.append([])
		return mz_matrix

	def evaluate(self):
		print '\n\n> Evaluating individual ## %s ################################' % (self.i)
		print('testk1')
		if not os.path.exists(self.output_dir):
			os.makedirs(self.output_dir)

		if not os.path.exists(self.filtered_dir):
			os.makedirs(self.filtered_dir)
		else:
			Utils.remove_files_in_dir(self.filtered_dir)
		
		evaluation_values = [[] for i in range(8)]
		evaluation_values = []
		filtered_file = os.path.join(self.filtered_dir, os.path.basename(self.random_mzML_path))
		mz_matrix = Individual_step1.individual_mz_matrix(self.individual, self.random_mzML_path, filtered_file)
		print('testk2')
		# In case all spectra are filterd out provide to this solution a very small fitness value.
		if not len(mz_matrix):
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			evaluation_values.append(0.001)
			return evaluation_values

		#This variable represents the number of MSMS spectra
		ms_ms2 = len(mz_matrix)
		(quantified, fully_quantified, reporter_ions_variation) = (0, 0, 0)
		maximum_value_per_peptide = []
		#print(str(self.i)+'Test message 1')
		t0 = time.time()
		for i, mz_values in enumerate(mz_matrix):
			if i % 2000 == 0:
				print '>> Spectra analysed: %s' % i

			(q, fq, variation) = quantify(mz_values,self.missing_values_raw)
			t1 = time.time()
			if q:
				quantified += 1
			if fq:
				fully_quantified += 1
			reporter_ions_variation += variation

			for line in self.reference_peptides_peaks:
				found = 0
				for val in line:
					#A window of 0.02 mz values is used to searche nearby the reference peaks.
					dictionary_result = count_intervals(
							mz_values, [
								val - 0.02,
								val + 0.02
							]
					)
					if dictionary_result[val + 0.02] >= 1:
						found += 1

				maximum_value_per_peptide.append((found / float(len(line))))

		print '<< Spectra analysis done in %.02fs' % (time.time() - t0)
		reference_peak_similarity = sum(maximum_value_per_peptide)/float(len(maximum_value_per_peptide))
		reporter_ions_variation   = reporter_ions_variation/float(ms_ms2)
		fully_quantified		  = fully_quantified/float(ms_ms2)
		quantified				= quantified/float(ms_ms2)


		################
		# xTandem Search
		################
		# empty input/output folders
		xt_input_dir = os.path.join(self.output_dir, 'xtandem_input')
		xt_output_dir = os.path.join(self.output_dir, 'xtandem_output')
		for path in [xt_input_dir, xt_output_dir]:
			if not os.path.exists(path):
				os.makedirs(path)
			else:
				Utils.remove_files_in_dir(path)

		xt_input_file = os.path.join(self.output_dir, 'input.xml')
		xt_output_dir = os.path.join(self.output_dir, 'xtandem_output')
		xt_output_file = os.path.join(xt_output_dir, 'output.xml')
		#create input.xml file
		#To check if name of file is correct
		print '\n >> xTandem Search'
		print 'file: %s' % (self.random_mzML_path)
		with open (xt_input_file, 'w') as input_xml_fid:
			input_xml_fid.write('<?xml version="1.0"?>\n')
			input_xml_fid.write('<bioml>\n')
			input_xml_fid.write('\t<note>\n')
			input_xml_fid.write('\tEach one of the parameters for x! tandem is entered as a labeled note node.\n')
			input_xml_fid.write('\tAny of the entries in the default_input_step1.xml file can be over-ridden by\n')
			input_xml_fid.write('\tadding a corresponding entry to this file. This file represents a minimum\n')
			input_xml_fid.write('\tinput file, with only entries for the default settings, the output file\n')
			input_xml_fid.write('\tand the input spectra file name.\n')
			input_xml_fid.write('\tSee the taxonomy.xml file for a description of how FASTA sequence list\n')
			input_xml_fid.write('\tfiles are linked to a taxon name.\n')
			input_xml_fid.write('\t</note>\n')
			input_xml_fid.write('\t<note type="input" label="list path, default parameters">default_input_step1.xml</note>\n')
			input_xml_fid.write('\t<note type="input" label="list path, taxonomy information">taxonomy.xml</note>\n')
			input_xml_fid.write('\t<note type="input" label="protein, taxon">human</note>\n')
			input_xml_fid.write('\t<note type="input" label="spectrum, path">%s</note>\n' % filtered_file)
			input_xml_fid.write('\t<note type="input" label="output, path">%s</note>\n' % xt_output_file)
			input_xml_fid.write('</bioml>\n')

		command = ['./tandem.exe', xt_input_file]
		result_p = Popen(command)

		#reading xTandem results
		result_p.communicate()

		(identified_xtandem, identified_xtandem2) = (0, 0)

		#parse  the xtandem result file
		base, ext = os.path.splitext(xt_output_file)
		xt_result_file = glob.glob(base + '*' +  ext)[0]
		proteins=list()
		with open(xt_result_file, 'r') as output_xml_fid:
			for line in output_xml_fid:
				xtandem1 = re.match(r'\s*<note.*label.*duplicate.peptide.ids[^>]*>([^<]+)</note>', line)
				xtandem2 = re.match(r'\s*<note.*label.*total.unique.assigned[^>]*>([^<]+)</note>', line)
				if line.find("<protein expect=")!=-1:
					words=line.split("label=\"")
					words2=words[1].split("\"")
					try:
						index_in_proteins=proteins.index(words2[0])
					except ValueError:
						index_in_proteins=-1
					if index_in_proteins==-1:
						proteins.append(words2[0])
				if xtandem1:
					identified_xtandem  = int(xtandem1.group(1))
				if xtandem2:
					identified_xtandem2 = int(xtandem2.group(1))


		identified_xtandem  = identified_xtandem  / float(ms_ms2)
		identified_xtandem2 = identified_xtandem2 / float(ms_ms2)
		

		goal1 = 1/(0.1+reporter_ions_variation)
		goal2 = reference_peak_similarity
		goal3 = quantified
		goal4 = fully_quantified
		goal5 = identified_xtandem
		goal6 = identified_xtandem2
		goal7 = len(proteins)/float(500) 

		print '--------------------------------------------------------------'
		print '|	 |  1  |  Reporter Ions Difference	 | %s ' % goal1
		print '|  G  |  2  |  Reference Peptides Peak	  | %s ' % goal2
		print '|  O  |  3  |  Quantified				   | %s ' % goal3
		print '|  A  |  4  |  Fully Quantified			 | %s ' % goal4
		print '|  L  |  5  |  Xtandem Identified Peptides  | %s ' % goal5
		print '|  S  |  6  |  Xtandem Matched Spectra	  | %s ' % goal6
		print '|	 |  7  |  Scaffold Identified Proteins | %s ' % goal7
		print '--------------------------------------------------------------'

		evaluation_values.append(goal1)
		evaluation_values.append(goal2)
		evaluation_values.append(goal3)
		evaluation_values.append(goal4)
		evaluation_values.append(goal5)
		evaluation_values.append(goal6)
		evaluation_values.append(goal7)
		evaluation_values.append(self.goal_significances[0] * goal1 +
								 self.goal_significances[1] * goal2 +
								 self.goal_significances[2] * goal3 +
								 self.goal_significances[3] * goal4 +
								 self.goal_significances[4] * goal5 +
								 self.goal_significances[5] * goal6 +
								 self.goal_significances[6] * goal7)

		return evaluation_values

def filter_experiments_new(mzML_file, individual):
	#Individual (list): Representing an individual solution
	#parsed_data (list): a list of MSExperiment storing parsed mzML files

	t0 = time.time()

	parsed_data_picked = MSExperiment()
	parsed_data		= MSExperiment()
	print('testn1')
	file = pyopenms.FileHandler()
	print(mzML_file)
	file.loadExperiment(mzML_file, parsed_data)
	print('testn2')
	if individual[0] > individual[6] and individual[6] > individual[9] and individual[9] > individual[15]:
		Filters.apply(
			['denoise', 'precursor_removal', 'normalizer', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[0] >= individual[6] and individual[6] >= individual[15] and individual[15] >= individual[9]:
		Filters.apply(
			['denoise', 'precursor_removal', 'peak_picking', 'normalizer', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[0] >= individual[9] and individual[9] >= individual[6] and individual[6] >= individual[15]:
		Filters.apply(
			['denoise', 'normalizer', 'precursor_removal', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[0] >= individual[9] and individual[9] >= individual[15] and individual[15] >= individual[6]:
		Filters.apply(
			['denoise', 'normalizer', 'peak_picking', 'precursor_removal', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[0] >= individual[15] and individual[15] >= individual[9] and individual[9] >= individual[6]:
		Filters.apply(
			['denoise', 'peak_picking', 'normalizer', 'precursor_removal', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[0] >= individual[15] and individual[15] >= individual[6] and individual[6] >= individual[9]:
		Filters.apply(
			['denoise', 'peak_picking', 'precursor_removal', 'normalizer', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[6] >= individual[0] and individual[0] >= individual[9] and individual[9] >= individual[15]:
		Filters.apply(
			['precursor_removal', 'denoise', 'normalizer', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[6] >= individual[0] and individual[0] >= individual[15] and individual[15] >= individual[9]:
		Filters.apply(
			['precursor_removal', 'denoise', 'peak_picking', 'normalizer', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[9] >= individual[0] and individual[0] >= individual[6] and individual[6] >= individual[15]:
		Filters.apply(
			['normalizer', 'denoise', 'precursor_removal', 'peak_picking', 'peak_filtering' ],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[9] >= individual[0] and individual[0] >= individual[15] and individual[15] >= individual[6]:
		Filters.apply(
			['normalizer', 'denoise', 'peak_picking', 'precursor_removal', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[6] >= individual[9] and individual[9] >= individual[0]:
		Filters.apply(
			['precursor_removal', 'normalizer', 'denoise', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	elif individual[9] >= individual[6] and individual[6] >= individual[0]:
		Filters.apply(
			['normalizer', 'precursor_removal', 'denoise', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	else:
		Filters.apply(
			['denoise', 'precursor_removal', 'normalizer', 'peak_picking', 'peak_filtering'],
			individual, parsed_data=parsed_data, parsed_data_picked=parsed_data_picked
		)
	print '<< Filters done in %.02fs\n' % (time.time() - t0)
	return parsed_data if individual[16] < 1 else parsed_data_picked


def apply_best_solution(input_folder,best_solution_file,quantified_spectra_thres, missing_values_filename, output_folder):
	best_solution_fid = open(best_solution_file,"r")
	number_of_lines=0
	best_solution=list()
	missing_values_fid	  = open(missing_values_filename,'r')
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
			filename_out = output_folder+"quantified_files/quantified_filtered_"+str(filename)
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
		print("test1")
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
		print("test2")
		print(individual)
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
		print("test3")
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

#Step 2 functions 
			
def quantification_and_spectra_alignment(input_folder,max_allowed_variance_in_peak_number,minimum_similarity_score,retention_time_tolerance,quantified_spectra_thres,abs_sim_thres,min_first_to_second_score_distance, missing_values_filename, output_folder):
	missing_values_fid	  = open(missing_values_filename,'r')
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
	maximum_number_of_spectra=60000
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
				for k in range((60000-parsed_file.size())):
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
				##Precursor Mass Filter start
				exp_spectra_list_precursor_mass=list()
				##Precursor Mass Filter end
				for sp in xrange(0,num_of_spectra-1):
					exp_spectra_list_rt.append(exp_spectra_list[sp].getRT())
					##Precursor Mass Filter start
					precursor=exp_spectra_list[sp].getPrecursors()
					exp_spectra_list_precursor_mass.append(precursor[0].getUnchargedMass())	
					##Precursor Mass Filter end
				for k in range(parsed_file.size()):
					parsed_file_k_MSLevel = parsed_file[k].getMSLevel()
					parsed_file_k_peaks = parsed_file[k].get_peaks()
					parsed_file_k_RT = parsed_file[k].getRT()
					##Precursor Mass Filter start
					if parsed_file[k].getMSLevel()==1:
						parsed_file_k_precursor_mass=-1
					else:
						precursors=parsed_file[k].getPrecursors()
						parsed_file_k_precursor_mass=precursors[0].getMZ()
					##Precursor Mass Filter end
					obj_searchparams.append(Searchparams(k,
						parsed_file_k_MSLevel,
						parsed_file_k_peaks,
						parsed_file_k_RT,
						##Precursor Mass Filter start
						parsed_file_k_precursor_mass,
						##Precursor Mass Filter end
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
						##Precursor Mass Filter start
						exp_spectra_list_precursor_mass,
						##Precursor Mass Filter end
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
		mapping_spectra_file_fid = open(output_folder+"mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
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
		similarities_mapping_spectra_file_fid = open(output_folder+"similarities_mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt","w")
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
		quantitative_values_fid = open(output_folder+"quantitative_values"+str(time.strftime("%Y_%m_%d"))+".txt","w")
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
		file.store(output_folder+"unified_spectrum_list_"+str(time.strftime("%Y_%m_%d"))+".mzML",exp_spectra_list)		
		print(str(number_of_parsed_files)+" mzML files were successfully filtered and stored")
	return [output_folder+"unified_spectrum_list_"+str(time.strftime("%Y_%m_%d"))+".mzML",output_folder+"mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt",output_folder+"quantitative_values"+str(time.strftime("%Y_%m_%d"))+".txt",output_folder+"similarities_mapping_spectra_file_"+str(time.strftime("%Y_%m_%d"))+".txt"]

def searchparams_spectrum_in_unified_list(searchparams):
	#print searchparams.parsed_file_k_MSLevel
	return searchparams.search_spectrum_in_unified_list()

class Searchparams():
	def __init__(self, k, parsed_file_k_MSLevel,parsed_file_k_peaks,parsed_file_k_RT, parsed_file_k_precursor_mass, num_of_spectra, retention_time_tolerance, max_allowed_variance_in_peak_number, intensities_sum, c, dot_product, minimum_similarity_score, spectra_found_list, arg_fabs, exp_spectra_list_rt, exp_spectra_list_precursor_mass, spectra_peaks, missing_values_raw,quantified_spectra_thres,abs_sim_thres, min_first_to_second_score_distance):
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
		##Precursor Mass Filter start
		self.exp_spectra_list_precursor_mass=exp_spectra_list_precursor_mass
		##Precursor Mass Filter end
		self.spectra_peaks = spectra_peaks
		self.parsed_file_k_MSLevel = parsed_file_k_MSLevel
		self.parsed_file_k_peaks = parsed_file_k_peaks
		self.parsed_file_k_RT = parsed_file_k_RT
		##Precursor Mass Filter start
		self.parsed_file_k_precursor_mass=parsed_file_k_precursor_mass
		##Precursor Mass Filter end
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
						if self.arg_fabs(self.exp_spectra_list_rt[sp]-RT_time_current_spectra)<self.retention_time_tolerance and ((self.arg_fabs(self.parsed_file_k_precursor_mass-self.exp_spectra_list_precursor_mass[sp])< self.minimum_similarity_score) or (self.arg_fabs(self.parsed_file_k_precursor_mass-((self.exp_spectra_list_precursor_mass[sp]+2*1.00783)/2))<2*self.minimum_similarity_score) or (self.arg_fabs(self.parsed_file_k_precursor_mass-((self.exp_spectra_list_precursor_mass[sp]+3*1.00783)/3))<2*self.minimum_similarity_score) or (self.arg_fabs(self.parsed_file_k_precursor_mass-((self.exp_spectra_list_precursor_mass[sp]+4*1.00783)/4))<2*self.minimum_similarity_score) or (self.arg_fabs(self.parsed_file_k_precursor_mass-((self.exp_spectra_list_precursor_mass[sp]+5*1.00783)/5))<2*self.minimum_similarity_score) or (self.arg_fabs(self.parsed_file_k_precursor_mass-((self.exp_spectra_list_precursor_mass[sp]+6*1.00783)/6))<2*self.minimum_similarity_score) ):
							if self.spectra_found_list[sp]<1:
								len2=len(self.spectra_peaks[sp][0])
								if len2>8 and self.arg_fabs(len1-len2)<self.max_allowed_variance_in_peak_number :
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
#Step 3 functions
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
	missing_values_fid2	  = open(missing_values_filename,'r')
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
	individuals=initialize_step3(population,semi_random_initialization_probability,min_values, max_values)
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
		evaluation_values=evaluate_individuals_step3(quantitative_spectra_dataset,quantitative_spectra_dataset_normalized,unified_spectra,individuals,goal_significances,database_name,spectraST_databases_filename,num_of_folds,used_spectra,classification_problems, number_of_used_tmt_channels, supervised_mode,output_folder)		
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
	input_xml_fid.write("\tAny of the entries in the default_input_step3.xml file can be over-ridden by\n")
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
	input_xml_fid.write("\t<note type=\"input\" label=\"list path, default parameters\">default_input_step3.xml</note>\n")
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
	return output_folder+"selected_normalized_missing_values_completed_"+str(time.strftime("%Y_%m_%d"))+".txt"

def initialize_step3(population,semi_random_initialization_probability,min_values, max_values):
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


def evaluate_individuals_step3(dataset,dataset_with_missing_values,unified_spectra,individuals,goal_significances,database_name,spectraST_databases_filename,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder):
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
		input_xml_fid.write("\tAny of the entries in the default_input_step3.xml file can be over-ridden by\n")
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
		input_xml_fid.write("\t<note type=\"input\" label=\"list path, default parameters\">default_input_step3.xml</note>\n")
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
		obj_individuals.append(Individual_step3(individual, i, dataset,dataset_with_missing_values,unified_spectra,goal_significances,num_of_folds,used_spectra,classification_problems,number_of_used_tmt_channels,supervised_mode,output_folder))
	pool = mp.Pool(processes=mp.cpu_count()-1)
	
	results = pool.map(evaluate_individual_step3, obj_individuals, chunksize=1)
	
	pool.close()
	pool.join()
	for i in range(len(results)):
		evaluation_values[6].append(float(results[i][0]))	
	for i in range(len(evaluation_values[0])):
		evaluation_values[7].append((float(evaluation_values[0][i])*float(goal_significances[0])+float(evaluation_values[1][i])*float(goal_significances[1])+float(evaluation_values[2][i])*float(goal_significances[2])+float(evaluation_values[3][i])*float(goal_significances[3])+float(evaluation_values[5][i])*float(goal_significances[5])+float(evaluation_values[6][i])*float(goal_significances[6]))/6.0)
	for i in range(len(evaluation_values[0])):
		print '-----------------Individual:'+str(i)+'----------------------------'
		print '--------------------------------------------------------------'
		print '|	 |  1  |  Minimizing selected by the classifier spetra | %s ' % evaluation_values[0][i]
		print '|	 |  2  |  Minimizing PTMs				   | %s ' % evaluation_values[1][i]
		print '|  G  |  3  |  Differentially Quantified Spectra			| %s ' % evaluation_values[2][i]
		print '|  O  |  4  |  Xtandem Identidied Peptides	  		   | %s ' % evaluation_values[3][i]
		print '|  A  |  5  |  Percentage of Validated Peptides Matches	   | %s ' % evaluation_values[4][i]
		print '|  L  |  6  |  Quantified Proteins			   | %s ' % evaluation_values[5][i]
		print '|  S  |  7  |  Classifiers Accuracy			   | %s ' % evaluation_values[6][i]
		print '|	 |  8  |  Weighted Sum of Goals 			   | %s ' % evaluation_values[7][i]
		print '--------------------------------------------------------------'	
	return evaluation_values
def evaluate_individual_step3(individual):
	return individual.evaluate()	
class Individual_step3():
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

#step4 functions

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
	return [output_folder+str(mode)+"/corrected_unified_spectrum_list_missing_values_imputed.mzML",output_folder+"new_quantified_files/"]

