#!/usr/bin/env python
# -*- coding: utf-8 -*-

####################################
#This code is the outcome of the project Innovative Processing of TMT Proteomics. This is a common project of InSyBio LTD and Nestle Institute of Health Sciences. 
#This program implements the step 2 code
####################################

import os
import random
import math
import copy
import time
import sys
import subprocess
import platform
import glob
from shutil import copyfile
import multiprocessing as mp
from pyopenms import *
from collections import defaultdict
from bisect import bisect_left
from subprocess import Popen
import warnings
import logging

###
import re
import csv
import pdb


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
                            random.gauss  (0.2,            0.1),
                            random.gauss  (10.0,           0.1*(max_values[3]-min_values[3])),
                            random.gauss  (11.0,           0.1*(max_values[4]-min_values[4])),
                            random.gauss  (4,              0.1*(max_values[5]-min_values[5])),
                            random.uniform(min_values[6],  max_values[6]),
                            random.uniform(min_values[7],  max_values[7]),
                            random.gauss  (2,              0.1*(max_values[5]-min_values[5])),
                            random.uniform(min_values[9],  max_values[9]),
                            random.uniform(min_values[10], max_values[10]),
                            random.gauss  (0.1,            0.1),
                            random.uniform(min_values[12], max_values[12]),
                            random.gauss  (0.05,           0.05),
                            random.gauss  (200,            0.1*(max_values[14]-min_values[14])),
                            random.uniform(min_values[15], max_values[15]),
                            random.uniform(min_values[16], max_values[16]),
                            random.gauss  (1,              0.1*(max_values[17]-min_values[17])),
                            random.gauss  (200,            0.1*(max_values[18]-min_values[18])),
                            random.gauss  (30,             0.1*(max_values[19]-min_values[19])),
                            random.gauss  (10,             0.1*(max_values[20]-min_values[20])),
                            random.gauss  (1,              0.1*(max_values[21]-min_values[21])),
                            random.gauss  (0.15,           0.15),
                            random.gauss  (0.25,           0.25) ]

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
            pp    = PeakPickerHiRes()
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
            pp    = PeakPickerCWT()
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
                                            max_values, goal_significances, initial_good_solutions, database_name, missing_values_filename ):
    # input_data_folder: (String) it is the data folder including the mzML which will be used during the evaluation phase
    # reference_peptides_filename: (String) The filename of the reference peptides peaks.
    #                              Its format should be tab delimited with a column being constituted from
    #                              the name of the peptides and the mz values for its reference peaks
    # population:  (integer) (default: 10) it is the number of individual solutions which are evolved on parallel
    # generations: (integer) (default: 10) it is the maximum nulber of generations which we allow for the population to get evolved
    # semi_random_initialization_probability: (float) (default: 0.8) it is the proportion (multiplied with 100) of the individuals
    #                                         which will be initialized with the semi automatic method. The rest will be initialized using random method
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
    individuals = init_individuals(population, semi_random_initialization_probability, min_values, max_values, initial_good_solutions)
    tstamp = time.strftime('%Y_%m_%d')
    output_folder='step1/'+str(tstamp)+'_'+str(time.time())+'/'
    if not os.path.exists(output_folder):
        os.makedirs(output_folder)
    
    best_performance_fid    = open(output_folder+'best_performance_'+str(tstamp)+'.txt', 'w')
    average_performance_fid = open(output_folder+'average_performance_'+str(tstamp)+'.txt', 'w')
    best_solutions_fid      = open(output_folder+'best_solutions_'+str(tstamp)+'.txt', 'w')

    best_solution_fid       = open(output_folder+'best_solution_'+str(tstamp)+'.txt', 'w')
    final_solutions_fid     = open(output_folder+'final_solutions_'+str(tstamp)+'.txt', 'w')
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

        evaluation_values = evaluate_individuals(input_data_folder, reference_peptides_filename,
                                                 individuals, goal_significances, database_name,missing_values_list)


        number_of_solutions = len(individuals)

        fronts = get_fronts(number_of_solutions, population, evaluation_values)

        print 'Pareto Frontiers are %s' % fronts

        #apply selection operator: Ranked base selection is used
        #find and write to file maximum and average performances
        eval_slice    = evaluation_values[-1][:population]
        max_eval      = max(eval_slice)
        max_position2 = (eval_slice).index(max_eval)
        sum_eval      = sum(eval_slice)

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
        max_eval     = max(eval_slice)
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


def init_individuals(population, semi_random_initialization_probability, min_values, max_values, initial_good_solutions):
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
    front    = 1

    eval_temp  = copy.deepcopy(evaluation_values)

    def _dominate(solution1, solution2):
        check     = 2
        ffs       = len(solution1)
        dominate1 = 1
        equal1    = 1
        f         = 0
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
            equal2    = 1
            f         = 0
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
            n          = 0
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

def evaluate_individuals(input_data_folder, reference_peptides_filename, individuals, goal_significances, database_name, missing_values_list):
    #input_data_folder:(String) it is the data folder including the mzML which will be used during the evaluation phase
    #reference_peptides_filename: (String) The filename of the reference peptides peaks.
    #                             Its format should be tab delimited with a column being constituted from
    #                             the name of the peptides and the mz values for its reference peaks
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
        if random_mzML_name==missing_values_list[i][0]:
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
        obj_individuals.append(Individual(individual, i,os.path.join(os.getcwd(), 'individuals') , reference_peptides_peaks, goal_significances,random_mzML_path,missing_values_experiment))


    evaluation_values = [[] for i in range(8)]
    pool = mp.Pool(processes=mp.cpu_count())
    results = pool.map(evaluate_individual, obj_individuals, chunksize=1)
    pool.close()
    pool.join()

    return [list(val) for val in zip(*results)]


def evaluate_individual(individual):
    return individual.evaluate()


class Individual():

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
        filtered_exp = filter_experiments(random_mzML_path, individual)
        file = pyopenms.FileHandler()
        file.storeExperiment(path_for_filtered_files, filtered_exp)

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
    
        if not os.path.exists(self.output_dir):
            os.makedirs(self.output_dir)

        if not os.path.exists(self.filtered_dir):
            os.makedirs(self.filtered_dir)
        else:
            Utils.remove_files_in_dir(self.filtered_dir)

        evaluation_values = [[] for i in range(8)]
        evaluation_values = []

       

        filtered_file = os.path.join(self.filtered_dir, os.path.basename(self.random_mzML_path))
        mz_matrix = Individual.individual_mz_matrix(self.individual, self.random_mzML_path, filtered_file)

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
        fully_quantified          = fully_quantified/float(ms_ms2)
        quantified                = quantified/float(ms_ms2)


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
            input_xml_fid.write('\tAny of the entries in the default_input.xml file can be over-ridden by\n')
            input_xml_fid.write('\tadding a corresponding entry to this file. This file represents a minimum\n')
            input_xml_fid.write('\tinput file, with only entries for the default settings, the output file\n')
            input_xml_fid.write('\tand the input spectra file name.\n')
            input_xml_fid.write('\tSee the taxonomy.xml file for a description of how FASTA sequence list\n')
            input_xml_fid.write('\tfiles are linked to a taxon name.\n')
            input_xml_fid.write('\t</note>\n')
            input_xml_fid.write('\t<note type="input" label="list path, default parameters">default_input.xml</note>\n')
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
        print '|     |  1  |  Reporter Ions Difference     | %s ' % goal1
        print '|  G  |  2  |  Reference Peptides Peak      | %s ' % goal2
        print '|  O  |  3  |  Quantified                   | %s ' % goal3
        print '|  A  |  4  |  Fully Quantified             | %s ' % goal4
        print '|  L  |  5  |  Xtandem Identified Peptides  | %s ' % goal5
        print '|  S  |  6  |  Xtandem Matched Spectra      | %s ' % goal6
        print '|     |  7  |  Scaffold Identified Proteins | %s ' % goal7
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

def filter_experiments(mzML_file, individual):
    #Individual (list): Representing an individual solution
    #parsed_data (list): a list of MSExperiment storing parsed mzML files

    t0 = time.time()

    parsed_data_picked = MSExperiment()
    parsed_data        = MSExperiment()

    file = pyopenms.FileHandler()
    file.loadExperiment(mzML_file, parsed_data)

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


if __name__ == '__main__':
    #Initialize random generator seed with the current local time in miliseconds
    random.seed()

    min_values = [0.0, 0.0, 0.1, 1.0,  5.0,  3.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.01, 0.0, 0.01,  50,    0.0, 0.0, 0.0, 1.0,   3.0,  1.0,  0.0, 0.0, 0.0]
    max_values = [1.0, 3.0, 2.0, 20.0, 20.0, 8.0, 1.0, 2.0, 5.0, 1.0, 5.0, 0.99, 3.0, 100.0, 300.0, 1.0, 3.0, 3.0, 400.0, 50.0, 20.0, 3.0, 1.0, 2.0]

    # sudo nohup python bs1_step1_v1.py
    #[ 1]   '/home/theofilatos/Documents/nihs_project/input_files/'
    #[ 2]   'reference_proteins_peptides_peaks.txt'
    #[ 3]   3
    #[ 4]   2
    #[ 5]   0.8
    #[ 6]   0.45
    #[ 7]   0.45
    #[ 8]   0.05
    #[ 9]   0.1
    #[10]   'initial_good_solutions.txt'
    #[11]   'goal_significances.txt'
    #[12]   'Swissprot_Human_BLA_2014_dec_08.fasta'
    #[13]   'missing_values.txt'


    
    #input_data_folder:(String) it is the data folder including the mzML which will be used during the evaluation phase
    input_data_folder                      = sys.argv[1]
    reference_peptides_filename            = sys.argv[2]
    population                             = int(sys.argv[3])
    generations                            = int(sys.argv[4])
    semi_random_initialization_probability = float(sys.argv[5])
    two_points_crossover_probability       = float(sys.argv[6])
    arithmetic_crossover_probability       = float(sys.argv[7])
    mutation_probability                   = float(sys.argv[8])
    gaussian_mutation_variance_proportion  = float(sys.argv[9])
    initial_good_solutions                 = sys.argv[10]
    goal_significances_filename            = sys.argv[11]
    database_name                          = sys.argv[12]
    missing_values_filename                = sys.argv[13]

    goal_significances  = []
    with open(goal_significances_filename) as tsv:
        for line in csv.reader(tsv, delimiter='\t'):
            for val in line:
                goal_significances.append(float(val))

    insybio_automatic_pipeline_construction(
            input_data_folder,
            reference_peptides_filename,
            population,
            generations,
            semi_random_initialization_probability,
            two_points_crossover_probability,
            arithmetic_crossover_probability,
            mutation_probability,
            gaussian_mutation_variance_proportion,
            min_values,
            max_values,
            goal_significances,
            initial_good_solutions,
            database_name,
            missing_values_filename)
