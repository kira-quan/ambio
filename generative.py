#!/usr/bin/env python
# This program takes in a set of reads and the template strand. For each read, where it aligns on the 
# placement strand, sequences of likely reads are produced from the template using a Gaussian Mixture
# Model. The read is then assigned to the position on the template where it is most likely

from __future__ import division
from scipy import stats
from read import Read, GenerativeRead, Alignment
from numpy import arange, array
from random import randint
from sklearn import mixture # Version 0.14.1
from ambio import get_base_num
import timeit
import numpy as np
import sys


def ambio(reads, training_reads, components):
	# VARIABLE DECLARATION
	# """
	# In genome changed index: 26 of base G to A
	# In genome changed index: 51 of base A to G
	# In genome changed index: 184 of base C to T
	# In repeat changed index: 22 of base C to T
	# In repeat changed index: 5 of base G to A

	# """
	# reads = []
	# template_list = [['G', 'G', 'C', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'A', 'A', 'T', 'A', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'G', 'C', 'A', 'A'], ['G', 'G', 'C', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'A', 'A', 'T', 'A', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'G', 'C', 'A', 'A'], ['A', 'G', 'G', 'C', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'A', 'G', 'G', 'C', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'G', 'G', 'A', 'G', 'G', 'C', 'A', 'C', 'G', 'G', 'C', 'G', 'A', 'A', 'G', 'A', 'A', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'C'], ['A', 'G', 'G', 'C', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'A', 'G', 'G', 'C', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'G', 'G', 'A', 'G', 'G', 'C', 'A', 'C', 'G', 'G', 'C', 'G', 'A', 'A', 'G', 'A', 'A', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'C'], ['C', 'C', 'A', 'C', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'C', 'G', 'C', 'G', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'C', 'G', 'A', 'C', 'A', 'G', 'G', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'C', 'G'], ['C', 'C', 'A', 'C', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'C', 'G', 'C', 'G', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'C', 'G', 'A', 'C', 'A', 'G', 'G', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'C', 'G'], ['G', 'A', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G', 'C', 'C', 'A', 'G', 'G', 'C', 'A', 'A', 'G', 'G', 'C', 'A', 'A', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'T', 'C', 'G', 'A', 'A', 'A', 'A', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'G', 'C'], ['G', 'A', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G', 'C', 'C', 'A', 'G', 'G', 'C', 'A', 'A', 'G', 'G', 'C', 'A', 'A', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'T', 'C', 'G', 'A', 'A', 'A', 'A', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'G', 'C']]
	# read_list = [['G', 'G', 'C', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'A', 'A', 'T', 'A', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'G', 'C', 'A', 'A'], ['G', 'G', 'C', 'A', 'G', 'A', 'A', 'A', 'G', 'A', 'C', 'G', 'G', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'A', 'A', 'T', 'A', 'A', 'C', 'A', 'C', 'A', 'C', 'C', 'G', 'A', 'C', 'G', 'G', 'C', 'A', 'G', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'G', 'C', 'A', 'A'], ['A', 'G', 'G', 'C', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'A', 'G', 'G', 'C', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'G', 'G', 'A', 'G', 'G', 'C', 'A', 'C', 'G', 'G', 'C', 'G', 'A', 'A', 'G', 'A', 'A', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'C'], ['A', 'G', 'G', 'C', 'G', 'G', 'A', 'A', 'C', 'A', 'C', 'A', 'G', 'G', 'C', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'A', 'C', 'C', 'G', 'G', 'A', 'G', 'G', 'C', 'A', 'C', 'G', 'G', 'C', 'G', 'A', 'A', 'G', 'A', 'A', 'C', 'G', 'A', 'G', 'A', 'C', 'G', 'C'], ['C', 'C', 'A', 'C', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'C', 'G', 'C', 'G', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'C', 'G', 'A', 'C', 'A', 'G', 'G', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'C', 'G'], ['C', 'C', 'A', 'C', 'C', 'A', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'A', 'A', 'A', 'G', 'C', 'A', 'C', 'A', 'A', 'C', 'C', 'G', 'C', 'G', 'C', 'C', 'A', 'G', 'A', 'A', 'C', 'C', 'G', 'A', 'C', 'A', 'G', 'G', 'A', 'C', 'C', 'A', 'A', 'A', 'C', 'C', 'G'], ['G', 'A', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G', 'C', 'C', 'A', 'G', 'G', 'C', 'A', 'A', 'G', 'G', 'C', 'A', 'A', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'T', 'C', 'G', 'A', 'A', 'A', 'A', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'G', 'C'], ['G', 'A', 'C', 'G', 'G', 'C', 'C', 'C', 'G', 'G', 'A', 'A', 'G', 'G', 'C', 'C', 'A', 'G', 'G', 'C', 'A', 'A', 'G', 'G', 'C', 'A', 'A', 'C', 'A', 'A', 'C', 'C', 'C', 'A', 'T', 'C', 'G', 'A', 'A', 'A', 'A', 'C', 'A', 'G', 'C', 'A', 'G', 'C', 'G', 'C']]

	# template_list.extend(read_list)

	# print len(template_list)
	# alleles = [[0, 0, 0, 0] for n in range(0,50)]

	# for r_index, r in enumerate(read_list):
	# 	if r_index % 2 == 0:
	# 		read = Read(r, template_list, [10, 5], 'CC@FFFFFHHGHGFH;EAEEIIIIFHIIFHGGIIIIHEIIGCGFIIHCAG', [alleles, alleles])
	# 	if r_index % 2 == 1:
	# 		read = Read(r, template_list, [5, 10], 'CC@FFFFFHHGHGFH;EAEEIIIIFHIIFHGGIIIIHEIIGCGFIIHCAG', [alleles,alleles])
	# 	reads.append(read)

	# Read in the data passed into the program
	#sample_read = Read("ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", ["ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATGTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATCCCCTTATACTCTCTATCTAAT", "ATATCGGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATTTCCGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCACCAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATTCCCTTATACTGTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTGAT"], [0.01, 0.01, 0.001, 0.1, 0.01, 0.1, 0.01, 0.1], 'CC@FFFFFHHGHGFH;EAEEIIIIFHIIFHGGIIIIHEIIGCGFIIHCAG')

	#sample_read = Read("ACTG", ["ATCG", "ACTG", "TCTG", "ACTA", "ACTA", "ACTG", "ATCG", "ATTG"], [0.1, 0.02, 0.1, 0.9, 0.5, 0.3, 0.5], "CC@F")

	#reads.append(sample_read)

	# For each read, generate the position assignment
	# for read in reads:
	# 	# TODO: Look at GMM and add in the bumps from the email
	# 	# read = find_position(read)
	# 	# read.print_read()
	# 	read_gmm(read)

	mult_reads_gmm(reads, training_reads, components)

	return 0

def mult_reads_gmm(reads, training_reads, components):
	"""
	Gaussian Mixture Model for multiple read
	"""

	prediction_zero_100 = 0
	prediction_one_100 = 0
	prediction_zero_200 = 0
	prediction_one_200 = 0

	base_opts = ['A', 'C', 'G', 'T']


	model = mixture.GMM(n_components=components, covariance_type='spherical')
	num_reads = len(reads)

	training_reads = [read.get_read().replace('\'', '') for read in training_reads]

	read_input = [read.get_read().replace('\'', '') for read in reads]
	# alignment_inputs = []
	# alignment_inputs.extend(read.get_alignments())

	# Generates observations
	# bases are converted to their ascii character values
	read_list = []
	for read in read_input:
		read_char = [convert_letter(c) for c in read]
		read_list.append(read_char)

	observations = []
	
	for alignment in training_reads:
		alignment_list = [convert_letter(c) for c in alignment] 
		observations.append( alignment_list )
	# for base_index, base in enumerate(read_main):
	# 	base_observations = [ord(base)]
	# 	for alignment in alignments:
	# 		base_observations.append(ord(alignment[base_index]))

	# 	observations.append(base_observations)

	model.fit(observations)
	means =  np.round(model.means_, 2)
	covars = np.round(model.covars_, 2)
	converted_means = []
	for num_list in means:
		# convert to nearest acceptable letter
		#char_means = [chr(int(n)) for n in num_list]
		char_means = [convert_to_letter(n) for n in num_list]
		converted_means.append(char_means)
	
	predictions = model.predict(read_list)

	read_predictions = []
	for index, prediction in enumerate(predictions):
		mapping = [prediction, reads[index]]
		read_predictions.append(mapping)
	

	for read_pr in read_predictions:
		
		prediction = read_pr[0]
		# def filt(x): return x[0] == prediction
		# matches = filter(filt, read_predictions)
		pr = prediction
		rps = int(float(read_pr[1].get_position()))
		# print '\n'
		# print prediction
		# print 'Converted Means: '
		# print ''.join(converted_means[prediction])
		# print 'Actual Read'
		# print read_pr[1].get_read()
		# print read_pr[1].get_position()
		# print 'Matches'
		# for m in matches:
		# 	print m[1].get_read() + ' Position: ' + m[1].get_position()
		# 	m[1].print_read()

		if pr == 0:
			if rps == 100:
				prediction_zero_100 = prediction_zero_100 + 1
				
			else:
				prediction_zero_200 = prediction_zero_200 + 1
				
		else:
			if rps == 100:
				prediction_one_100 = prediction_one_100 + 1
				
			else:
				prediction_one_200 = prediction_one_200 + 1
				

	print '\n-------------Predictions---------------------'
	print 'Prediction: 0 Position: 100 Num: ' + str(prediction_zero_100)
	print 'Prediction: 1 Position: 100 Num: ' + str(prediction_one_100)
	print 'Prediction: 0 Position: 200 Num: ' + str(prediction_zero_200)
	print 'Prediction: 1 Position: 200 Num: ' + str(prediction_one_200)

	print '\n------Means: -----------'
	for mean in converted_means:
		print ''.join(mean) 

	# for index, prediction in enumerate(predictions):
	# 	print 'Read: '
	# 	print reads[index].get_read()
	# 	print 'Prediction: '
	# 	print prediction
	# 	print converted_means[prediction]
	# 	print 'Means: '
	# 	print means[prediction]
	# 	print covars[prediction]
	# 	print '----------------------------------------\n'


	# posteriors = model.predict_proba(read_list)
	# print model.get_params(deep=True)
	# sample = model.sample()
	# print [convert_to_letter(n) for n in sample[0]]

 

def convert_to_letter(num):
		if num <= 50:
			return 'A'
		elif num <= 150:
			return 'C'
		elif num <= 250:
			return 'G'
		elif num <= 350:
			return 'T'
		else:
			return 'N'


def convert_letter(char):
	base_values = {
	'A' : 0,
	'C' : 100,
	'G' : 200,
	'T' : 300, 
	'N' : 400
	}
	# base_values = {
	# 'A' : 0,
	# 'C' : 0.25,
	# 'G' : 0.5,
	# 'T' : 0.75, 
	# '\'' : 1
	# }

	return base_values[char]


def read_gmm(read):
	"""
	Gaussian Mixture Model for the read
	"""
	model = mixture.GMM(n_components=8)
	read_main = read.get_read()
	alignments = read.get_alignments()

	# Generates observations
	# bases are converted to their ascii character values
	read_list = [ord(c) for c in read_main]
	observations = [ ]
	for alignment in alignments:
		alignment_list = [ord(c) for c in alignment] 
		observations.append( alignment_list )
	# for base_index, base in enumerate(read_main):
	# 	base_observations = [ord(base)]
	# 	for alignment in alignments:
	# 		base_observations.append(ord(alignment[base_index]))

	# 	observations.append(base_observations)

	print model.fit(observations)
	print np.round(model.means_, 2)
	
	print model.predict([read_list])




def find_position(read, alpha=0.01, beta=0.01, k=0.6):
	"""
	This finds the optimal position for the read based on a GMM that generates probable templates based on the
	the read. From the generated reads, it finds the one with the
	highest score and assigns that as the optimal position
	"""
	# VARIABLE DECLARATION
	# alpha = insertion probability
	# beta = deletion probability
	# k = exploratory
	pi = 0 # index of the base of interest
	not_alpha = 1 - alpha
	not_beta = 1 - beta
	read_main = read.get_read()

	# the possibilites for a base, the four nucleotide bases, N for not enough information, R for
	# the removal of a base +_ for the insertion of a base
	# TODO: Get rid of N, fix it, check the init_bases loop below and change
	base_opts = ['A', 'C', 'G', 'T', 'R', '+A', '+C', '+G', '+T']

	# TODO: Move this so only calculate once
	# initialize with uniform probability for each option
	init_base_opts_prob = [0.1 for b in base_opts]

	# multiply by alpha or ~alpha for insertion and beta or ~beta for deletion
	for index, base_prob in enumerate(init_base_opts_prob):
		# the first four bases are the regular nucleotides with the 5th being a base without info
		if index <= 3:
			init_base_opts_prob[index] = base_prob *  not_alpha * not_beta
		# the 6th base option is a deletion
		if index is 4:
			init_base_opts_prob[index] = base_prob * not_alpha * beta
		# the final four bases are insertions
		if index > 4:
			init_base_opts_prob[index] = base_prob * alpha * not_beta

	# for read, generate possible alignments
	generated_templates = []

	for base_index, base in enumerate(read_main):
		# calculate prior based on the quality score of the base and the quality score of the alignment
		# there is a higher probability of the base in the read if the quality score of the base is high

		# TODO: Get the alignment distribution for each base index, higher probability if the base is one of
		# the previously found
		prior = read.get_base_quality_score(base_index) * k
		not_prior = 1 - prior
		base_choice = base_opts.index(base)
		base_opts_prob = init_base_opts_prob[:]
		base_prob_list = read.get_base_probs(base_index)

		for index, base_prob in enumerate(base_opts_prob):
			if index is base_choice:
				base_opts_prob[index] = base_prob * prior * base_prob_list[index] 
			else:
				base_opts_prob[index] = base_prob * not_prior * base_prob_list[index]
		

		# Normalize base_opts_prob
		prob_sum = sum(base_opts_prob)
		final_base_opts_prob = [prob/prob_sum for prob in base_opts_prob]
		# Generate Draws
		xk = arange(9)
		discrete = stats.rv_discrete(name="discrete",values=(xk, final_base_opts_prob))
		draws = discrete.rvs(size=10000)

		# If it is the first base, add a new list with the base draw
		if not generated_templates:
			for base_draw in draws:
				generated_templates.append([base_opts[base_draw]])

		else:
			# Append draw to generated reads
			for read_index, base_draw in enumerate(draws):
				generated_templates[read_index].append(base_opts[base_draw])

	# Find the matching generated templates and get score
	scores = find_generated_templates(read.get_alignments(), generated_templates)

	# Set score
	for alignment_index, score in enumerate(scores):
		read.set_alignment_probability(alignment_index, score)

	# Set the position for the read to the highest position
	read.find_highest_position()

	# After generating a score for each alignment, return the read
	return read
			

def find_generated_templates(alignments, generated_templates):
	"""
	Finds the matching generated templates and returns that probability score for the associated template
	It will return 0 if not found
	"""
	# initialize to a score of zero
	template_scores = [0 for a in alignments]
	template_strings = [''.join(t) for t in generated_templates]
	num_generated = len(generated_templates)
	
	for alignment_index, alignment in enumerate(alignments):
		count = template_strings.count(alignment)
		template_scores[alignment_index] = count / num_generated

	# Think about normalizing all the scores for the alignments
	return template_scores
def convert_to_phred(score):
	"""
	Converts Fastq phred scores into values
	"""
	# TODO: might need to change based on illumina
	score_list = []
	for value in score:
		numeral = ord(value) - 33
		pe_score = pow(10, (numeral/-10))
		score_list.append(1 - pe_score)
	
	return score_list

def read_in_file():
	"""
	Reads in the files passed in to the script
	"""
	# Declare variables
	reads = []

	# Get command line arguments
	arguments = sys.argv
	arguments_length = len(arguments)

	# Read file is the first argument
	read_file_name = arguments[1]

	# Process read file 
	read_file = open(read_file_name, 'r')
	for line in read_file:
		read_info = line.split()
		read_string = read_info[2].replace('\'', '')
		new_read = GenerativeRead(read_string, [], read_info[5], read_info[3], None, [], read_info[0], read_info[1], read_info[4]) 
		reads.append(new_read)
	read_file.close()

	# Repeat regions file in the second argument
	repeat_file_name = arguments[2]

	# Process repeat file
	repeat_file = open(repeat_file_name, 'r')
	alignments = [[]]
	alignment_index = -1
	previous_line = ''


	for line in repeat_file:
		alignment_info = line.split()

		# This consists of a tuple of alignment string, alignment start position and alignment chromosome
		#new_align = alignment_info[2], alignment_info[4], alignment_info[3]

		new_align = Alignment(alignment_info[2], None, alignment_info[4], alignment_info[3])

		if previous_line != alignment_info[0]:
			# It is not a repeat
			alignment_index = alignment_index + 1
			alignments.append([])
			previous_line = alignment_info[0]

		alignments[alignment_index].append(new_align)

	repeat_file.close()

	# Associate each read with the other alignments
	for read in reads:
		# Find the other alignments
		pos = read.get_position()
		found = False
		found_index = -1

		for a_index, alignment_lists in enumerate(alignments):
			# find matching alignments
			# TODO: Don't add alignment already have
			# TODO: Make functional with filter
			for align in alignment_lists:
				if align.get_position() == pos:
					found = True
					found_index = a_index
					break

			if found is True:
				break

		if found is True:
			for new_align in alignments[found_index]:
				read.add_alignment(new_align)
			


	# SNP files are the remaining ones
	snp_file_names = [arguments[file_id] for file_id in range(3, arguments_length) ]

	# Process SNP files
	for file_name in snp_file_names:
		snp_file = open(file_name, 'r')

		for line in snp_file:
			snp_info = line.split()
			snps = snp_info[3].split('/')
			snp_pos = int(float(snp_info[2]))

			# Ignore alleles that are longer than one base

			
			if all(len(x) < 2 for x in snps):

				# Iterate through reads and determine whether or not it contains this SNP
				pos_low = snp_pos - 49
			

				for read in reads:
					positions = read.get_alignment_positions()

					for p_index, p in enumerate(positions):
						p = int(float(p))
						if p >= pos_low and p <= snp_pos:
							# Get index of snp
							offset = snp_pos - p
							calls = [0, 0, 0, 0]
							for snp in snps:
								call_index = get_base_num(snp)
								calls[call_index] = 1

							# Add the SNP to the read
							read.add_snp(p_index, offset, calls)
							
		snp_file.close()
	return reads

def filter_reads(reads, positions):
	"""
	Returns only the reads that have an overlap of over half with the given positions
	"""

	filtered_reads = []
	read_array = array(reads)

	for position in positions:
		low_position = position - 25
		high_position = position + 25
		p = []

		for r in reads:
			pos = int(float(r.get_position()))
			if pos > low_position and pos < high_position:
				p.append(r) 

		filtered_reads.extend(p)

	return filtered_reads

def combine_reads(filtered_reads, positions):
	"""
	Combines reads so that they all have the desired start position
	"""

	combined_reads = []
	true_reads = []

	for r in filtered_reads:
		# Find associated position
		r_position = float(r.get_position())
		desired_start = -1

		for p in positions:
			low_position = p - 25
			high_position = p + 25
			if r_position > low_position and r_position < high_position:
				desired_start = p
				break

		if desired_start is not -1:
			# Find another read that overlaps
			if r_position < desired_start:
				offset = desired_start - r_position
				for r2 in filtered_reads:
					r2_position = float(r2.get_position())
					if r2_position > desired_start and r2_position <= desired_start + offset and r2_position != r_position:
						fuse_read = r2
						break
			elif r_position == desired_start:
				fuse_read = None
			else:
				offset = r_position - desired_start
				for r2 in filtered_reads:
					r2_position = float(r2.get_position())
					r2_end = r2_position + 49
					if r2_end + 49 > desired_start and r2_end >= r_position - 1 and r2_position != r_position:
						fuse_read = r2 
						break

			if fuse_read is None:
				if r_position == desired_start:
					true_reads.append(r)
			else:
				r.fuse_read(fuse_read, desired_start)
				combined_reads.append(r)
				

	def f(x): return len(x.get_read()) == 50
	combined_reads = filter(f, combined_reads)

	# for c in combined_reads:
	# 	c.print_read()
	# 	print '\n'
	
	return (combined_reads, true_reads)

def create_training_set(reads, positions, snps, num_reads):
	training_reads = []
	num_given = len(reads)
	base_opts = ['A', 'C', 'G', 'T']

	read_snps = []
	for r_ind, r in enumerate(reads):
		alignment_snps = snps[r_ind]

		# Determine where the SNPs are
		known_snps = []

		for index, base in enumerate(alignment_snps):
			if base[0] != -1:
				known_snps.append([index, base])

		read_snps.append(known_snps)


	for t in range(0, num_reads):
		# Randomly choose a read
		r = randint(0,num_given-1)

		# Get alignment
		test_alignment = reads[r]
		test_position = positions[r]
		generated_template = list(test_alignment)
		align_snp = read_snps[r]
		snp_length = len(align_snp)

		# Randomly choose known snps to implement
		if snp_length > 1:
			num_snps = randint(0, snp_length-1)

			for s in range(num_snps):
				# choose index to mutate
				ind = randint(0, snp_length-1)
				base_ind = randint(0,3)

				# get which index is non zero
				call_ind = align_snp[ind][0]
				base_snp = align_snp[ind][1]

				while base_snp[base_ind] is 0:
					base_ind = (base_ind + 1) % 4
				

				generated_template[call_ind] = base_opts[base_ind] 


		# Choose number of mutation
		num_mutations = randint(0,5)

		for mutation in range(0,num_mutations):
			# choose index to mutate
			ind = randint(0,49)
			base_ind = randint(0,3)
			generated_template[ind] = base_opts[base_ind]

		generated_string = ''.join(generated_template) 

		# Create New Read
		test_read = GenerativeRead(generated_string, [test_alignment], [], [], [], [test_position], [], test_position)
		training_reads.append(test_read)


	return training_reads 

if __name__ == '__main__':
	reads = read_in_file()

	#reads = filter_reads(reads, [49937899, 22382490, 100790155])
	reads = filter_reads(reads, [100, 200])
	
	#returned_combined = combine_reads(reads, [49937899, 22382490, 100790155])
	returned_combined = combine_reads(reads, [100, 200])
	reads = returned_combined[0]
	true_reads = returned_combined[1]

	templates = true_reads[0].get_alignments()
	template_positions = true_reads[0].get_alignment_positions()
	template_snps = true_reads[0].get_all_alignment_alleles()

	print '------------Templates----------------'
	print templates

	for test in range(0, 10):	
		# Create a new set of training reads
		training_reads = create_training_set(templates, template_positions, template_snps, 1000)

		# Run the program
		ambio(true_reads, training_reads, 2)

