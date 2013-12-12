#!/usr/bin/env python
# This program takes in a set of reads and the template strand. For each read, where it aligns on the 
# placement strand, sequences of likely reads are produced from the template using a Gaussian Mixture
# Model. The read is then assigned to the position on the template where it is most likely

from __future__ import division
from scipy import stats
from read import Read
from numpy import arange
from random import randint
from sklearn import mixture # Version 0.14.1
import timeit
import numpy as np
import sys


def ambio(reads):
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

	mult_reads_gmm(reads)

	return 0

def mult_reads_gmm(reads):
	"""
	Gaussian Mixture Model for multiple read
	"""
	base_opts = ['A', 'C', 'G', 'T']


	model = mixture.GMM(n_components=25)
	num_reads = len(reads)

	# Train on 10% of the reads
	train_num = int(num_reads * 0.5)
	
	training_reads = []

	for train in range(0, train_num):
		ind = randint(0, (num_reads-1))
		training_reads.append(reads[ind].get_read())

	read_input = [read.get_read() for read in reads]
	# alignment_inputs = []
	# alignment_inputs.extend(read.get_alignments())

	# Generates observations
	# bases are converted to their ascii character values
	read_list = []
	for read in read_input:
		read_char = [convert_letter(c) for c in read]
		read_list.append(read_char)

	observations = []
	print training_reads
	for alignment in training_reads:
		alignment_list = [convert_letter(c) for c in alignment] 
		observations.append( alignment_list )
	# for base_index, base in enumerate(read_main):
	# 	base_observations = [ord(base)]
	# 	for alignment in alignments:
	# 		base_observations.append(ord(alignment[base_index]))

	# 	observations.append(base_observations)

	print model.fit(observations)
	means =  np.round(model.means_, 2)
	converted_means = []
	for num_list in means:
		# convert to nearest acceptable letter
		char_means = []
		#char_means = [chr(int(n)) for n in num_list]

		for n in num_list:
			if n < 0:
				char_means.append('A')
			elif n < 150:
				char_means.append('C')
			elif n < 250:
				char_means.append('G')
			elif n < 350:
				char_means.append('T')
			else:
				char_means.append('\'')


		converted_means.append(char_means)
	
	predictions = model.predict(read_list)
	
	for index, prediction in enumerate(predictions):
		print 'Read: '
		print reads[index].get_read()
		print 'Prediction: '
		print prediction
		print converted_means[prediction]
		print 'Means: '
		print means[prediction]
		print '----------------------------------------\n'

def convert_letter(char):
	base_values = {
	'A' : 0,
	'C' : 100,
	'G' : 200,
	'T' : 300, 
	'\'' : 400
	}

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
		new_read = Read(read_info[2], [read_info[2]], read_info[5], read_info[3], None, [read_info[1]], read_info[0], read_info[4] ) 
		reads.append(new_read)
	read_file.close()

	return reads

if __name__ == '__main__':
	reads = read_in_file()
	ambio(reads)







