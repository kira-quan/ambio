#!/usr/bin/env python
# This program takes in a set of reads and the template strand. For each read, where it aligns on the 
# placement strand, sequences of likely reads are produced from the template using a Gaussian Mixture
# Model. The read is then assigned to the position on the template where it is most likely

from __future__ import division
from scipy import stats
from read import Read, Alignment
from numpy import arange
from sklearn import mixture # Version 0.14.1
import timeit
import numpy as np
import sys


def ambio():
	# VARIABLE DECLARATION
	reads = [] # a list of the reads from the subject

	# Read in the data passed into the program
	#sample_read = Read("ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", ["ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATGTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATCCCCTTATACTCTCTATCTAAT", "ATATCGGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATTTCCGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCACCAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATTCCCTTATACTGTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTGAT"], [0.01, 0.01, 0.001, 0.1, 0.01, 0.1, 0.01, 0.1], 'CC@FFFFFHHGHGFH;EAEEIIIIFHIIFHGGIIIIHEIIGCGFIIHCAG')

	alleles = [0.05, 0, 0, 0.95], [0, 0.5, 0, 0.5], [0, 0, 0, 1], [0.1, 0, 0.9, 0]
	changed_alleles = [1, 0, 0, 0], [0, 0.5, 0, 0.5], [0, 0, 0, 1], [0.5, 0, 0.5, 0]
	sample_read = Read("ACTG", ["ATCG", "ACAG", "TCTG", "ACTA", "ACTA", "ACTA", "ATCG", "ATTG"], [0.1, 0.02, 0.1, 0.9, 0.5, 0.3, 0.5], "CC@F", [alleles, alleles, alleles, changed_alleles, alleles, alleles, alleles, alleles])

	reads.append(sample_read)

	# For each read, generate the position assignment
	for read in reads:
		# TODO: Look at GMM and add in the bumps from the email
		if len(read.get_alignments()) > 1:
			read = find_position(read)
			read.print_read()
		else:
			# There is only one alignment possibility
			read.set_alignment_probability(0, 1)
			read.print_read()
		
	return 0

def find_position(read,t=1,a=1,p=1):
	"""
	This finds the optimal position for the read based on a GMM that generates probable templates based on the
	the read. From the generated reads, it finds the one with the
	highest score and assigns that as the optimal position
	t, a, and p are all weighting scores for transitions, allele frequencies, and phred scores respectively
	"""
	# VARIABLE DECLARATION
	read_main = read.get_read()

	# the possibilites for a base, the four nucleotide bases, N for not enough information, R for
	# the removal of a base +_ for the insertion of a base
	base_opts = ['A', 'C', 'G', 'T']

	# Transition probabilities A <-> G and C <-> T are twice as likely as any other conversion
	# -1 denotes the index of the given base key
	transitions = {
		'A' : [-1, 0.5, 1, 0.5],
		'C' : [0.5, -1, 0.5, 1],
		'G' : [1, 0.5, -1, 0.5],
		'T' : [0.5, 1, 0.5, -1],
		'N'	: [1, 1, 1, 1]
	}

	# initialize with uniform probability for each option
	init_base_opts_prob = [0.1 for b in base_opts]

	for base_index, base in enumerate(read_main):
		# calculate prior based on the quality score of the base and the quality score of the alignment
		# there is a higher probability of the base in the read if the quality score of the base is high

		# Calculate likelihood
		# Phred quality scores
		phred_score = read.get_base_quality_score(base_index)

		# Allele frequencies
		# Moved to the read class
		# base_allele_freq = read.get_allele_frequencies(base_index)

		likelihood = [0.25 for r in range(0,4)]
		base_transition = transitions[base]

		# keep likelihoods uniform if N
		if base is not 'N':
			base_choice = base_opts.index(base)
			for index, l_prob in enumerate(likelihood):
				if base_choice is index:
					likelihood[index] = l_prob * (phred_score * p)
				else:
					# Transition versus Transversion added here
					likelihood[index] = l_prob * ((1 - phred_score) * p) * (base_transition[index] * t)
			
		base_opts_prob = init_base_opts_prob[:]
		
		for index, base_prob in enumerate(base_opts_prob):
			base_opts_prob[index] = base_prob * likelihood[index]

		# Normalize base_opts_prob
		prob_sum = sum(base_opts_prob)
		final_base_opts_prob = [np.nextafter(prob/prob_sum, 1) for prob in base_opts_prob]

		# Update alignment scores
		read.update_alignment_probability_with_list(base_index, final_base_opts_prob, a)
		
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
		new_read = Read(read_string, [], read_info[5], read_info[3], None, [], read_info[0], read_info[1], read_info[4]) 
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
			print '\n'
			print len(alignments[found_index])
			for new_align in alignments[found_index]:
				read.add_alignment(new_align)
				print new_align.get_position()
			
			read.print_read()


	# SNP files are the remaining ones
	snp_file_names = [arguments[file_id] for file_id in range(3, arguments_length) ]

	# Process SNP files
	for file_name in snp_file_names:
		snp_file = open(file_name, 'r')
		snp_file.close()

if __name__ == '__main__':
	read_in_file()
	#print timeit.timeit(ambio, number=1)







