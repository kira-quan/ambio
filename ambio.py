#!/usr/bin/env python
# This program takes in a set of reads and the template strand. For each read, where it aligns on the 
# placement strand, sequences of likely reads are produced from the template using a Gaussian Mixture
# Model. The read is then assigned to the position on the template where it is most likely

from __future__ import division
from scipy import stats
from read import Read

def ambio():
	# VARIABLE DECLARATION
	reads = [] # a list of the reads from the subject

	# Read in the data passed into the program

	# For each read, generate the position assignment
	for read in reads:
		read = find_position(read)

	return 0

def find_position(read, alpha=0, beta=0):
	"""
	This finds the optimal position for the read based on a GMM that generates probable reads based on the
	template strands that this read could be assigned to. From the generated reads, it finds the one with the
	highest score and assigns that as the optimal position
	"""
	# VARIABLE DECLARATION
	# alpha = insertion probability
	# beta = deletion probability
	pi = 0 # index of the base of interest
	not_alpha = 1 - alpha
	not_beta = 1 - beta

	# the possibilites for a base, the four nucleotide bases, N for not enough information, R for
	# the removal of a base +_ for the insertion of a base
	base_opts = ['A', 'C', 'G', 'T', 'N', 'R', '+A', '+C', '+G', '+T']

	# TODO: Move this so only calculate once
	# initialize with uniform probability for each option
	init_bases_opts_prob = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 0.1]

	# multiply by alpha or ~alpha for insertion and beta or ~beta for deletion
	for base_prob, index in enumerate(init_bases_opts_prob):
		# the first four bases are the regular nucleotides with the 5th being a base without info
		if index <= 4:
			base_prob = base_prob *  not_alpha * not_beta
		# the 6th base option is a deletion
		if index is 5:
			base_prob = base_prob * not_alpha * beta
		# the final four bases are insertions
		if index > 5:
			base_prob = base_prob * alpha * not_beta


	# for each alignment, generate read draws
	for alignment, alignment_index in enumerate(read.alignments):
		alignment_quality_score = read.get_alignment_quality_score(alignment_index)
		generated_reads = []

		# intialize with initial base probabilties
		base_opts_prob = init_bases_opts_prob

		for base, base_index in enumerate(alignment):
			# calculate prior based on the quality score of the base and the quality score of the alignment
			# there is a higher probability of the base in the read if the quality score of the base is high
			# and the alignment quality score is high
			prior = read.get_bases_quality_score(base_index) * alignment_quality_score
			not_prior = 1 - prior
			base_choice = base_opts.find(base)

			for base_prob, index in base_opts_prob:
				if index is base_choice:
					base_prob = base_prob * prior
				else:
					base_prob = base_prob * not_prior

			# Normalize base_opts_prob
			prob_sum = sum(base_opts_prob)
			final_base_opts_prob = [prob/prob_sum for prob in base_opts_prob]

			# Generate Draws
			draws = stats.rv_discrete(0, 9, values=(range(10), final_base_opts_prob)).rvs(size=100)

			# If it is the first base, add a new list with the base draw
			if generated_reads is empty:
				for base_draw in draws:
					generated_reads.append([base_opts[base_draw]])

			else:
				# Append draw to generated reads
				for base_draw, read_index in enumerate(draws):
					generated_reads[read_index].append(base_opts[base_draw])

			# Find the matching generated read and get score
			score = find_generated_read(read.get_read(), generated_reads)

			# Set score
			read.set_alignment_probability(alignment_index, score)

	# Set the position for the read to the highest position
	read.find_highest_position()

	# After generating a score for each alignment, return the read
	return read
			

def find_generated_read(read, generated_reads):
	"""
	Finds the matching generated read and returns that probability score
	It will return -1 if not found
	"""
	# convert to set for comparison purpose
	read_set = set(read)
	generated_read_score = -1

	for generated_read  in generated_reads:
		if read_set == set(generated_read):
			# T0D0: Figure out how to calculate probability score of each generated read
			# This is a placeholder for now
			generated_read_score = 0
			break

	return generated_read_score







