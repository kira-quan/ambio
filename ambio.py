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
	sample_read = Read("ACGT", ["ACGT", "CCGT"], [0.01, 0.01], [0.1, 0.001, 0.0001, 0.01])
	reads.append(sample_read)

	# For each read, generate the position assignment
	for read in reads:
		# TODO: Look at GMM and add in the bumps from the email
		read = find_position(read)
		print 

	return 0

def find_position(read, alpha=0.01, beta=0.01):
	"""
	This finds the optimal position for the read based on a GMM that generates probable templates based on the
	the read. From the generated reads, it finds the one with the
	highest score and assigns that as the optimal position
	"""
	# VARIABLE DECLARATION
	# alpha = insertion probability
	# beta = deletion probability
	pi = 0 # index of the base of interest
	not_alpha = 1 - alpha
	not_beta = 1 - beta
	read_main = read.get_read()

	# the possibilites for a base, the four nucleotide bases, N for not enough information, R for
	# the removal of a base +_ for the insertion of a base
	# TODO: Get rid of N, fix it, check the init_bases loop below and change
	base_opts = ['A', 'C', 'G', 'T', 'N', 'R', '+A', '+C', '+G', '+T']

	# TODO: Move this so only calculate once
	# initialize with uniform probability for each option
	base_opts_prob = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 0.1]

	# multiply by alpha or ~alpha for insertion and beta or ~beta for deletion
	for index, base_prob in enumerate(base_opts_prob):
		# the first four bases are the regular nucleotides with the 5th being a base without info
		if index <= 4:
			base_prob = base_prob *  not_alpha * not_beta
		# the 6th base option is a deletion
		if index is 5:
			base_prob = base_prob * not_alpha * beta
		# the final four bases are insertions
		if index > 5:
			base_prob = base_prob * alpha * not_beta


	# for read, generate possible alignments		
	generated_templates = []

	for base_index, base in enumerate(read_main):
		# calculate prior based on the quality score of the base and the quality score of the alignment
		# there is a higher probability of the base in the read if the quality score of the base is high
		# and the alignment quality score is high
		prior = read.get_bases_quality_score(base_index)
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
		draws = stats.rv_discrete(0, 9, values=(range(10), final_base_opts_prob)).rvs(size=10)

		# If it is the first base, add a new list with the base draw
		if generated_templates is empty:
			for base_draw in draws:
				generated_reads.append([base_opts[base_draw]])

		else:
			# Append draw to generated reads
			for read_index, base_draw in enumerate(draws):
				generated_templates[read_index].append(base_opts[base_draw])

		# Find the matching generated templates and get score
		scores = find_generated_templates(read.get_alignments(), generated_templates)
		print "Generated Templates: "
		print generated_templates

		# Set score
		for alignment_index, score in scores:
			read.set_alignment_probability(alignment_index, score)

	# Set the position for the read to the highest position
	read.find_highest_position()

	# After generating a score for each alignment, return the read
	return read
			

def find_generated_read(alignments, generated_templates):
	"""
	Finds the matching generated templates and returns that probability score for the associated template
	It will return 0 if not found
	"""
	# initialize to a score of zero
	template_scores = [0 for a in alignments]
	num_generated = len(generated_templates)
	
	for alignment_index, alignment in enumerate(alignments):
		count = generated_templates.count(alignment)
		template_scores[alignment_index] = count / num_generated

	return template_scores







