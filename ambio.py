#!/usr/bin/env python
# This program takes in a set of reads and the template strand. For each read, where it aligns on the 
# placement strand, sequences of likely reads are produced from the template using a Gaussian Mixture
# Model. The read is then assigned to the position on the template where it is most likely

from __future__ import division
from scipy import stats
from read import Read
from numpy import arange
import timeit

def ambio():
	# VARIABLE DECLARATION
	reads = [] # a list of the reads from the subject

	# Read in the data passed into the program
	sample_read = Read("ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", ["ATATCCCTACCAATCTATCCCCAAAAATTCCCTTATACTCTCTATCTAAT", "ATATCCGTACCAATCTATCCCCAACAATTCCCTTATACTCTCTATCTAAT"], [0.01, 0.01], [0.12117885924899019, 0.33415160102483166, 0.10632170228284443, 0.14686061309953347, 0.36473687807903743, 0.01788901020746836, 0.4556917988559128, 0.05398658415199342, 0.3141159989393558, 0.5188456039136414, 0.09759694223390225, 0.6238380031950486, 0.3518641378282821, 0.044791323688791684, 0.9623702159874176, 0.010160983349749797, 0.8957657921549259, 0.6371182123240026, 0.6692419521899688, 0.9793392110614689, 0.8495260682941772, 0.8060781219027288, 0.7690154226048116, 0.21056506302357447, 0.1139836759312679, 0.9014064122332206, 0.19196431719394835, 0.28746715704141923, 0.9962959228799395, 0.8927351709880836, 0.3583972951261465, 0.7808689587097891, 0.9727105799239993, 0.2927723520315926, 0.4530579140475315, 0.3647337645210673, 0.28903500743577293, 0.8009388450973709, 0.33743469460700337, 0.7943477102002011, 0.0015752540267214288, 0.7880426755460618, 0.27760313914221435, 0.33296403513830775, 0.09086776155049259, 0.41120800537971347, 0.9670421855602295, 0.3672821866743671, 0.23757120448243074, 0.6763406432571134]
)
	reads.append(sample_read)

	# For each read, generate the position assignment
	for read in reads:
		# TODO: Look at GMM and add in the bumps from the email
		read = find_position(read)
		#read.print_read()

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
		# and the alignment quality score is high
		prior = read.get_base_quality_score(base_index)
		not_prior = 1 - prior
		base_choice = base_opts.index(base)
		base_opts_prob = init_base_opts_prob[:]

		for index, base_prob in enumerate(base_opts_prob):
			if index is base_choice:
				base_opts_prob[index] = base_prob * prior
			else:
				base_opts_prob[index] = base_prob * not_prior
		

		# Normalize base_opts_prob
		prob_sum = sum(base_opts_prob)
		final_base_opts_prob = [prob/prob_sum for prob in base_opts_prob]
		# Generate Draws
		xk = arange(9)
		discrete = stats.rv_discrete(name="discrete",values=(xk, final_base_opts_prob))
		draws = discrete.rvs(size=1000)

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

if __name__ == '__main__':
	print timeit.timeit(ambio, number=10000)





