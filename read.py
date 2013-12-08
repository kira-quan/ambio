#!/usr/bin/env python
# This file defines the Read class for a read from a DNA alignment. This class contains the read,
# associated alignments, scores for those alignments and bases, and the position to align the 
# read

from __future__ import division
from numpy import nextafter

class Read:
	def __init__(self, read, alignments, alignment_quality_scores, bases_quality_score, allele_frequencies):
		# the read
		self.read = read
		# a list of the possible alignments for a read
		# the first alignment is the original alignment
		self.alignments = alignments 
		# a list of the corresponding quality scores of each alignment, same order as alginments
		self.alignment_quality_scores = alignment_quality_scores
		# a list of the quality score for each base (in order) in a read
		self.bases_quality_score = self.convert_to_phred(bases_quality_score)
		# the best position for the read, where the position refers to the alignment index
		self.position = 0
		# a list of generate probability score for each alignment
		self.alignment_probability_scores = [-1 for a in alignments]
		# alternative base calls
		self.base_call_possibilities = self.find_alternative_alignments()
		# a list of allele frequencies for each base where each entry in the list is a
		# a list of frequencies in the following order [A, C, G, T]
		self.allele_frequencies = [[]for a in alignments]
		for align_index, align  in enumerate(allele_frequencies):
			for base_index, base in enumerate(align):
				self.allele_frequencies[align_index].append(base)
				# for call in base:
				# 	if call is 0:
				# 		# Don't want probabilities to be zero for other options
				# 		self.allele_frequencies[align_index][base_index] = nextafter(call, 1)
		print self.allele_frequencies[0]
		# bases
		self.base_opts = ['A', 'C', 'G', 'T']

	def get_allele_frequencies(self, base_index):
		"""
		Return the entry in the allele frequencies list at the given base index
		"""
		return self.allele_frequencies[base_index]

	def find_alternative_alignments(self):
		"""
		Generates other possibilities for a base call based on the given alignments
		"""
		base_opts = ['A', 'C', 'G', 'T', 'R', '+A', '+C', '+G', '+T']
		call_possibilities = []
		num_alignments = len(self.alignments)
		for index, base_call in enumerate(self.read):
			prob_base = [0 for o in base_opts]
			# find other alignments
			for a in self.alignments:
				call = a[index]
				opt_index = base_opts.index(call)
				prob_base[opt_index] = prob_base[opt_index] + 1
				
			prob_base = [base / num_alignments for base in prob_base]
			call_possibilities.append(prob_base)
		return call_possibilities

	def print_read(self):
		print "Read is: " + self.read
		print "Alignments are: "
		for alignment_index, alignment in enumerate(self.alignments):
			print alignment
			print "Score: " + str(self.alignment_probability_scores[alignment_index])
		print "Position is: " + str(self.position)

	def get_base_probs(self, index):
		return self.base_call_possibilities[index]

	def get_alignments(self):
		return self.alignments

	def get_alignment_quality_score(self, alignment_index):
		return self.alignment_quality_scores[alignment_index]

	def get_base_quality_score(self, base_index):
		return self.bases_quality_score[base_index]

	def get_position(self):
		return self. position

	def set_position(self, new_position):
		self.position = new_position

	def set_alignment_probability(self, alignment_index, alignment_prob):
		# TODO: check that doesn't set all of them
		self.alignment_probability_scores[alignment_index] = alignment_prob

	def update_alignment_probability(self, alignment_index, new_prob):
		"""
		If the current probability is -1 that means it hasn't been set yet and
		it will set the alignment probability to new_prob, if not then it will multiply
		the alignment probability by new_prob
		"""
		current_score = self.alignment_probability_scores[alignment_index]
		if current_score is -1:
			self.alignment_probability_scores[alignment_index] = new_prob
		else:
			self.alignment_probability_scores[alignment_index] = current_score * new_prob 

	def update_alignment_probability_with_list(self, base_index, likelihoods):
		"""
		For the given index and likelihood list, updates each alignment with the probability
		"""
		base_opts = self.base_opts
		for align_index, alignment in enumerate(self.alignments):
			# Find probability associated with the base call in that alignment
			
			alignment_base_index = base_opts.index(alignment[base_index])
			if align_index is 1:
				print alignment_base_index
				print alignment[base_index]

			# Include alignment allele frequency
			print align_index
			print self.allele_frequencies[align_index]
			base_prob = likelihoods[alignment_base_index] * self.allele_frequencies[align_index][base_index][alignment_base_index]
			self.update_alignment_probability(align_index, base_prob)

		print 'updated alignment scores'
		print self.alignment_probability_scores

	def get_read(self):
		return self.read

	def find_highest_position(self):
		"""
		Sets the position to the alignment with the highest score
		"""
		highest_index = 0
		highest_score = -1

		# Normalize base_opts_prob
		prob_sum = sum(self.alignment_probability_scores)
		print prob_sum
		self.alignment_probability_scores = [(prob/prob_sum) for prob in self.alignment_probability_scores]
		print self.alignment_probability_scores

		for index, score in enumerate(self.alignment_probability_scores):
			if score > highest_score:
				highest_score = score
				highest_index = index
			if score is highest_score:
				# in case of tie, then assign to alignment with highest alignment score
				if self.alignment_quality_scores[index] > self.alignment_quality_scores[highest_index]:
					highest_score = score
					highest_index = index

		self.set_position(highest_index)
	
	def convert_to_phred(self, score):
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