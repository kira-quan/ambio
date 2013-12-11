#!/usr/bin/env python
# This file defines the Read class for a read from a DNA alignment. This class contains the read,
# associated alignments, scores for those alignments and bases, and the position to align the 
# read

from __future__ import division
from numpy import nextafter

# class Alignment:
# 	def __init__(self, alignment_sequence, alignment_quality_score, allele_frequencies, position):
# 		self.alignment = alignment_sequence
# 		self.quality_score = alignment_quality_score
# 		self.probability_score = -1
# 		self.allele_frequencies = []
# 		self.position = position
# 		self.allele_frequencies.append(base)
# 		for call_index, call in enumerate(base):
# 			if call is 0:
# 				# Don't want probabilities to be zero for other options
# 				self.allele_frequencies[base_index][call_index] = nextafter(call, 1)

# 	def get_allele_frequencies(self, base_index):
# 		"""
# 		Return the entry in the allele frequencies list at the given base index
# 		"""
# 		return self.allele_frequencies[base_index]


class Read:
	def __init__(self, read, alignments, mapq, bases_quality_score, allele_frequencies, alignment_positions=None, original_chromosome=None, cigar='50M'):
		# the read
		self.read = read
		# a list of the possible alignments for a read
		# the first alignment is the original alignment
		self.alignments = alignments 
		# a list of the corresponding quality scores of each alignment, same order as alginments
		#self.alignment_quality_scores = alignment_quality_scores
		# a list of the quality score for each base (in order) in a read

		self.mapq = mapq
		self.cigar = cigar

		self.bases_quality_score = self.convert_to_phred(bases_quality_score)
		# the best position for the read, where the position refers to the alignment index
		self.position = 0
		# a list of generate probability score for each alignment
		self.alignment_probability_scores = [-1 for a in alignments]
		# alternative base calls
		# self.base_call_possibilities = self.find_alternative_alignments()


		self.alignment_positions = alignment_positions
		self.original_chromosome = original_chromosome

		# a list of allele frequencies for each base where each entry in the list is a
		# a list of frequencies in the following order [A, C, G, T]
		if allele_frequencies is not None:
			self.allele_frequencies = [[]for a in alignments]
			for align_index, align  in enumerate(allele_frequencies):
				for base_index, base in enumerate(align):
					self.allele_frequencies[align_index].append(base)
					for call_index, call in enumerate(base):
						if call is 0:
							# Don't want probabilities to be zero for other options
							self.allele_frequencies[align_index][base_index][call_index] = nextafter(call, 1)
		else:
			self.allele_frequencies = None
		# bases
		self.base_opts = ['A', 'C', 'G', 'T']

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

	def print_all_read_info(self):
		print '\n----------------------------------------------\n'
		print "Read is: " + self.read
		print "Alignments are: "
		for alignment_index, alignment in enumerate(self.alignments):
			print alignment
			print "Score: " + str(self.alignment_probability_scores[alignment_index])
			print 'Position: ' + str(self.alignment_positions[alignment_index])
		print "Position is: " + str(self.position)
		print 'CIGAR is: ' + str(self.cigar) + 'MAPQ is: ' + str(self.mapq)
		print 'Original Chromosome' + str(self.original_chromosome)
		print 'Phred: ' 
		print self.bases_quality_score

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

	def update_alignment_probability_with_list(self, base_index, likelihoods, weight):
		"""
		For the given index and likelihood list, updates each alignment with the probability and 
		the given weight for allele_frequencies
		"""
		base_opts = self.base_opts
		for align_index, alignment in enumerate(self.alignments):
			# Find probability associated with the base call in that alignment
			
			alignment_base_index = base_opts.index(alignment[base_index])

			# Include alignment allele frequency
			if self.allele_frequencies is not None:
				base_allele_frequency = self.allele_frequencies[align_index][base_index][alignment_base_index]
			else:
				base_allele_frequency = 1

			read_base = self.read[base_index]

			# Don't penalize for unknown SNP if phred score is high
			if self.get_base_quality_score(base_index) > 0.95 and read_base != alignment[base_index] and base_allele_frequency is nextafter(0,1):
				base_prob = likelihoods[alignment_base_index] 
			else:
				base_prob = likelihoods[alignment_base_index] * base_allele_frequency * weight

			self.update_alignment_probability(align_index, base_prob)

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
		self.alignment_probability_scores = [(prob/prob_sum) for prob in self.alignment_probability_scores]

		for index, score in enumerate(self.alignment_probability_scores):
			if score > highest_score:
				highest_score = score
				highest_index = index

			# TODO: Add in tie breaker
			# if score is highest_score:
			# 	# in case of tie, then assign to alignment with highest alignment score
			# 	if self.alignment_quality_scores[index] > self.alignment_quality_scores[highest_index]:
			# 		highest_score = score
			# 		highest_index = index

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