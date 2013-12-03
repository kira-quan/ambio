#!/usr/bin/env python
# This file defines the Read class for a read from a DNA alignment. This class contains the read,
# associated alignments, scores for those alignments and bases, and the position to align the 
# read

class Read:
	def __init__(self, read, alignments, alignment_quality_scores, bases_quality_score):
		# the read
		self.read = read
		# a list of the possible alignments for a read
		self.alignments = alignments 
		# a list of the corresponding quality scores of each alignment, same order as alginments
		self.alignment_quality_scores = alignment_quality_scores
		# a list of the quality score for each base (in order) in a read
		self.bases_quality_score = bases_quality_score
		# the best position for the read, where the position refers to the alignment index
		self.position = 0
		# a list of generate probability score for each alignment
		self.alignment_probability_scores = [0 for a in alignments]

	def print_read(self):
		print "Read is: " + self.read
		print "Alignments are: "
		for alignment_index, alignment in enumerate(self.alignments):
			print alignment
			print "Score: " + str(self.alignment_probability_scores[alignment_index])
		print "Position is: " + str(self.position)

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

	def get_read(self):
		return self.read

	def find_highest_position(self):
		"""
		Sets the position to the alignment with the highest score
		"""
		highest_index = 0
		highest_score = -1
		for index, score in enumerate(self.alignment_probability_scores):
			if score > highest_score:
				highest_score = score
				highest_index = index
		self.set_position(highest_index)