#!/usr/bin/env python
"""
Creates a synthetic data set to use to test ambio.py
"""
from numpy import random

class Read:
	def __init__(self, nucleotides, score):
		"""
		nucloetides: a list of nucleotides
		score: score associated with the read
		"""
		self.nucleotides = nucleotides
		self.score = score

class Nucleotide:
	def __init__(self, base=None, alleles=None, score=0, read=False):
		"""
		base: the nucleotide base ATCG
		score: the score associated with the base reading
		alleles: a list of known alleles
		read: True if the nucleotide is from a read, False if the nucleotide is from a genome
		"""
		self.score = score
		self.alleles = alleles
		self.base = base

	def define_alleles(self, alleles):
		"Pass in a list of alleles for this nucleotide"
		self.alleles = alleles

	def get_alleles(self):
		return self.alleles

	def set_score(self, score):
		"Pass in a float representing the score for this nucleotide"
		self.score = score

	def get_score(self):
		return self.score

	def set_base(self):
		self.base = base

	def get_base(self):
		return self.base

	def print_nucleotide(self):
		print "Nucleotide: " + self.base + " " + str(self.score) 
		print self.alleles


def sampler(read_length, genome_length):
	repeat_region_length = 200
	repeat_region_positions = []
	nucleotides = ['A', 'T', 'C', 'G']
	synthetic_genome = []
	reads = []

	for n in range(genome_length):

		# Random draw for nucleotide attributes over a uniform distribution
		nuc_number = random.randint(0,3)
		allele_number = random.randint(0,2)

		# Generate New Nucleotide
		alleles = []
		for a in range(allele_number):
			allele_number = (nuc_number + allele_number) % 3
			alleles.append(nucleotides[allele_number])

		new_nucleotide = Nucleotide(nucleotides[nuc_number], alleles)
		synthetic_genome.append(new_nucleotide)

	






