#!/usr/bin/env python
"""
Creates a synthetic data set to use to test ambio.py
"""
from numpy import random

# sys.stdout.write('%s\n', foo())

# class GeneratedRead:
# 	def __init__(self, nucleotides, score):
# 		"""
# 		nucloetides: a list of nucleotides
# 		score: score associated with the read
# 		"""
# 		self.nucleotides = nucleotides
# 		self.score = score

# class Nucleotide:
# 	def __init__(self, base=None, alleles=None, score=0, read=False):
# 		"""
# 		base: the nucleotide base ATCG
# 		score: the score associated with the base reading
# 		alleles: a list of known alleles
# 		read: True if the nucleotide is from a read, False if the nucleotide is from a genome
# 		"""
# 		self.score = score
# 		self.alleles = alleles
# 		self.base = base

# 	def define_alleles(self, alleles):
# 		"Pass in a list of alleles for this nucleotide"
# 		self.alleles = alleles

# 	def get_alleles(self):
# 		return self.alleles

# 	def set_score(self, score):
# 		"Pass in a float representing the score for this nucleotide"
# 		self.score = score

# 	def get_score(self):
# 		return self.score

# 	def set_base(self):
# 		self.base = base

# 	def get_base(self):
# 		return self.base

# 	def print_nucleotide(self):
# 		print "Nucleotide: " + self.base + " " + str(self.score) 
# 		print self.alleles


def sampler(read_length, genome_length):
	repeat_region = []
	repeat_region_positions = []
	nucleotide = ['A', 'C', 'G', 'T']
	synthetic_genome = []
	reads = []
	sample_genome = []
	sample_repeat_region = []
	templates = []

	# Generate Genome
	for n in range(genome_length):

		# Random draw for nucleotide attributes over a uniform distribution
		nuc_number = random.randint(0,3)
		# allele_number = random.randint(0,1)

		# # Generate New Nucleotide
		# alleles = []
		# for a in range(allele_number):
		# 	allele_number = (nuc_number + allele_number) % 3
		# 	alleles.append(nucleotides[allele_number])

		#new_nucleotide = Nucleotide(nucleotides[nuc_number], alleles)
		synthetic_genome.append(nucleotide[nuc_number])

	# Generate Repetitive Region
	repeat_region = synthetic_genome
	for n in range(0, genome_length, 30):

		# Create an error
		nuc_number = random.randint(0,3)
		repeat_region[n] = nucleotide[nuc_number]

	# Generate Sample Genome
	sample_genome = synthetic_genome
	sample_repeat_region = repeat_region

	# Create SNP into sample
	snps = []
	snps.append(random.randint(0, genome_length)) # in genome, known
	snps.append(random.randint(0, genome_length)) # in genome, known
	snps.append(random.randint(0, genome_length)) # in genome, unknown
	snps.append(random.randint(0, genome_length)) # in repeat, known
	snps.append(random.randint(0, genome_length)) # in repeat, unknown


	# Add SNPs into sample
	for index, s in enumerate(snps):
		if index <= 2:
			original = sample_genome[s]
			typ = 'genome'
			sample_genome[s] = change_snp(original)
			new = sample_genome[s]
		else:
			original = sample_repeat_region[s]
			typ = 'repeat'
			sample_repeat_region[s] = change_snp(original)
			new = sample_repeat_region[s]

		print 'In ' + typ + ' changed index: ' + str(s) + ' of base ' + original + ' to ' +  new



	# Generate Reads 
	# 50 bp from genome, 50 bp from repeat not overlapping
	for n in range(0,genome_length,50):
		read_genome = sample_genome[n:n+50]
		read_repeat = sample_repeat_region[n:n+50]
		template_genome = synthetic_genome[n:n+50]
		template_repeat = repeat_region[n:n+50]

		reads.append(read_genome)
		reads.append(read_repeat)
		templates.append(template_genome)
		templates.append(template_repeat)

	print 'Results'
	print 'Template Genome'
	print synthetic_genome
	print 'Repeat Template'
	print repeat_region
	print 'Sample Genome'
	print sample_genome
	print 'Sample Repeat'
	print sample_repeat_region
	print 'Templates'
	print templates
	return reads


def change_snp(nuc):
	if nuc is 'A':
		return 'G'
	if nuc is 'C':
		return 'T'
	if nuc is 'G':
		return 'A'
	if nuc is 'T':
		return 'C'

print sampler(50, 200)







