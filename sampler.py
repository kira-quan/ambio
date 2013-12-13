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


def make_read_file(num_reads, genome_length):

	read_file = open('generated_reads.txt', 'w')
	repeat_file = open('generated_regions.txt', 'w')
	snp_file = open('generate_snps.txt', 'w')

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
	repeat_region = synthetic_genome[:]
	for n in range(0, genome_length, 30):
		# Create an error
		nuc_number = random.randint(0,3)
		repeat_region[n] = nucleotide[nuc_number]

	# Generate Sample Genome
	sample_genome = synthetic_genome[:]
	sample_repeat_region = repeat_region[:]

	# Create SNP into sample
	snps = []
	snps.append(random.randint(0, genome_length)) # in genome, known
	snps.append(random.randint(0, genome_length)) # in genome, known
	snps.append(random.randint(0, genome_length)) # in genome, unknown
	snps.append(random.randint(0, genome_length)) # in repeat, known
	snps.append(random.randint(0, genome_length)) # in repeat, unknown
	snps.append(random.randint(0, genome_length)) # in repeat, unknown
	

	known_snps = [snps[0], snps[1], snps[3]]
	known_snp_info = []


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

		if index is 0 or index is 1 or index is 3:
			snp_gen = original + '/' + new
			known_snp_info.append(snp_gen)

	# Generate reads
	for generated_read_num in range(0, num_reads):
		chromosome = 'chr11'

		# Get region
		region = generated_read_num % 2

		if region is 0:
			# get from sample genome
			position = 100
			read_list = sample_genome[:]
			base_quality = 'JJIHFJIGJIJJJIIIIJJIHGHJJIJIHAHFAJJJJHHHHHFFFFFCCC'
			cigar = '50M'
			quality = 255

			# add in base miscalls
			num_miscalls = random.randint(1,4)

			for n in range(0,num_miscalls):
				index = random.randint(0, 49)
				new_base = nucleotide[random.randint(0,3)]
				read_list[index] = new_base
				
			sequence = ''.join(read_list)
			

		else:
			# get from repeat genome
			position = 200
			read_list = sample_repeat_region[:]
			base_quality = 'JJIHFJIGJIJJJIIIIJJIHGHJJIJIHAHFAJJJJHHHHHFFFFFCCC'
			cigar = '50M'
			quality = 255

			# add in base miscalls
			num_miscalls = random.randint(1,2)

			for n in range(0,num_miscalls):
				index = random.randint(0, 49)
				new_base = nucleotide[random.randint(0,3)]
				read_list[index] = new_base

			sequence = ''.join(read_list)

		write_read = chromosome + ' ' + str(position) + ' ' + sequence + ' ' + base_quality + ' ' + cigar + ' ' + str(quality) + '\n'
		read_file.write(write_read)

	# Generate repeats
	repeat_score = 100
	repeat_first = ''.join(synthetic_genome)
	repeat_second = ''.join(repeat_region)
	chromosome = 'chr11'

	write_repeat_genome = repeat_first + ' ' + str(repeat_score) + ' ' + repeat_first + ' ' + chromosome + ' ' + str(100) + '\n'
	write_repeat_region = repeat_first + ' ' + str(repeat_score) + ' ' + repeat_second + ' ' + chromosome + ' ' +str(200)
	repeat_file.write(write_repeat_genome)
	repeat_file.write(write_repeat_region)

	# Generate SNP file
	genome_rs = 'rs100'
	repeat_rs = 'rs200'
	
	known_snp_one = genome_rs + ' ' + str(11) + ' ' + str(100 + known_snps[0]) + ' ' + known_snp_info[0] + '\n'
	known_snp_two = genome_rs + ' ' + str(11) + ' ' + str(100 + known_snps[1]) + ' ' + known_snp_info[1] + '\n'
	known_snp_three = repeat_rs + ' ' + str(11) + ' ' + str(200 + known_snps[2]) + ' ' + known_snp_info[2]

	snp_file.write(known_snp_one)
	snp_file.write(known_snp_two)
	snp_file.write(known_snp_three)	

	read_file.close()
	repeat_file.close()
	snp_file.close()

make_read_file(1000, 50)