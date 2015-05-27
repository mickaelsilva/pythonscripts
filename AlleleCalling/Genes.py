#! /usr/bin/env python

class Gene:
	"""Represents a gene"""
	def __init__(self, tag, description, name, tigr4, seq):
		self.tag = tag
		self.description = description
		self.name = name
		self.tigr4 = tigr4
		self.seq = seq.upper()
		self.length = len(seq)
		
	def toString(self):
		print self.tag, self.description, self.name, self.tigr4, self.seq, self.length

class SetOfGenes:
	"""Represents a set of genes"""
	def __init__(self):
		self.genesDict = {}
		self.nGenes = 0

# old version		
#	def addGene(self, gene):
#		self.genesDict[gene.tag] = gene
#		self.nGenes += 1
	
	def addGene(self, gene):
		self.genesDict[gene.seq] = gene
		self.nGenes += 1

# old version	
#	def addUniqueGene(self, gene):
#		if self.nGenes != 0:
#			for key, value in self.genesDict.iteritems():
#				if value.tag == gene.tag and value.description == gene.description and value.seq == gene.seq:
#					return
#		
#		self.genesDict[gene.tag] = gene
#		self.nGenes += 1

	def addUniqueGene(self, gene):
		if self.nGenes != 0:
			if gene.seq in self.genesDict:
				return

		self.genesDict[gene.seq] = gene
		self.nGenes += 1

# old version - no use		
#	def getGene(self, key):
#		return self.genesDict[key]
		
	def genesToFile(self, fName):
		f = open(fName, 'w')
		for key, gene in self.genesDict.iteritems():
			f.write("> "+gene.tag+' '+gene.description+"\n"+gene.seq+"\n")
		f.close()
