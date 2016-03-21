#!/usr/bin/python
import HTSeq
from Bio.Seq import Seq
import os.path
import argparse


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq,transTable):
	seq=DNASeq
	tableid=transTable
	reversedSeq=False
	try:
		myseq= Seq(seq)
		protseq=Seq.translate(myseq, table=tableid,cds=True)
	except:
		reversedSeq=True
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			protseq=Seq.translate(myseq, table=tableid,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				protseq=Seq.translate(myseq, table=tableid,cds=True)
			except:
				reversedSeq=False
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					protseq=Seq.translate(myseq, table=tableid,cds=True)
				except Exception as e:

					raise ValueError(e)
	return protseq,seq,reversedSeq

def main():
	
	parser = argparse.ArgumentParser(description="This program downloads sequencing runs given the sra RUN ID in a list to a selected directory")
	parser.add_argument('-i', nargs='?', type=str, help='list genes', required=True)
	parser.add_argument('-r', nargs='?', type=bool, help='Return values', required=False)
	
	args=parser.parse_args()

	
	genes = args.i
	
	try:
		ReturnValues=bool(args.r)
	except:
		ReturnValues=False
		pass
		
	analyzeCDS(genes,ReturnValues)

def analyzeCDS(genes,transTable,ReturnValues):
	
	gene_fp = open( genes, 'r')
	
	stopc=0
	notStart=0
	notMultiple=0
	totalalleles=0
	
	
	
	statsPerGene={}
	
	for gene in gene_fp:
		
		listStopc=[]
		listnotStart=[]
		listnotMultiple=[]
		
		
		print "####################"
		print str(os.path.basename(gene))
		
		k=0
		gene = gene.rstrip('\n')
		multiple=True
		gene_fp2 = HTSeq.FastaReader(gene)
		
		# translate each allele and report the error if unable to translate
		for allele in gene_fp2: 
			
			k+=1
			# if allele is not multiple of 3 it's useless to try to translate
			if (len(allele.seq) % 3 != 0):
				multiple=False
				listnotMultiple.append(str(k))
				print "allele "+str(k)+" is not multiple of 3"
				pass
			else:
				try:
					protseq,seq,reversedSeq=translateSeq(allele.seq, transTable)
					
				except Exception, err:
					if "Extra in frame stop codon found" in str(err):
						stopc+=1
						listStopc.append(str(k))
					elif "is not a start codon" in str(err):
						notStart+=1
						listnotStart.append(str(k))
					else:
						print err
					print "allele "+str(k)+" is not translating"
					pass
		statsPerGene[gene]=listnotMultiple,listStopc,listnotStart,k
		totalalleles+=k
		
	print str(stopc) + " alleles have stop codons inside"
	print str(notStart) + " alleles don't have start codons"
	print "total of alleles : " + str(totalalleles)
	
	if not ReturnValues:
		with open("CheckCDSResults.txt", "wb") as f:
			f.write("Alleles with stop codons inside: \n")
			for item in listStopc:
				f.write(item)
				f.write("\n")
			f.write("\nAlleles without start codon: \n")
			for item in listnotStart:
				f.write(item)
				f.write("\n")
	else:
		return statsPerGene
	
if __name__ == "__main__":
    main()

