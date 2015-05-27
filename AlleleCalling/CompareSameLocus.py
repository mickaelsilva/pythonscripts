import HTSeq
import argparse
import os.path
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from Bio.Blast.Applications import NcbiblastnCommandline
import shutil

def concat_genes(listgenes, outFileName):
	fp = open(listgenes, 'r')
	concatGenes=''
	for genefile in fp:
		genefile = genefile.rstrip('\n')
		biggestallename=''
		biggestallelestr=''
		biggestallelelen=0
		gene_fp = HTSeq.FastaReader(genefile)
		for contig in gene_fp:
			contigLen=int(len(contig.seq))
			if contigLen>biggestallelelen:
				#biggestallename=contig.name+" "+contig.descr
				biggestallelestr = contig.seq
		#concatGenes+=">"+str(os.path.basename(genefile))+" | "+biggestallename+"\n"+biggestallelestr+"\n"
		concatGenes+=">"+genefile+"\n"+biggestallelestr+"\n"
	
	with open(outFileName, "wb") as f:
		f.write(concatGenes)

def alignHasGoodMatch(alignment,geneDict,LocusID,blast_record,alreadyUsed):
	sameAllele=0
	for match in alignment.hsps:
					bmAlleleLen=0
					bmAllele=''
					for seq, alleleid in geneDict.iteritems():
						if alleleid == alignment.hit_def:
							bmAllele=seq
							bmAlleleLen= len(bmAllele)
							break
						
					
					idratio=float(match.identities)/float(bmAlleleLen)
					
					sizeratio=float(match.align_length)/float(bmAlleleLen)
					
					
					if sizeratio>0.8 and sizeratio<1.2 and idratio>0.8 and alignment.hit_def not in alreadyUsed :
						LocusID+=1
						
						genename=alignment.hit_def.split("/")
						genename=genename[len(genename)-1]
						
						newpath='./sameLocus/L'+str(LocusID)+'_'+genename
						
						shutil.copy2(alignment.hit_def, newpath)
						
						#genename=blast_record.query.split("/")
						#genename=genename[len(genename)-1]
						#newpath='./sameLocus/L'+str(LocusID)+'_'+genename
						#shutil.copy2(blast_record.query, newpath)
						
						alreadyUsed.append(alignment.hit_def)
						
						isGood=True
						sameAllele+=1
	return 	sameAllele,alignment.hit_def,LocusID,alreadyUsed

def main():

	parser = argparse.ArgumentParser(description="Given two list of genes, creates a folder with paired files when located on the same locus")
	parser.add_argument('-i', nargs='?', type=str, help='1st list of genes files to compare', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='2nd list of genes files to compare', required=True)
	
	args = parser.parse_args()
	geneFiles1 = args.i
	geneFiles2 = args.g
	
		
	name1="concat1.fasta"
	name2="concat2.fasta"
		
	concat_genes(geneFiles1, name1)
	concat_genes(geneFiles2, name2)
	
	#orderedAlleleNames=[]

	geneDict={}
	gene_fp = HTSeq.FastaReader(name1)
	alleleI=0
	for allele in gene_fp:
		#if allele.seq in geneDict:
		#	print "\nWARNING: this file contains a repeated allele, it should be checked. Ignoring it now!\n", geneFile
		#else:
			#orderedAlleleNames.append(allele.name)
		geneDict[ allele.seq ] = allele.name
		alleleI += 1
	
	gene_fp = HTSeq.FastaReader(name1)
	geneFile = os.path.abspath( name1 )
	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 , False)
	
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'

	# list of results - the output of the function
	resultsList = []

					# ------------------------------ RUNNING BLAST ------------------------------ #

	cline = NcbiblastnCommandline(query=name2, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline, blast_out_file, name2)
	samelocus=0
	alreadyUsed=[]
	nomatch=0
	small=0
	if not os.path.exists("./sameLocus"):
		os.makedirs("./sameLocus")
	LocusID=0
	for blast_record in blast_records:
		try:
			alignment=blast_record.alignments[1]
			#print blast_record.query
			#print alignment.num_alignments
			
			try:
					#print alleleLength, alignment.length
				i=0
				align=blast_record.alignments[i]	
				while i<len(blast_record.alignments):
					if align.hit_def:
						result,allelename2,LocusID,alreadyUsed=alignHasGoodMatch(align,geneDict,LocusID,blast_record,alreadyUsed)
						if result>0 and allelename2:
							samelocus+=result
							i+=999
						else:
							small+=1
							i+=999
							alreadyUsed.append(allelename2)
					elif allelename :
						#alreadyUsed.append(allelename)
						result,allelename,LocusID,alreadyUsed=alignHasGoodMatch(align,geneDict,LocusID,blast_record,alreadyUsed)
						if result>0:
							samelocus+=result
							i+=999
						else:
							small+=1
							i+=999
							#alreadyUsed.append(allelename2)
					else :
						nomatch+=1
					#print align.length, alleleleng
					
					i+=1
			except Exception as e:
				print e
					#print "lkjh"
				pass
		except:
			try:
				alignment=blast_record.alignments[0]
				#print blast_record.query
				
				result,allelename,LocusID,alreadyUsed=alignHasGoodMatch(alignment,geneDict,LocusID,blast_record,alreadyUsed)
				if result>0 and allelename:
					samelocus+=result
				else :
					small+=1
				#alreadyUsed.append(allelename)
				#alreadyUsed.append(alignment.hit_def)
			except:
				nomatch+=1
				
	
	print "%s are within same locus, %s had no match and %s had a bigger than 0.2 ratio size difference or less than 0.8 similarity ratio" % (samelocus,nomatch, small)
	
	os.remove(name1)
	os.remove(name2)
	shutil.rmtree('./blastdbs')
	
if __name__ == "__main__":
	main()
