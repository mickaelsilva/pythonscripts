#!/usr/bin/python
import HTSeq
import sys
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
from Counter import Counter 		# --- Counter.py from python3 is needed to run this script --- #
#from collections import Counter
import os
import os.path
import string
import argparse
import subprocess
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from Genes import Gene
from Genes import SetOfGenes
from multiprocessing import Pool
import time
import shutil



def runProdigal(contigsFasta):

        # ------------ #
        # RUN PRODIGAL #
        # ------------ #  #/scratch/NGStools/prodigal-2.50/prodigal -i /scratch/spneumoniae_600/ERR067978/velvet_assembly/contigs.fa -c -m -g 11 -p single -f sco -q > test_all.txt
	prodigal_path='prodigal'
	proc = subprocess.Popen([prodigal_path, '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q'], stdout=subprocess.PIPE)
	
	cdsDict = {}
	tempList = []
	line = ' '
	while line != '':
		
		# when it finds a contig tag
		if "seqhdr" in line:
		# add contig to cdsDict and start new entry
			
			if len(tempList) > 0:

				 # --- brute force parsing of the contig tag - better solution is advisable --- #
				
				i=0
				for l in contigTag:
					if l == ' ':
						break
					i+=1
				contigTag=contigTag[:i]
				
				cdsDict[contigTag.replace("\r","")] = tempList
				tempList = []
			
			contigTag = line.split('"')[-2]
			
                # when it finds a line with cds indexes
		elif line[0] == '>':
			
                        # parsing
			cdsL = line.split('_')

                        # --- each element of this list is a pair of indices - the start and the end of a CDS --- #

			tempList.append([ int(cdsL[1]) - 1 , int(cdsL[2]) ])		# start index correction needed because prodigal indexes start in 1 instead of 0
			#print tempList
                # reads the stdout from 'prodigal'
		line = proc.stdout.readline()

	# ADD LAST
	if len(tempList) > 0:
		
		 # --- brute force parsing of the contig tag - better solution is advisable --- #

		i=0
		for l in contigTag:
			if l == ' ':
				break
			i+=1
		contigTag=contigTag[:i]
		cdsDict[contigTag.replace("\r","")] =tempList
	#print cdsDict.keys()
	return cdsDict

def getBlastScoreRatios(genefile):
	
	gene_fp = HTSeq.FastaReader(genefile)
	alleleI=0
	allelescores=[]
	alleleProt=''
	for allele in gene_fp: #new db for each allele to blast it against himself
		alleleI+=1
		genome=-1
		alleleProt+=">"+str(alleleI)+"\n"+str(translateSeq(allele.seq)+"\n")
	basepath="./blastdbs/temp"+str(os.path.basename(genefile))
	if not os.path.exists(basepath):
		os.makedirs(basepath)
	with open(basepath+'/protein.fasta', "wb") as f:
		f.write(alleleProt)
	Gene_Blast_DB_name = Create_Blastdb( basepath+'/protein.fasta', 1, True )
		# --- get BLAST score ratio --- #
	cline = NcbiblastpCommandline(query=basepath+'/protein.fasta', db=Gene_Blast_DB_name, evalue=0.001, out=basepath+'protein.xml', outfmt=5)
		#print cline
	allelescore=0
	blast_records = runBlastParser(cline,basepath+'protein.xml', alleleProt)
	
	for blast_record in blast_records:
		found=False 
		for alignment in blast_record.alignments:
			if found is False:
				#print blast_record.query, alignment.hit_def
				for match in alignment.hsps:
					
						#print "---------------------"
					if(int(alignment.hit_def)== int(blast_record.query)):
						#print match
						allelescores.append(int(match.score))
						found=True
						break
			else:
				break
	#print allelescores
	return alleleI,allelescores,Gene_Blast_DB_name

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	try:
		myseq= Seq(seq)
		#print myseq
		protseq=Seq.translate(myseq, table=11,cds=True)
	except:
		try:
			seq=reverseComplement(seq)
			myseq= Seq(seq)
			#print myseq
			protseq=Seq.translate(myseq, table=11,cds=True)
						
		except:
			try:
				seq=seq[::-1]
				myseq= Seq(seq)
				#print myseq
				protseq=Seq.translate(myseq, table=11,cds=True)
			except:
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					#print myseq
					protseq=Seq.translate(myseq, table=11,cds=True)
				except:
					raise
	return protseq
	

def printinfo(genome, gene):
	print "_______________________________________________________"
	print "Genome : "+ str(os.path.basename(genome)) 
	print "Locus : "+ str(os.path.basename(gene)) 

# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def callAlleles(argumentList):
	
	geneFile = argumentList[0]
	genomesList = argumentList[1]
	listOfProts = argumentList[2]
	listAllCDS = argumentList[3]
	#print geneFile
	gene_fp = HTSeq.FastaReader(geneFile)
	alleleI = 0
	#inverted=False
	#orderedAlleleNames=[]
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	bestmatch=[0,0,False,'',0] #score, score ratio, perfectmatch, key name of the DNA sequence string, allele ID
	allelescores=[]
	
	
	alleleI,allelescores,Gene_Blast_DB_name=getBlastScoreRatios(geneFile)
	
	genome=-1	
	
	for protList in listOfProts:

		#alleleI = 0
		#alleleProt=''
		#for allele in gene_fp: #new db for each allele to blast it against himself
		#	alleleI+=1
		#	alleleProt+=">"+str(alleleI)+"\n"+str(translateSeq(allele.seq)+"\n")
		basepath="./blastdbs/temp"+str(os.path.basename(geneFile))
		#if not os.path.exists(basepath):
		#	os.makedirs(basepath)
		#with open(basepath+'/protein.fasta', "wb") as f:
		#	f.write(alleleProt)
		#Gene_Blast_DB_name = Create_Blastdb( basepath+'/protein.fasta', 1, True )
		genome+=1
		with open(basepath+'/proteinList.fasta', "wb") as f:
			f.write(protList)
		#Gene_Blast_DB_name = Create_Blastdb( './temp/proteinList.fasta', 1, True )
		cline = NcbiblastpCommandline(query=basepath+'/proteinList.fasta', db=Gene_Blast_DB_name, evalue=0.001, out=basepath+'proteinList.xml', outfmt=5)
		#print cline
		blast_records = runBlastParser(cline, basepath+'proteinList.xml', basepath+'/proteinList.fasta')
		for blast_record in blast_records:
				
			for alignment in blast_record.alignments:
				#print alignment
					#print alignment.hsps
				#print alignment.hit_id
				#print alignment.hit_def
					#print alignment.title
				for match in alignment.hsps:
					#print blast_record.query
					#print match
					#print alleleI, len(allelescores)
					scoreRatio=float(match.score)/float(allelescores[int(alignment.hit_def)-1])
					#print scoreRatio
					#print alignment.hit_def
					cdsStrName=blast_record.query
					if(scoreRatio == 1 and bestmatch[2] is False):
						bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alignment.hit_def)]
						#print alignment
						#print match
					elif(scoreRatio == 1 and match.score>bestmatch[0]):
						bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alignment.hit_def)]
						#print match
					elif(match.score>bestmatch[0] and scoreRatio>0.4 and scoreRatio>bestmatch[1] and bestmatch[2] is False):
						#print match.query
						#print match.sbjct
						#print allelescores
						bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alignment.hit_def)]
						#print match
		#print bestmatch
				
		if bestmatch[0]==0:
					#if no best match was found
					
					###################
					# LOCUS NOT FOUND #
					###################
						
			resultsList.append('LNF3:-1')            # append result to the list of results
			perfectMatchIdAllele.append('LNF')
			#printinfo(genomeFile,geneFile)
			print "Locus not found, no matches \n"
			
				
		elif bestmatch[2] is True:
						
					#if a perfect match was found
					
						################################################
						# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
						################################################
					
			perfectMatchIdAllele.append(str(bestmatch[4]))
			resultsList.append('EXC:' + str(bestmatch[4]) )
				
		else:
					#######################
					# ADD INFERRED ALLELE #		# a new allele 
					#######################
					
											
				#print "infered allele has location : "+(CDSType)
				#printinfo(genomeFile,geneFile) 
			tagAux='INF'
			perfectMatchIdAllele.append( tagAux +"-"+str(alleleI+1))
			print "New allele! Adding allele "+ tagAux + str(alleleI+1) +" to the database\n"
																				
			resultsList.append( tagAux + str(alleleI+1) )

					#orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] +"_" +str(os.path.basename(genomeFile)))	
										# --- add the new allele to the gene fasta --- #

			fG = open( geneFile, 'a' )
			fG.write('>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n')
					#print alleleStr
				
			listOfCDS=listAllCDS[genome]
			#print listOfCDS
			fG.write( listOfCDS[">"+bestmatch[3]] + '\n')
			fG.close()
					#alleleI += 1
					# --- remake blast DB --- #
			Gene_Blast_DB_name = Create_Blastdb( basepath+'/protein.fasta', 1, True )
			alleleI,allelescores,Gene_Blast_DB_name=getBlastScoreRatios(geneFile)
	#x=y
	shutil.rmtree(basepath)

	
	final =	(resultsList,perfectMatchIdAllele)	
	#return (resultsList)
	return final



# ================================================ MAIN ================================================ #

def main():
	
	
	
	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('-p', nargs='?', type=str, help="True to give a phyloviz output file type, false is predefined", required=False)

	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	phylovizinput=False
	if(args.p):
		phylovizinput=args.p
	# ------------------------------------------------- #
	#           RUN PRODIGAL OVER ALL GENOMES           #
	# ------------------------------------------------- #

	print ("Starting Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	listOfCDSDicts = []
	listOfGenomes = []
	listOfGenomesDict = []

	fp = open(genomeFiles, 'r')

	for genomeFile in fp:

		genomeFile = genomeFile.rstrip('\n')
		genomeFile = genomeFile.rstrip('\r')
		listOfGenomes.append( genomeFile )
		genomeDict = {}
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:

			genomeDict[ contig.name ] = (contig.seq)

		listOfGenomesDict.append( genomeDict )

	fp.close()


	print ("Starting Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	poolJobs = Pool()

	totgenomes= len(listOfGenomes)
	
	it = poolJobs.imap(runProdigal, listOfGenomes)
	pc1=0
        while 1:
		try:
			pc1+=1
			print (str(pc1)+"/"+str(totgenomes))
			listOfCDSDicts.append( it.next() )
			
		except StopIteration:
			break

	poolJobs.close()
	poolJobs.join()

	print ("Finishing Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	# ----------------------------- #
	# Each gene has a different job #
	# ----------------------------- #
	
	#---CDS to protein---#
	listOFProt=[]
	listAllCDS=[]
	
	i = 0
	j=0
	for genomeFile in listOfGenomes:
		listOfCDS={}
		genomeProts=""
		currentCDSDict = listOfCDSDicts[i]
		#print currentCDSDict
		currentGenomeDict = listOfGenomesDict[i]
		i+=1
		for contigTag,value in currentCDSDict.iteritems():
			#print contigTag,value
			for protein in value:
				genomeProts+=">"+contigTag+"|protein"+str(j)+"\n"
				
				#print contigTag, protein[0], protein[1]
				#print currentGenomeDict[ contigTag ]
				seq= currentGenomeDict[ contigTag ][ protein[0]:protein[1] ].upper()
				listOfCDS[">"+contigTag+"|protein"+str(j)]=seq
				protseq=translateSeq(seq)
				genomeProts+=str(protseq)+"\n"
				j+=1
		listAllCDS.append(listOfCDS)
		listOFProt.append(genomeProts)
			#seq=listOfGenomesDict[ contigTag ][ value[0]:value[1] ].upper()
	#print listOFProt
	
	
	argumentsList = []
	lGenesFiles = []
	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	gene_fp = open( genes, 'r')
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		lGenesFiles.append( gene )
		argumentsList.append( [gene, listOfGenomes, listOFProt,listAllCDS] )

#		callAlleles([gene, listOfGenomes, listOfCDSDicts, listOfGenomesDict])
	gene_fp.close()

	#output=callAlleles([gene, listOfGenomes, listOfCDSDicts, listOfGenomesDict])
	#print output
#	raise SystemExit

	totloci= len(argumentsList)
	
	poolJobs = Pool()
	pc2=0
	it = poolJobs.imap(callAlleles, argumentsList)
	output =[]

	while 1:
		try:
			pc2+=1
			print (str(pc2)+"/"+str(totloci))
			output.append( it.next() )
			

		except StopIteration:
			break

	poolJobs.close()
        poolJobs.join()

	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	

	print "##################################################\n %s genomes used for %s loci" % (len(output[0][0]),len(output) )
	numberexactmatches=0
	for gene in output:
		for gAllele in gene[0]:
			if("EXC:" in gAllele):
				numberexactmatches+=1
					
				
	print "\n %s exact matches found out of %s" % (numberexactmatches,(len(output[0][0])*len(output)) )	
	print "\n %s percent of exact matches \n##################################################" % (float((numberexactmatches*100)/(len(output[0][0])*len(output))) )	
		
	print "\nWriting output files\n"
	args.o = '/' + args.o
	
	if(phylovizinput is False):
		genesI=0
		for geneOut in output:
			i=0
			for gAllele in geneOut[0]:
				
				currentGenome = listOfGenomes[i]
				currentGenome=currentGenome.split("/")
				currentGenome=currentGenome[len(currentGenome)-1].split(".")
				gOutFile = os.path.dirname( "./" )
				finalname=(args.o).split("/")
				gOutFile += "/"+str(currentGenome[0])+finalname[1]				
				if not os.path.isfile( gOutFile )or (i==0 and genesI==0):

					aux = 'w'
				else:
					aux = 'a'
					gAllele = '\n' + gAllele
					
					
				f = open(gOutFile, aux)
				f.write(gAllele + ':' + lGenesFiles[genesI])
				f.close()
				i+=1
					
			genesI+=1
	
	else:
		try:
			phylovout=[]
			genesnames=[]
			statistics=[]
			
			for gene in lGenesFiles:
				
				genename=gene.split("/")
				#genename=genename[len(genename)-1].split(".")
				genename=genename[len(genename)-1]
				genesnames.append(genename)
			for geneOut in output:
				gene=0
				alleleschema=[]
				while gene<len(output[0][0]): 

					genename=(geneOut[1][gene]).split("_")

					if(len(genename)!=1):
						alleleschema.append(genename[1])
					else:
						alleleschema.append(genename[0])
					gene+=1
				phylovout.append(alleleschema)
			
			genome=0
			finalphylovinput= "FILE"+ "\t" 
			for geneid in genesnames:
				finalphylovinput+= str(geneid)+ "\t"
				
			
			while genome<len(listOfGenomes):
				currentGenome = os.path.basename(listOfGenomes[genome])
				statsaux=[0]*6 # EXC INF LNF LOT incomplete SAC
				finalphylovinput+= "\n" + currentGenome + "\t"
				for gene in phylovout:
					
					val= str(gene[genome])
					finalphylovinput+= val + "\t"
					if "INF" in val:
						statsaux[1]+=1
					elif "LNF" in val:
						statsaux[2]+=1
					elif "LOT" in val:
						statsaux[3]+=1
					elif "incomplete" in val:
						statsaux[4]+=1 
					elif "small" in val:
						statsaux[5]+=1
					else:
						statsaux[0]+=1
					
				genome+=1
				statistics.append(statsaux)
			gOutFile = os.path.dirname( "./")
			gOutFile += args.o
			statswrite='Stats:\tEXC\tINF\tLNF\tLOT\tincomplete\tsmall'
			i=0
			genome=0
			while genome<len(listOfGenomes):
				currentGenome = os.path.basename(listOfGenomes[genome])
				statsaux=[0]*6 # EXC NA INF LNF LOT incomplete SAC
				statswrite+= "\n" + currentGenome + "\t"
				for k in statistics[i]:
					statswrite+= str(k) + "\t"
				i+=1	
				genome+=1
					
			print statswrite
			with open(gOutFile, 'w') as f:
				f.write(finalphylovinput)
			statoutfile=os.path.dirname( "./")
			with open("stastics.txt", 'w') as f:
				f.write(str(statswrite))
		except Exception as e:
			print e
			gOutFile = os.path.dirname( "./")
			gOutFile += args.o
			with open(gOutFile, 'w') as f:
				f.write(str(output))
			
		

	print ("Finished Script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

if __name__ == "__main__":
    main()
