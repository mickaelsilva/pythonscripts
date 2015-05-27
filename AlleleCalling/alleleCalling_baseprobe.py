#!/usr/bin/python
import HTSeq
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
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


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]


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
	listOfGenomesDict = argumentList[2]
	
	gene_fp = HTSeq.FastaReader(geneFile)
	geneDict = {}
	alleleI = 1
	orderedAlleleNames=[]
	biggestAllelelen=0
	smallestAllelelen=99999
	for allele in gene_fp:
		if allele.seq in geneDict:
			print "\nWARNING: this file contains a repeated allele, it should be checked. Ignoring it now!\n", geneFile
		else:
			if len(allele.seq)>biggestAllelelen:
				biggestAllelelen=len(allele.seq)
			if len(allele.seq)<smallestAllelelen:
				smallestAllelelen=len(allele.seq)
			orderedAlleleNames.append(allele.name)
			geneDict[ allele.seq ] = alleleI
			alleleI += 1
	#print geneDict
	#print orderedAlleleNames

	# --- make 1st blast DB --- #

	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'

	# list of results - the output of the function
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	
	for genomeFile in genomesList:
		#print geneDict
		currentGenomeDict = listOfGenomesDict[i]
		
		#print genomeFile
		#print resultsList
		#print geneDict
		#print orderedAlleleNames
		i+=1		# it has to be incremented here
		if genomeFile[-1] == '\n':
			genomeFile = genomeFile[:-1]


                # ------------------------------ RUNNING BLAST ------------------------------ #

		cline = NcbiblastnCommandline(query=genomeFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		blast_records = runBlastParser(cline, blast_out_file, genomeFile)
		

		# ------ DETERMINING BEST MATCH ------ #

		# bestMatch = ['rec.query','hsp', lenRatio]
		bestMatch = ['','', 0]
		bestMatchContig=''
		bestMatchContigLen=''
		bestalignlen=0
		perfectMatch=False
		bmAlleleLen2=0
		bmAllele=''
		#noAlignment=False
		for blast_record in blast_records:
			# --- the LNF cases are now called outside de loop --- #
			#print blast_record
			if perfectMatch==True:
				break
			try:
				#print blast_record.alignments
				hspC = blast_record.alignments[0]
				
				if bestMatch[0] == '' and bestMatch[1] == '':
					bestMatch[0] = blast_record.query
					bestMatch[1] = hspC
			except IndexError:
				continue


			# --- the contig tag is used in the progigal function --- #

			contigTag = blast_record.query
			# --- brute force parsing of the contig tag - better solution is advisable --- #			

			j=0
			for l in contigTag:
				if l == ' ':
					break
				j+=1

			contigTag = contigTag[:j]

			contigLen = blast_record.query_letters
			#print blast_record.query_id
			
			# --- iterating over all the results to determine the best match --- #
			for alignment in blast_record.alignments:

				index=orderedAlleleNames.index(alignment.hit_def)
				#print alignment.hit_def
				for k, v in geneDict.iteritems():
					if v == index+1:
						bmAlleleLen2= len(k)
					
				if perfectMatch:
					break
				for match in alignment.hsps:
					#print match
					scoreRatio = float(match.score) / float(bmAlleleLen2)
					
					#print alignment.hit_def
					#print match.identities
					#print bmAlleleLen2
					#print
					#print match.identities
					#print len(match.sbjct)
					
					#if #identities is the same as the length of the allele and it has no gaps or N's
					if (int(match.identities)==int(bmAlleleLen2) and int(match.identities)==int(len(match.sbjct)) and "N" not in match.query ): 
						index=orderedAlleleNames.index(alignment.hit_def)
						
						for seq, alleleid in geneDict.iteritems():
							if alleleid == index+1:
								bmAllele=seq
								break
						bmAlleleLen= len(bmAllele)
						
						lenratio=float(len(match.query))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, alignment.hit_def,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						perfectMatch=True
						index=orderedAlleleNames.index(alignment.hit_def)
						bmAlleleLen= len(geneDict.keys()[index])
						break
						
					#choose the match with the best score ratio score/length of allele	
					elif scoreRatio > bestMatch[2]:
						index=orderedAlleleNames.index(alignment.hit_def)
						#print orderedAlleleNames
						#print geneDict
						#print orderedAlleleNames
						#print alignment.hit_def
						#print index
						#print geneDict
						for seq, alleleid in geneDict.iteritems():
							if alleleid == index+1:
								bmAllele=seq
								break
						bmAlleleLen= len(bmAllele)
						#print bmAllele
						lenratio=float(len(match.query))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, alignment.hit_def,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						bestMatchContigLen=blast_record.query_letters
						
						#print match.query
						bestalignlen=alignment.length
						
						
					
					if perfectMatch==True:
						break


		# ---------- ALLELE CALLING AFTER DETERMINING BEST MATCH ---------- #

		
		
		try:
			#print bestMatch[0]
			match = bestMatch[1]
			#print match
			#print match.sbjct
			geneLen = bestMatch[5]
			alleleStr = match.query
			nIdentities = match.identities
			idPercent = float(nIdentities) / float(geneLen)
			scoreRatio = bestMatch[2]
			lenRatio = bestMatch[4]
			#print perfectMatch
			#print "\nContig best/exact match is :"
			#print bestMatchContig +"\n"
		
		except:
			#if no best match was found
			
			###################
			# LOCUS NOT FOUND #
			###################
			
			resultsList.append('LNF:-1')            # append result to the list of results
			perfectMatchIdAllele.append('LNF')
			printinfo(genomeFile,geneFile)
			print "Locus not found \n"
			continue
		
		
		if perfectMatch is True:
			
			#if a perfect match was found
			
			if match.sbjct_start > match.sbjct_end: #reverse the order if needed
				alleleStr = reverseComplement(alleleStr)
			alleleNumber = geneDict[ alleleStr ]
			
					################################################
					# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
					################################################
			if "_" in bestMatch[3]:
				a=bestMatch[3].split("_")
				perfectMatchIdAllele.append(a[1])
			else:
				perfectMatchIdAllele.append(bestMatch[3])
			resultsList.append('EXC:' + str(alleleNumber) )
								
			

		elif bestMatch[0] != '' and perfectMatch is not True:
						
			#if a best match was found but it's not an exact match	

					###########################
					# LOCUS ON THE CONTIG TIP #
					###########################
						
			
			if match.query_start ==1 and len(match.query) < geneLen:
			
				resultsList.append('LOT5:-1')
				perfectMatchIdAllele.append('LOT5')
				printinfo(genomeFile,geneFile)
				
				print "Locus is on the 5' tip of the contig \n"
			
			
			elif match.query_end == bestMatchContigLen and len(match.query) < bestMatchContigLen:
				
				resultsList.append('LOT3:-1')
				perfectMatchIdAllele.append('LOT3')
				printinfo(genomeFile,geneFile)

				print "Locus is on the 3' tip of the contig \n"
			
			elif bestMatchContigLen <= geneLen:
				
				resultsList.append('LOTSC:-1')
				perfectMatchIdAllele.append('LOTSC')
				printinfo(genomeFile,geneFile)
				#print match.query_start
				print "Locus is bigger than the contig \n"
								
				
			elif 'N' in alleleStr:
				#TODO gravar para ficheiro
					#####################
					# ALLELE NOT FOUND  #		# N base found!
					#####################
				
				geneFile2= os.path.splitext(geneFile)[0] + "LNFN.fasta"
				print geneFile2
				with open(geneFile2, 'a') as f:
					f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+"\n")
					f.write((alleleStr) +"\n")
				resultsList.append('LNFN:-1')
				perfectMatchIdAllele.append('LNFN')
				printinfo(genomeFile,geneFile) 
				print "LNFN, contains N bases! \n"
			
			
			
			else:
				
				#removing gaps
					
				alleleStr = alleleStr.replace('-', '')
				#lenExtraThresh=int(biggestAllelelen*0.1)
				lenExtraThresh=50
			
				#print alleleStr
				# --- it might be needed to obtain the reverse complement of the allele string --- #
				if match.sbjct_start > match.sbjct_end:
					alleleStr = reverseComplement(alleleStr)
					
				#if alleleStr in geneDict:  #if best match without gaps is already defined, example: best match allele was already defined but without gaps it's equal to a NA added
				#	alleleNumber = geneDict[ alleleStr ]
					
					################################################
					# EXACT MATCH --- MATCH == GENE --- GENE FOUND #						
					################################################
				#	perfectMatchIdAllele.append("EXC2-"+str(alleleNumber))
				#	resultsList.append('EXC2:' + str(alleleNumber) )
						

				#else: #check if best match without gaps are contained inside an already defined allele

				isContainedDefinedAllele = False	
				#print geneDict.keys()[0]
				definedAllele=''
				definedAlleleName=''
				for k in geneDict.keys():
					if alleleStr in k:
						definedAllele=k
						#print alleleStr
						isContainedDefinedAllele = True
						definedAlleleName=geneDict.get(k)
						break
						
				if isContainedDefinedAllele  and int(len(match.query))<=int(len(definedAllele))+lenExtraThresh and int(len(match.query))>=int(len(definedAllele))-lenExtraThresh :
					#allele without gaps is contained in a defined allele
					#best match with gaps has same size +1/-1 base as the defined allele
					
					#print int(len(definedAllele)), int(len(match.sbjct))
					
						
					if int(len(alleleStr))==int(len(definedAllele)): # if match without gaps has same size as the defined allele 
						tagAux = 'NA1:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA1-"+str(alleleI))
						
					elif int(len(alleleStr))==int(len(definedAllele))-1 : # if match without gaps has minus one base than the defined allele
						
						tagAux = 'NA2:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA2-"+str(alleleI))
					#elif int(len(alleleStr))==int(len(definedAllele))+1 : # if match without gaps has plus one base than the defined allele 
					#	tagAux = 'NA3:'
					#	printinfo(genomeFile,geneFile) 
					#	perfectMatchIdAllele.append("NA3-"+str(alleleI))
						
					else:												# if match without gaps has more than one base missing comparing to the defined allele 
						tagAux = 'NA4:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA4-"+str(alleleI))
					#TODO catch +1 and others
							
					print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
					geneDict[alleleStr] = alleleI
						
					resultsList.append( tagAux + str(alleleI) )
						
					orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)))									
					# --- add the new allele to the gene fasta --- #
						
					fG = open( geneFile, 'a' )
					fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)) + '\n')
					fG.write( alleleStr + '\n')
					fG.close()
					alleleI += 1
					Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
					
				#if best match is not contained in an already defined allele, check if it has similar size with the match allele and has 0.8 similarity
				elif not isContainedDefinedAllele and idPercent > 0.8 and int(len(match.query))<=int(geneLen)+lenExtraThresh and int(len(match.query))>=int(geneLen)-lenExtraThresh :
					#best match with gaps has 80% identity
					#best match with gaps is the same size or +1/-1 as the defined allele
					
					ratio=float(len(alleleStr)) / float(geneLen)
					
					if ratio>=0.8 and ratio<=1.2: # if match without gaps has same size as the best match allele and 80%similarity
						
						tagAux = ''
						extraleft=0
						extraright=0
						tS=0
						tE=0
						#print int(geneLen), len(match.sbjct)
						#print match.sbjct
						#print match
						handle = open(genomeFile, "rU")
						record_dict = SeqIO.to_dict(SeqIO.parse(handle, "fasta"))
						handle.close()
						record= record_dict[bestMatchContig]
						#print match.sbjct
						#if(int(len(alleleStr))<int(len(match.query)) and int(len(match.query))<int(geneLen)) and int(geneLen)==int(match.sbjct_start): #if best match allele has missing bases, the tips would be cut
						
						#if len(match.sbjct)<geneLen and "-" not in match.sbjct:  #if the allele is not fully used against the match, compensate the tips
						
						if (1<int(match.sbjct_start) and 1<int(match.sbjct_end)):
							
							if match.sbjct_start > match.sbjct_end:
								extraleft=match.sbjct_end-1
								
							else:
								extraleft=match.sbjct_start-1
								
								
						if (int(geneLen)>int(match.sbjct_start) and int(geneLen)>int(match.sbjct_end) ): # if 3' tip bases of the allele are missing on the match
							
							
							if match.sbjct_start > match.sbjct_end:
								extraright=geneLen-match.sbjct_start
								
							else:
								extraright=geneLen-match.sbjct_end
								
						#print 	extraleft, 	extraright
						
						
						if match.sbjct_start > match.sbjct_end:
							tS=match.query_start-extraright-1
							tE=match.query_end+extraleft
							alleleStr=str(record.seq[tS:tE])
							alleleStr = reverseComplement(alleleStr)
						else:
							tS=match.query_start-extraleft-1
							tE=match.query_end+extraright
							alleleStr=str(record.seq[tS:tE])
						
						tagAux = 'NA5:'
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("NA5-"+str(alleleI))
							
						print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
						geneDict[alleleStr] = alleleI
							
						resultsList.append( tagAux + str(alleleI) )
							
						orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)))									
						# --- add the new allele to the gene fasta --- #
							
						fG = open( geneFile, 'a' )
						fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] +'_' + str(os.path.basename(genomeFile)) + '\n')
						fG.write( alleleStr + '\n')
						fG.close()
						alleleI += 1
						Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
						
					else:
						##################
						# LNF WTFFF #
						##################
						geneFile2= os.path.splitext(geneFile)[0] + "LNF3.fasta"
						print geneFile2
						with open(geneFile2, 'a') as f:
							f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
							f.write((alleleStr) +"\n")
							f.write(">Allele\n")
							f.write((bmAllele)+"\n")
						resultsList.append('LNF3')
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append("LNF3")
						print "No allele found"
						
				elif isContainedDefinedAllele:
							####################
						# UNDEFINED ALLELE #		# it is contained in another allele
						####################
						
					alleleStr=match.query
					#if match.sbjct_start > match.sbjct_end:    #### - error
						#alleleStr = reverseComplement(alleleStr)
					resultsList.append('UND:-1')
					perfectMatchIdAllele.append("undefined allele")
					printinfo(genomeFile,geneFile) 
					print "Undefined allele \n"
					
					geneFile2= os.path.splitext(geneFile)[0] + "undefined.fasta"
					print geneFile2
					with open(geneFile2, 'a') as f:
						f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
						f.write((alleleStr) +"\n")
						#f.write(">BlastBestMatch"+str(definedAlleleName)+"\n")
						#f.write((alleleStr)+"\n")
						f.write(">Allele"+str(definedAlleleName)+"\n")
						f.write((definedAllele)+"\n")
					
				
							
				elif lenRatio < 0.5:
						
					###############
					# SMALL MATCH #
					###############
								
					resultsList.append('SAC:-1')		# don't know what 'SAC' stands for
					perfectMatchIdAllele.append('small match')
					printinfo(genomeFile,geneFile) 
					print "lower than 50% match \n"	
							
				elif lenRatio < 0.8 and idPercent < 0.5:
						#####################
					# INCOMPLETE ALLELE #		# it was not possible to extend it to at least 80% of the length of the gene
					#####################
					resultsList.append('INC:-1')
					perfectMatchIdAllele.append('allele incomplete')
					printinfo(genomeFile,geneFile)
					print "Incomplete allele\n"
						
				else:	
					##################
					# LNF WTFFF #
					##################
					geneFile2= os.path.splitext(geneFile)[0] + "LNF2.fasta"
					print geneFile2
					with open(geneFile2, 'a') as f:
						f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
						f.write((alleleStr) +"\n")
						f.write(">Allele\n")
						f.write((bmAllele)+"\n")
					resultsList.append('LNF2')
					printinfo(genomeFile,geneFile) 
					perfectMatchIdAllele.append("LNF2")
					print "No allele found"

						
							
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

	listOfGenomes = []
	listOfGenomesDict = []

	fp = open(genomeFiles, 'r')

	for genomeFile in fp:

		genomeFile = genomeFile.rstrip('\n')
		listOfGenomes.append( genomeFile )
		genomeDict = {}
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:

			genomeDict[ contig.name ] = contig.seq

		listOfGenomesDict.append( genomeDict )

	fp.close()


	
	# ----------------------------- #
	# Each gene has a different job #
	# ----------------------------- #
	
	argumentsList = []
	lGenesFiles = []
	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

	
	gene_fp = open( genes, 'r')
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		lGenesFiles.append( gene )
		argumentsList.append( [gene, listOfGenomes, listOfGenomesDict] )

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
				statsaux=[0]*7 # EXC NA undef LNF LOT incomplete SAC 
				finalphylovinput+= "\n" + currentGenome + "\t"
				for gene in phylovout:
					
					val= str(gene[genome])
					finalphylovinput+= val + "\t"
					if "NA" in val:
						statsaux[1]+=1
					elif "undefined" in val:
						statsaux[2]+=1
					elif "LNF" in val:
						statsaux[3]+=1
					elif "LOT" in val:
						statsaux[4]+=1
					elif "incomplete" in val:
						statsaux[5]+=1 
					elif "small" in val:
						statsaux[6]+=1
					else:
						statsaux[0]+=1
					
				genome+=1
				statistics.append(statsaux)
			gOutFile = os.path.dirname( "./")
			gOutFile += args.o
			statswrite='Stats:\tEXC\tNA\tundefined\tLNF\tLOT\tincomplete\tsmall'
			i=0
			genome=0
			while genome<len(listOfGenomes):
				currentGenome = os.path.basename(listOfGenomes[genome])
				statsaux=[0]*7 # EXC NA undef LNF LOT incomplete SAC
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
