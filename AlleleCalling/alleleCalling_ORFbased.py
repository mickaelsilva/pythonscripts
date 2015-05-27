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


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]



def extendCDS(contigTag, cdsDict, matchStart, matchStop, contigsDict, maxsize,minsize,sizeratio):
	
		# ---------------------- #
        # Trying to extending it #
        # ---------------------- #
	CDStype=''
	cdsString = ''

        # --- tests if the cds contains the gene --- #
	if contigTag not in cdsDict:

		print "\nCheck if this contig contains CDSs:", contigTag
		#print cdsDict.keys()

		return False, '', CDStype

	listCdsInds = cdsDict[ contigTag ]

	# --- test if it is needed to invert the sequence --- #

	
	if matchStart > matchStop:
		
		aux = matchStart
		matchStart = matchStop
		matchStop = aux
		
	cdsIntersectMatch=[]
	for j in listCdsInds:
		
              # cds start
		tS = j[0]
               # cds stop
		tE = j[1]

		cdslen=tE-tS
		#print tS , matchStart
		#print tE, matchStop
		#print cdslen
		if(matchStart>= tS and (matchStop<=tE or matchStart<tE)): #if CDS start is before match start and the CDS is either after or before the match end
			cdsIntersectMatch.append(j)
			
		elif (matchStart<= tS and tS<=matchStop): #if CDS start inside of the match (doesnt matter where it ends)
			cdsIntersectMatch.append(j)	
		
		elif (matchStop< tS and matchStop< tE): #if CDS start after and finish after the match ends cycle
			break
	#print cdsIntersectMatch
	biggestlenInsideMatch=0
	lenInsideMatch=0
	for j in cdsIntersectMatch:
	     # cds start
		tS = j[0]
         # cds stop
		tE = j[1]

		cdslen=tE-tS
		maxratio=float(cdslen) / float(maxsize)
		minratio=float(cdslen) / float(minsize)
		
		#get length of the CDS against the match
		if tE<=matchStop and tS<matchStart: #if CDS starts before and ends inside match
			lenInsideMatch=tE-matchStart
			#   ___bm____     
			#__________
			#   _______
			
		elif tE<=matchStop and tS>=matchStart: # if CDS starts inside and finishes inside the match
			lenInsideMatch=tE-tS
			#   ___bm____     
			#     ____
			
		elif tE>matchStop and tS<matchStart: # if CDS is bigger than match
			lenInsideMatch=matchStop-matchStart
			#   ___bm____    
			# _____________
			
		elif tE>=matchStop and tS>matchStart: # if CDS starts inside match and ends outside
			lenInsideMatch=matchStop-tS
			#   ___bm____    
			#     __________
			#     _______
			
		else : #same size
			lenInsideMatch=tE-matchStart
			#  ___bm____    
			#  _________
		
		if(biggestlenInsideMatch<lenInsideMatch and ((maxratio<1+sizeratio and maxratio>=1) or (minratio>1-sizeratio and minratio<=1))):
			biggestlenInsideMatch=lenInsideMatch
			if (matchStart==tS and tE == matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "same size as match"
			#  ___bm____    
			#  _________
			
			elif (matchStart>tS and tE == matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "stop codon in match end"
			#   ___bm____    
			# ___________
			elif (matchStart==tS and tE > matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "start codon in match beggining"
			#   ___bm____     
			#   ___________
			
			elif (matchStart==tS and tE > matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "bigger than match"
			#   ___bm____    
			# _____________
			
			elif (matchStart<=tS and tE < matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "cds inside match"
			#   ___bm____     
			#     ____
			
			elif (matchStart<=tS and tE > matchStop):
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "start codon inside match"
			#   ___bm____    
			#     __________
			
			else : # two cds in same match??
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "stop codon inside match"
			#   ___bm____     
			#__________
			
			
		
	if cdsString =='' :
		return False, cdsString, CDStype
	else:
		return True, cdsString, CDStype

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
	listOfCDSDicts = argumentList[2]
	listOfGenomesDict = argumentList[3]
	
	gene_fp = HTSeq.FastaReader(geneFile)
	geneDict = {}
	alleleI = 1
	inverted=False
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

	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1, False )
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'

	# list of results - the output of the function
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	for genomeFile in genomesList:
		currentCDSDict = listOfCDSDicts[i]
		currentGenomeDict = listOfGenomesDict[i]
		
		#print genomeFile
		#print resultsList
		#print geneDict
		#print orderedAlleleNames
		i+=1		# it has to be incremented here
		if genomeFile[-1] == '\n':
			genomeFile = genomeFile[:-1]

                # ------------------------------ RUNNING BLAST ------------------------------ #
		#print Gene_Blast_DB_name
		#cline = NcbiblastnCommandline(query=genomeFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		cline = NcbiblastnCommandline(query=genomeFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		#print cline
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
						for seq, alleleid in geneDict.iteritems():
							if alleleid == index+1:
								bmAllele=seq
								break
						bmAlleleLen= len(bmAllele)
						lenratio=float(len(match.query))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, alignment.hit_def,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						bestMatchContigLen=blast_record.query_letters
						if match.sbjct_start > match.sbjct_end:
							inverted=True
						#print match.query
						bestalignlen=alignment.length
						#print match
						#print bmAlleleLen, bestMatchContig
						
						
					
					if perfectMatch==True:
						break


		# ---------- ALLELE CALLING AFTER DETERMINING BEST MATCH ---------- #

		
		
		try:
			#print bestMatch[0]
			match = bestMatch[1]
			#print match.query
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
			
			resultsList.append('LNF3:-1')            # append result to the list of results
			perfectMatchIdAllele.append('LNF')
			printinfo(genomeFile,geneFile)
			print "Locus not found, no matches \n"
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
				# ------------------------------------------------------------------------------------------------------- #
				#                                                                                                         #
				#                                   USING PRODIGAL TO TRY TO EXTEND CDS                                   #
				#                                                                                                         #
				# ------------------------------------------------------------------------------------------------------- #
				
				CDSType=''
				sizeratio=0.2
				ORFFoundInMatch, strCDS, CDSType = extendCDS(bestMatchContig, currentCDSDict, match.query_start, match.query_end, currentGenomeDict, biggestAllelelen, smallestAllelelen,sizeratio)
				# --- if it was possible to extend it using prodigal --- #
				
				#print ORFFoundInMatch
				#print strCDS
				#print CDSType
				isContainedDefinedAllele = False
				
					
				if ORFFoundInMatch :
					alleleStr = strCDS
					if match.sbjct_start > match.sbjct_end: #reverse the order if needed
						alleleStr = reverseComplement(alleleStr)
					
					lenRatio = float(len(strCDS)) / float(geneLen)
					defAllele=[]
					if alleleStr in geneDict:  #if ORF found is already defined
						alleleNumber = geneDict[ alleleStr ]
						
						################################################
						# EXACT MATCH --- MATCH == GENE --- GENE FOUND #						
						################################################
						perfectMatchIdAllele.append(alleleNumber)
						resultsList.append('EXC2:' + str(alleleNumber) )


						
					else:
								#######################
								# ADD INFERRED ALLELE #		# a new allele that was extended with prodigal
								#######################
						if(CDSType=='stop codon in match end'):
							tagAux = 'INF1:'
						elif(CDSType=='start codon in match beggining'):
							tagAux = 'INF2:'
						elif(CDSType=='bigger than match'):
							tagAux = 'INF3:'
						elif(CDSType=='same size as match'):
							tagAux = 'INF4:'
						elif(CDSType=='cds inside match'):
							tagAux = 'INF5:'
						elif(CDSType=='start codon inside match'):
							tagAux = 'INF6:'
						else:
							tagAux = 'INF7:'
								
						print "infered allele has location : "+(CDSType)
						printinfo(genomeFile,geneFile) 
						perfectMatchIdAllele.append( tagAux +"-"+str(alleleI))
						print "New allele Infered with prodigal! Adding allele "+ tagAux + str(alleleI) +" to the database\n"
							
							
								
						geneDict[alleleStr] = alleleI
								
						resultsList.append( tagAux + str(alleleI) )

						orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] +"_" +str(os.path.basename(genomeFile)))	
							# --- add the new allele to the gene fasta --- #

						fG = open( geneFile, 'a' )
						fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomeFile)) + '\n')
						#print alleleStr
						fG.write( alleleStr + '\n')
						fG.close()
						alleleI += 1
							

							# --- remake blast DB --- #
							
						Gene_Blast_DB_name = Create_Blastdb( geneFile, 1,False )
						
						
						
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
					print "CDS not found"

						
							
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

			genomeDict[ contig.name ] = contig.seq

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
	
	argumentsList = []
	lGenesFiles = []
	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	gene_fp = open( genes, 'r')
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		lGenesFiles.append( gene )
		argumentsList.append( [gene, listOfGenomes, listOfCDSDicts, listOfGenomesDict] )

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
