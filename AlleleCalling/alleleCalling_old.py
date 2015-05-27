#!/usr/bin/python
import HTSeq
import sys
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

                                cdsDict[contigTag] = tempList
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

                cdsDict[contigTag] = tempList
	
	return cdsDict


def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]



def extendCDS(contigTag, cdsDict, matchStart, matchStop, contigsDict, maxsize):
	
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
	for j in listCdsInds:
		
              # cds start
		tS = j[0]
               # cds stop
		tE = j[1]

		cdslen=tE-tS
		#print tS , matchStart
		#print tE, matchStop
		#print cdslen
		lenratio=float(cdslen) / float(maxsize)
		#print lenratio
               # CDS found
		if(lenratio>0.8 and lenratio <1.2 and cdslen>len(cdsString)):
			#print cdsString
			
			if(cdslen==maxsize and matchStart==tS and tE == matchStop): # se do mesmo tamanho que o alelo e comeca e acaba no mesmo sitio que o match
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "same size as allele"
				break
			
			elif (len(cdsString)<(tE-tS) and matchStart >= tS and tE >= matchStop): # se maior que a cds anterior, comecar antes do match e acabar depois do match

				#if contigTag not in contigsDict:
					#print contigTag, 'not in contigsDict'
				#	raise SystemExit
				
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "larger than match"
			
				
			elif(len(cdsString)<(tE-tS) and matchStart < tS and tE >=matchStop and ((maxsize/2)> (tS-matchStart))): # se maior que a cds anterior, comecar depois do match e acabar depois do match
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "start codon inside match"
				
			elif(len(cdsString)<(tE-tS) and matchStart >= tS and tE <=matchStop and ((maxsize/2)> (matchStop-tE))): # se maior que a cds anterior, comecar antes do match e acabar antes do match
				
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "early stop codon in match"
				
				
			elif (len(cdsString)<(tE-tS) and matchStart < tS and tE <matchStop):# se maior que a cds anterior, comecar depois do match e acabar antes do match
				cdsString = contigsDict[ contigTag ][ tS:tE ].upper()
				CDStype = "start and end inside match"
				
		elif (matchStop< tS and matchStop< tE and cdsString !=''):
			return True, cdsString, CDStype
		elif (matchStop< tS and matchStop< tE and cdsString ==''):
			return False, cdsString, CDStype	
			
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
	for allele in gene_fp:
		if allele.seq in geneDict:
			print "\nWARNING: this file contains a repeated allele, it should be checked. Ignoring it now!\n", geneFile
		else:
			if len(allele.seq)>biggestAllelelen:
				biggestAllelelen=len(allele.seq)
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

					
					if (int(match.identities)==int(bmAlleleLen2) and int(match.identities)==int(len(match.sbjct)) and "N" not in match.query ): 
						index=orderedAlleleNames.index(alignment.hit_def)
						bmAlleleLen= len(geneDict.keys()[index])
						
						lenratio=float(len(match.query))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, alignment.hit_def,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						perfectMatch=True
						index=orderedAlleleNames.index(alignment.hit_def)
						bmAlleleLen= len(geneDict.keys()[index])
						break
						
					elif scoreRatio > bestMatch[2]:
						index=orderedAlleleNames.index(alignment.hit_def)
						bmAllele=geneDict.keys()[index]
						bmAlleleLen= len(bmAllele)
						lenratio=float(len(match.query))/float(bmAlleleLen)
						bestMatch = [blast_record.query, match, scoreRatio, alignment.hit_def,lenratio,bmAlleleLen]
						bestMatchContig=contigTag
						bestMatchContigLen=blast_record.query_letters
						if match.sbjct_start > match.sbjct_end:
							inverted=True
						#print match.query
						bestalignlen=alignment.length
						
						
					
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
			resultsList.append('LNF3:-1')            # append result to the list of results
			perfectMatchIdAllele.append('LNF')
			printinfo(genomeFile,geneFile)
			print "Locus not found \n"
			continue
		
		#TODO check identities >0.8
		
		if perfectMatch is True:
			#TODO perfect match to top
			if match.sbjct_start > match.sbjct_end:
				alleleStr = reverseComplement(alleleStr)
			#TODO test replace -
			#alleleStr = alleleStr.replace('-', '')
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
		
						###################
						# LOCUS NOT FOUND #
						###################
			
		#elif bestMatch[0] == '':
		#	resultsList.append('LNF:-1')            # append result to the list of results
		#	perfectMatchIdAllele.append('LNF')
		#	printinfo(genomeFile,geneFile)
		#	print "Locus not found \n"

		elif bestMatch[0] != '' and perfectMatch is not True:
						
				

					###########################
					# LOCUS ON THE CONTIG TIP #
					###########################
			
			
			#if match.query_start == 1 or bestMatchContigLen <= match.query_end:
			## TODO-
			## 1 -  LOT5 match.query_start ==1 and match.length < match.subj.length (allele length) alignement length
			## 2 - LOT 3' match.query_end == match.query.length (contig length) and match.length < contig length (allele length??)
			## 3 - LOT SC bestMatchContigLen <= allele length
			
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
				extended, strCDS, CDSType = extendCDS(bestMatchContig, currentCDSDict, match.query_start, match.query_end, currentGenomeDict, geneLen)
				# --- if it was possible to extend it using prodigal --- #
				
				#print extended
				#print strCDS
				#print CDSType
				

				if extended :
					alleleStr = strCDS
					
					lenRatio = float(len(strCDS)) / float(geneLen)
					#print alleleStr
					#print lenRatio
				elif not extended and biggestAllelelen > geneLen:
					extended, strCDS, CDSType = extendCDS(bestMatchContig, currentCDSDict, match.query_start, match.query_end, currentGenomeDict, biggestAllelelen)
					if extended :
						alleleStr = strCDS
					
						lenRatio = float(len(strCDS)) / float(geneLen)
					else:
						alleleStr = alleleStr.replace('-', '')
				
				
				else:
					# --- removing gaps '-' --- #
				#print alleleStr
					
					alleleStr = alleleStr.replace('-', '')

				# --- continuing the allele calling --- #

					
			
					#print geneDict
					#print alleleStr
					# --- it might be needed to obtain the reverse complement of the allele string --- #
				if match.sbjct_start > match.sbjct_end:
					alleleStr = reverseComplement(alleleStr)
					
				if alleleStr in geneDict:
					alleleNumber = geneDict[ alleleStr ]
					
					################################################
					# EXACT MATCH --- MATCH == GENE --- GENE FOUND #						
					################################################
					perfectMatchIdAllele.append(alleleNumber)
					resultsList.append('EXC:' + str(alleleNumber) )
						

				else:

					isUndefined = False	
					#print geneDict.keys()[0]
					defAllele=''
					defAlleleName=''
					for k in geneDict.keys():
						if alleleStr in k:
							defAllele=k
							#print alleleStr
							isUndefined = True
							defAlleleName=geneDict.get(k)
							break

						
					if extended and isUndefined and idPercent > 0.8 and ((int(len(match.query))==int(len(defAllele)) or int(len(match.query))==int(len(defAllele))+1 or int(len(match.query))==int(len(defAllele))-1)) :
						#extended allele to compare may be different from the allele to compare from bm	
						alleleStr=match.query
							
						alleleStr = alleleStr.replace('-', '')
							
						if match.sbjct_start > match.sbjct_end:    #### - error??
							alleleStr = reverseComplement(alleleStr)
							
						if int(len(alleleStr))==int(len(defAllele)): # se o match for do mesmo tamanho que o alello
							tagAux = 'NA1:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA1-"+str(alleleI))
							
						elif int(len(alleleStr))==int(len(defAllele))-1 : # se o match tiver uma base a mais que o alelo
							tagAux = 'NA2:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA2-"+str(alleleI))
							
						else:												#se o match tiver uma base a menos que o alelo
							tagAux = 'NA3:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA3-"+str(alleleI))
							
								
						print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
						geneDict[alleleStr] = alleleI
							
						resultsList.append( tagAux + str(alleleI) )
							
						orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)))									
						# --- add the new allele to the gene fasta --- #
							
						fG = open( geneFile, 'a' )
						fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)) + '\n')
						fG.write( alleleStr + '\n')
						fG.close()
						alleleI += 1
						Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
					
					elif not extended and idPercent > 0.8 and ((int(len(match.query))==int(geneLen) or int(len(match.query))==int(geneLen)+1 or int(len(match.query))==int(geneLen)-1)) :
							
						alleleStr=match.query
							
						alleleStr = alleleStr.replace('-', '')
							
						if match.sbjct_start > match.sbjct_end:    #### - error??
							alleleStr = reverseComplement(alleleStr)
							
						if int(len(alleleStr))==int(geneLen): # se o match for do mesmo tamanho que o alello
							tagAux = 'NA4:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA4-"+str(alleleI))
							
						elif int(len(alleleStr))==int(geneLen)-1 : # se o match tiver uma base a mais que o alelo
							tagAux = 'NA5:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA5-"+str(alleleI))
							
						else:												#se o match tiver uma base a menos que o alelo
							tagAux = 'NA6:'
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append("NA6-"+str(alleleI))
							
								
						print "New allele found! Adding allele "+ tagAux + str(alleleI) +" to the database"
						geneDict[alleleStr] = alleleI
							
						resultsList.append( tagAux + str(alleleI) )
							
						orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)))									
						# --- add the new allele to the gene fasta --- #
							
						fG = open( geneFile, 'a' )
						fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)) + '\n')
						fG.write( alleleStr + '\n')
						fG.close()
						alleleI += 1
						Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
					
							
					elif isUndefined:

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
							#f.write(">BlastBestMatch"+str(defAlleleName)+"\n")
							#f.write((alleleStr)+"\n")
							f.write(">Allele"+str(defAlleleName)+"\n")
							f.write((defAllele)+"\n")
						
					
						
						
						
					else:
						if not extended :
							
								
							if lenRatio < 0.5:
							
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
								print "Not extended and no allele found"

						else:

								#######################
								# ADD INFERRED ALLELE #		# a new allele that was extended with prodigal
								#######################
							if(CDSType=='larger than match'):
								tagAux = 'INF1:'
							elif(CDSType=='start codon inside match'):
								tagAux = 'INF2:'
							elif(CDSType=='early stop codon in match'):
								tagAux = 'INF3:'
							elif(CDSType=='same size as allele'):
								tagAux = 'INF4:'
							else:
								tagAux = 'INF5:'
								
							print "infered allele has location : "+(CDSType)
							printinfo(genomeFile,geneFile) 
							perfectMatchIdAllele.append( tagAux +"-"+str(alleleI))
							print "New allele Infered with prodigal! Adding allele "+ tagAux + str(alleleI) +" to the database\n"
							
							
								
							geneDict[alleleStr] = alleleI
								
							resultsList.append( tagAux + str(alleleI) )
								
							orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)))	
							# --- add the new allele to the gene fasta --- #

							fG = open( geneFile, 'a' )
							fG.write('>allele_' + str(alleleI) + '_' + tagAux[:-1] + str(os.path.basename(genomeFile)) + '\n')
							#print alleleStr
							fG.write( alleleStr + '\n')
							fG.close()
							alleleI += 1
							

							# --- remake blast DB --- #
							
							Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
							
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
				statsaux=[0]*8 # EXC NA INF undef LNF LOT incomplete SAC
				finalphylovinput+= "\n" + currentGenome + "\t"
				for gene in phylovout:
					
					val= str(gene[genome])
					finalphylovinput+= val + "\t"
					if "NA" in val:
						statsaux[1]+=1
					elif "INF" in val:
						statsaux[2]+=1
					elif "undefined" in val:
						statsaux[3]+=1
					elif "LNF" in val:
						statsaux[4]+=1
					elif "LOT" in val:
						statsaux[5]+=1
					elif "incomplete" in val:
						statsaux[6]+=1 
					elif "small" in val:
						statsaux[7]+=1
					else:
						statsaux[0]+=1
					
				genome+=1
				statistics.append(statsaux)
			gOutFile = os.path.dirname( "./")
			gOutFile += args.o
			statswrite='Stats:\tEXC\tNA\tINF\tundefined\tLNF\tLOT\tincomplete\tsmall'
			i=0
			genome=0
			while genome<len(listOfGenomes):
				currentGenome = os.path.basename(listOfGenomes[genome])
				statsaux=[0]*8 # EXC NA INF undef LNF LOT incomplete SAC
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
