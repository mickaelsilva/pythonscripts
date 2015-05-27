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
import time
import pickle

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
		
		#if(biggestlenInsideMatch<lenInsideMatch and ((maxratio<1+sizeratio and maxratio>=1) or (minratio>1-sizeratio and minratio<=1))):
		if(biggestlenInsideMatch<lenInsideMatch and (maxratio<1+sizeratio and minratio>1-sizeratio)):
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
def main():
	
	try:
		input_file = sys.argv[1]
		temppath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"

	argumentList=[]
	with open(input_file,'rb') as f:
		argumentList = pickle.load(f)
	
	geneFile = argumentList[0]
	genomesList = argumentList[1]
	#listOfCDSDicts = argumentList[2]
	
	basepath=temppath
	
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
	geneF = os.path.basename(geneFile)
	blast_out_file = os.path.dirname(geneFile)+"/blastdbs/"+geneF + '.xml'

	# list of results - the output of the function
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	genomeDict = {}
	for genomeFile in genomesList:
		#currentCDSDict = listOfCDSDicts[i]
		
		filepath=os.path.join(basepath,str(os.path.basename(genomeFile))+"_ORF.txt")
		with open(filepath,'rb') as f:
			currentCDSDict = pickle.load(f)
		
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			genomeDict[ contig.name ] = sequence
		
		currentGenomeDict = genomeDict
		
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

		print ("Finished Blast at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))

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
					if (int(match.identities)==int(bmAlleleLen2) and int(match.identities)==int(len(match.sbjct)) and "N" not in match.query and "K" not in match.query and "R" not in match.query ): 
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
		print ("Finished choosing best match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		
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
			
			
			elif match.query_end == bestMatchContigLen and len(match.query) < geneLen:
				
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
				print ("Finished extension at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				#print ORFFoundInMatch
				#print strCDS
				#print CDSType
				isContainedDefinedAllele = False
				
				try:	
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
				
				except:
					if ORFFoundInMatch :
						alleleStr = strCDS
					##################
					# LNF WTFFF #
					##################
					geneFile2= os.path.splitext(geneFile)[0] + "LNF99.fasta"
					print geneFile2
					with open(geneFile2, 'a') as f:
						f.write(">"+ (str(os.path.basename(genomeFile)))+"|"+(str(os.path.basename(geneFile)))+" | "+str(bestMatchContig)+"\n")
						f.write((alleleStr) +"\n")
						f.write(">Allele\n")
						f.write((bmAllele)+"\n")
					resultsList.append('LNF99')
					printinfo(genomeFile,geneFile) 
					perfectMatchIdAllele.append("LNF99")
					print "A problem occurred"
						
							
	final =	(resultsList,perfectMatchIdAllele)	
	#return (resultsList)
	print ("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(basepath , os.path.basename(geneFile)+"_result.txt")
	#print filepath
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	return True
	
if __name__ == "__main__":
    main()
