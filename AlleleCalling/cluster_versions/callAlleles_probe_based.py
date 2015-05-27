#!/usr/bin/python

import HTSeq
import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna
from Bio.Blast import NCBIXML
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline

import os
import os.path
import string
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
import time
import pickle

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
					if (int(match.identities)==int(bmAlleleLen2) and int(match.identities)==int(len(match.sbjct)) and "N" not in match.query and "K" not in match.query and "Y" not in match.query and "R" not in match.query ): 
						
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
						bestalignlen=alignment.length
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
					
					
				

			elif 'N' in alleleStr or "K" in alleleStr or "R" in alleleStr or "Y" in alleleStr:
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
				print "LNFN, contains strange (N,K,R) bases! \n"
			
			
			
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
					Gene_Blast_DB_name = Create_Blastdb( geneFile, 1,False )
					
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
						try:
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
							Gene_Blast_DB_name = Create_Blastdb( geneFile, 1,False )
						except:
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
	print ("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(basepath , os.path.basename(geneFile)+"_result.txt")
	#print filepath
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	return True
	
if __name__ == "__main__":
    main()
