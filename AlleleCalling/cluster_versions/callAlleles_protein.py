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
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import runBlastParser
import time
import pickle
import shutil

def getBlastScoreRatios(genefile,basepath):
	
	gene_fp = HTSeq.FastaReader(genefile)
	alleleI=0
	allelescores=[]
	alleleProt=''
	alleleList=[]
	for allele in gene_fp: #new db for each allele to blast it against himself
		alleleI+=1
		genome=-1
		alleleList.append(allele.seq)
		translatedSequence,x,y=translateSeq(allele.seq)
		alleleProt+=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")
	#basepath="./blastdbs/temp"+str(os.path.basename(genefile))
	#if not os.path.exists(basepath):
	#	os.makedirs(basepath)
	proteinfastaPath=os.path.join(basepath,str(os.path.basename(genefile)+'_protein.fasta'))
	blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(genefile) + '.xml')
	with open(proteinfastaPath, "wb") as f:
		f.write(alleleProt)
	
	print ("Starting Blast alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	Gene_Blast_DB_name = Create_Blastdb( proteinfastaPath, 1, True )
	print proteinfastaPath
		# --- get BLAST score ratio --- #
	cline = NcbiblastpCommandline(query=proteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		#print cline
	allelescore=0
	
	print ("Parse bsr blast at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	blast_records = runBlastParser(cline,blast_out_file, alleleProt)
	
	print ("Blasted alleles at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
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
	return alleleI,allelescores,Gene_Blast_DB_name,alleleList
	
def reDogetBlastScoreRatios(genefile,basepath,alleleI,allelescores2,newGene_Blast_DB_name,alleleList2):
	
	gene_fp = HTSeq.FastaReader(genefile)
	#alleleI=0
	#allelescores=[]
	alleleProt=''
	#alleleList=[]
	"""for allele in gene_fp: #new db for each allele to blast it against himself
		print allele
		alleleI+=1
		genome=-1
		alleleList2.append(allele.seq)
		translatedSequence,x,y=translateSeq(allele.seq)
		print translatedSequence
		alleleProt+=">"+str(alleleI)+"\n"+str(translatedSequence+"\n")"""
	
	alleleI+=1
		
	proteinfastaPath=genefile
	print proteinfastaPath
	blast_out_file2 = os.path.join(basepath,'blastdbs/temp.xml')
	#with open(proteinfastaPath, "wb") as f:
	#	f.write(alleleProt)
	#Gene_Blast_DB_name = Create_Blastdb( proteinfastaPath, 1, True )
		# --- get BLAST score ratio --- #
	cline = NcbiblastpCommandline(query=proteinfastaPath, db=newGene_Blast_DB_name, evalue=0.001, out=blast_out_file2, outfmt=5)
		#print cline
	allelescore=0
	blast_records = runBlastParser(cline,blast_out_file2, proteinfastaPath)
	found =False
	for blast_record in blast_records:
		found=False 
		print blast_record
		#print blast_record.header
		#print blast_record.alignments[0]
		for alignment in blast_record.alignments:
			print alignment,alignment.hsps
			if found is False:
				#print blast_record.query, alignment.hit_def
				for match in alignment.hsps:
					print match
					
						#print "---------------------"
					if(int(alignment.hit_def)== int(blast_record.query)):
						#print match
						allelescores2.append(int(match.score))
						found=True
						break
			else:
				break
	#print allelescores2, alleleList2
	return alleleI,allelescores2,alleleList2

def reverseComplement(strDNA):

	basecomplement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
        strDNArevC = ''
        for l in strDNA:

        	strDNArevC += basecomplement[l]

        return strDNArevC[::-1]

def translateSeq(DNASeq):
	seq=DNASeq
	reversedSeq=False
	try:
		myseq= Seq(seq)
		#print myseq
		protseq=Seq.translate(myseq, table=11,cds=True)
	except:
		reversedSeq=True
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
				reversedSeq=False
				try:
					seq=seq[::-1]							
					seq=reverseComplement(seq)
					myseq= Seq(seq)
					#print myseq
					protseq=Seq.translate(myseq, table=11,cds=True)
				except:
					print "translated error"
					raise
	return protseq,seq,reversedSeq
	

# ======================================================== #
#            Allele calling and classification             #
# ======================================================== #
def main():
	print ("Starting script at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
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
	
	basepath=os.path.join(temppath,os.path.splitext(geneFile)[0])
	if not os.path.exists(basepath):
			os.makedirs(basepath)
	#print geneFile
	gene_fp = HTSeq.FastaReader(geneFile)
	alleleI = 0
	#inverted=False
	#orderedAlleleNames=[]
	resultsList = []
	i = 0
	perfectMatchIdAllele=[]
	perfectMatchIdAllele2=[]
	allelescores=[]
	
	print ("Getting BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	alleleI,allelescores,Gene_Blast_DB_name,alleleList=getBlastScoreRatios(geneFile,basepath)
	print ("Finished BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	genome=-1	
	
	genomeDict = {}
	print ("starting allele call at: "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	for genomeFile in genomesList:
		print genomeFile
		bestmatch=[0,0,False,'',0] #score, score ratio, perfectmatch, key name of the DNA sequence string, allele ID
		currentGenomeDict={}
		currentCDSDict={}
		#currentCDSDict = listOfCDSDicts[i]
		
		filepath=os.path.join(temppath,str(os.path.basename(genomeFile))+"_ORF_Protein.txt")
		with open(filepath,'rb') as f:
			currentCDSDict = pickle.load(f)
		
		g_fp = HTSeq.FastaReader( genomeFile )
		for contig in g_fp:
			sequence=str(contig.seq)
			genomeDict[ contig.name ] = sequence
		
		currentGenomeDict = genomeDict
		#print currentGenomeDict
		#alleleI = 0
		#alleleProt=''
		#for allele in gene_fp: #new db for each allele to blast it against himself
		#	alleleI+=1
		#	alleleProt+=">"+str(alleleI)+"\n"+str(translateSeq(allele.seq)+"\n")
		#basepath="./blastdbs/temp"+str(os.path.basename(geneFile))
		#if not os.path.exists(basepath):
		#	os.makedirs(basepath)
		#with open(basepath+'/protein.fasta', "wb") as f:
		#	f.write(alleleProt)
		#Gene_Blast_DB_name = Create_Blastdb( basepath+'/protein.fasta', 1, True )
		genome+=1
		listOfCDS=currentCDSDict
		genomeProteinfastaPath=os.path.join(temppath,str(os.path.basename(genomeFile)+'_Protein.fasta'))
		
		
		print ("Blasting alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		blast_out_file = os.path.join(basepath,"blastdbs/"+os.path.basename(geneFile)+ '_List.xml')

		#with open(basepath+'/proteinList.fasta', "wb") as f:
		#	f.write(protList)
		#Gene_Blast_DB_name = Create_Blastdb( './temp/proteinList.fasta', 1, True )
		cline = NcbiblastpCommandline(query=genomeProteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		#print cline
		#try:
		
		print ("Parse bsr blast at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		blast_records = runBlastParser(cline, blast_out_file, genomeProteinfastaPath)
		#except:
		#	cline = NcbiblastpCommandline(query=genomeProteinfastaPath, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
		
		print ("Blasted alleles on genome at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
		
		for blast_record in blast_records:
				
			
			for alignment in blast_record.alignments:
				
				for match in alignment.hsps:
					#print blast_record.query
					#print match
					#print alleleI, len(allelescores)
					scoreRatio=float(match.score)/float(allelescores[int(alignment.hit_def)-1])
					#print scoreRatio
					#print alignment.hit_def
					cdsStrName=blast_record.query
					DNAstr=listOfCDS[">"+cdsStrName]

					AlleleDNAstr=alleleList[int(alignment.hit_def)-1]
					compare=False
					if DNAstr==AlleleDNAstr is False:
						try:
							DNAstr=reverseComplement(DNAstr)
							if DNAstr==AlleleDNAstr is False:
								pass
							else:
								compare=True
						except:
							pass
					else:
						compare=True
						
					if "N" in DNAstr or "K" in DNAstr or "R" in DNAstr:
						pass
						
					elif(scoreRatio == 1 and bestmatch[2] is False and compare is True):
						bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alignment.hit_def),match,len(AlleleDNAstr),blast_record.query_letters]
						#print alignment
						#print match
					elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is True):
						bestmatch=[match.score,scoreRatio,True,cdsStrName,int(alignment.hit_def),match,len(AlleleDNAstr),blast_record.query_letters]
						#print match
					elif(scoreRatio == 1 and bestmatch[2] is False and compare is False):
						bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alignment.hit_def),match,len(AlleleDNAstr),blast_record.query_letters]
						#print alignment
						#print match
					
					elif(scoreRatio == 1 and match.score>bestmatch[0] and compare is False):
						bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alignment.hit_def),match,len(AlleleDNAstr),blast_record.query_letters]
						#print match
					elif(match.score>bestmatch[0] and scoreRatio>0.6 and scoreRatio>bestmatch[1] and bestmatch[2] is False):
						#print match.query
						#print match.sbjct
						#print allelescores
						bestmatch=[match.score,scoreRatio,False,cdsStrName,int(alignment.hit_def),match,len(AlleleDNAstr),blast_record.query_letters]
						#print match
		#print bestmatch
		
		print ("Classifying the match at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))		
		if bestmatch[0]==0 or "N" in AlleleDNAstr or "K" in AlleleDNAstr or "R" in AlleleDNAstr :
					#if no best match was found
					
					###################
					# LOCUS NOT FOUND #
					###################
			if 	bestmatch[0]==0:		
				resultsList.append('LNF3:-1')            # append result to the list of results
				perfectMatchIdAllele.append('LNF')
				perfectMatchIdAllele2.append('LNF')
				#printinfo(genomeFile,geneFile)
				print "Locus not found, no matches \n"
			else:
				resultsList.append('LNFN:-1')            # append result to the list of results
				perfectMatchIdAllele.append('LNF')
				perfectMatchIdAllele2.append('LNF')
				#printinfo(genomeFile,geneFile)
				print "Locus has strange base (N, K or R) \n"
		
		
		elif bestmatch[2] is True:
			contigname=bestmatch[3]	
			
			contigname=contigname.split("&")
			matchLocation=contigname[2]	
			#matchLocation=matchLocation.split("-")
			contigname=contigname[0]	
			
			alleleStr=listOfCDS[">"+bestmatch[3]]
			protSeq,alleleStr,Reversed=translateSeq(alleleStr)
					#if a perfect match was found
					
						################################################
						# EXACT MATCH --- MATCH == GENE --- GENE FOUND #
						################################################
					
			perfectMatchIdAllele.append(str(bestmatch[4]))
			if not Reversed:
				perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"+")
			else:
				perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation)+"&"+"-")
			resultsList.append('EXC:' + str(bestmatch[4]) )
			
		else:
			
			match=bestmatch[5]
			#print match
			geneLen=bestmatch[6]

			contigname=bestmatch[3]	
			#print contigname
			
			contigname=contigname.split("&")
			matchLocation=contigname[2]	
			matchLocation=matchLocation.split("-")
			contigname=contigname[0]
			
			seq=currentGenomeDict[ contigname ]
			bestMatchContigLen=len(seq)
			
			alleleStr=listOfCDS[">"+bestmatch[3]]
			protSeq,alleleStr,Reversed=translateSeq(alleleStr)
			
			
			print match
			print matchLocation
			print bestMatchContigLen
			# get extra space to the right and left between the contig and match 
			rightmatchContig=bestMatchContigLen-int(matchLocation[1])	
			leftmatchContig=int(matchLocation[0])
			
			if Reversed:
				aux=rightmatchContig
				rightmatchContig=leftmatchContig
				leftmatchContig=aux
			"""else:
				rightmatchContig=bestMatchContigLen-int(matchLocation[0])	
				leftmatchContig=int(matchLocation[1])"""
			
			print rightmatchContig,leftmatchContig
			
			
			# get extra space to the right and left between the allele and match
			
			rightmatchAllele=geneLen-(int(match.sbjct_end)*3)	
			leftmatchAllele=(int(match.sbjct_start)*3)
			"""if Reversed: 
				aux=rightmatchAllele
				rightmatchAllele=leftmatchAllele
				leftmatchAllele=aux"""
				
			"""else:
				rightmatchAllele=geneLen-(int(match.sbjct_start	)*3)
				leftmatchAllele=(int(match.sbjct_end)*3)"""
			
			#if a best match was found but it's not an exact match	

					###########################
					# LOCUS ON THE CONTIG TIP #
					###########################
			
			print rightmatchAllele, leftmatchAllele
			print geneLen
			
			
			
			#if bestMatchContigLen <= geneLen:
			if leftmatchContig<leftmatchAllele and 	rightmatchContig < rightmatchAllele:
			
				resultsList.append('LOTSC:-1')
				perfectMatchIdAllele.append('LOTSC')
				perfectMatchIdAllele2.append('LOTSC')
				#printinfo(genomeFile,geneFile)
				#print match.query_start
				print "Locus is bigger than the contig \n"
			
			#if match.query_start ==1 and len(match.query) < geneLen:		
			elif leftmatchContig<leftmatchAllele:
				
				
				resultsList.append('LOT3:-1')
				perfectMatchIdAllele.append('LOT3')
				perfectMatchIdAllele2.append('LOT3')
				
				print "Locus is on the 3' tip of the contig \n"
			
			
			#elif match.query_end == bestMatchContigLen and len(match.query) < geneLen:
			elif 	rightmatchContig < rightmatchAllele:
				
				resultsList.append('LOT5:-1')
				perfectMatchIdAllele.append('LOT5')
				perfectMatchIdAllele2.append('LOT5')

				print "Locus is on the 5' tip of the contig \n"
			
			
						
				
		
				
			else:
						#######################
						# ADD INFERRED ALLELE #		# a new allele 
						#######################
						
												
					#print "infered allele has location : "+(CDSType)
					#printinfo(genomeFile,geneFile) 
				tagAux='INF'
				perfectMatchIdAllele.append( tagAux +"-"+str(alleleI+1))
				#perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1]))
				
				if not Reversed:
					perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"+")
				else:
					perfectMatchIdAllele2.append(str(contigname)+"&"+str(matchLocation[0])+"-"+str(matchLocation[1])+"&"+"-")
				
				
				print "New allele! Adding allele "+ tagAux + str(alleleI+1) +" to the database\n"
																					
				resultsList.append( tagAux + str(alleleI+1) )

						#orderedAlleleNames.append('allele_' + str(alleleI) + '_' + tagAux[:-1] +"_" +str(os.path.basename(genomeFile)))	
											# --- add the new allele to the gene fasta --- #
				
				
				appendAllele='>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n'
				fG = open( geneFile, 'a' )
				#fG.write('>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n')
				fG.write(appendAllele)
						#print alleleStr
					
				
				#print listOfCDS
				#alleleStr=listOfCDS[">"+bestmatch[3]]
				#match=bestmatch[5]
				#reverse the order if needed
				#if match.sbjct_start > match.sbjct_end: 
				#	alleleStr = reverseComplement(alleleStr)
				fG.write( alleleStr + '\n')
				fG.close()
				
				fG = open( os.path.join(basepath,str(os.path.basename(geneFile)+'_protein2.fasta')), 'w' )
				#fG.write('>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n')
				fG.write('>'+str(alleleI+1)+'\n'+str(protSeq) + '\n')
				fG.close()
				fG = open( os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta')), 'a' )
				#fG.write('>allele_' + str(alleleI+1) + '_' + tagAux[:-1] +"_" + str(os.path.basename(genomesList[genome])) + '\n')
				fG.write('>'+str(alleleI+1)+'\n'+str(protSeq) + '\n')
						#print alleleStr
				fG.close()	
				
				#print listOfCDS
				#alleleStr=listOfCDS[">"+bestmatch[3]]
				match=bestmatch[5]
				
						#alleleI += 1
						# --- remake blast DB --- #
				alleleList.append(alleleStr)
				Gene_Blast_DB_name = Create_Blastdb( os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta')), 1, True )
				print os.path.join(basepath,str(os.path.basename(geneFile)+'_protein.fasta'))
				genefile2= os.path.join(basepath,str(os.path.basename(geneFile)+'_protein2.fasta'))
				Gene_Blast_DB_name2 = Create_Blastdb( genefile2, 1, True )
				print ("Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
				alleleI,allelescores,alleleList=reDogetBlastScoreRatios(genefile2,basepath,alleleI,allelescores,Gene_Blast_DB_name2,alleleList)
				print allelescores
				print ("Done Re-calculating BSR at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	#x=y
	#shutil.rmtree(basepath)

	
	final =	(resultsList,perfectMatchIdAllele)	
	#return (resultsList)
	print ("Finished allele calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	filepath=os.path.join(temppath , os.path.basename(geneFile)+"_result.txt")
	filepath2=os.path.join(temppath , os.path.basename(geneFile)+"_result2.txt")
	#print filepath
	with open(filepath, 'wb') as f:
		pickle.dump(final, f)
	with open(filepath2, 'wb') as f:
		pickle.dump(perfectMatchIdAllele2, f)
	shutil.rmtree(basepath)
	return True
	
if __name__ == "__main__":
    main()
