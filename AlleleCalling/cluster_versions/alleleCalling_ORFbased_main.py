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
import glob
import drmaa
import pickle
import shutil

# ================================================ MAIN ================================================ #

def main():
	
	#time ./alleleCalling_ORFbased_main.py -i asd.txt -g allffn.txt -o out_all_fnn_spades.txt -p True
	
	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='List of genes (fasta)', required=True)
	parser.add_argument('-o', nargs='?', type=str, help="Name of the output files", required=True)
	parser.add_argument('-p', nargs='?', type=str, help="True to give a phyloviz output file type, false is predefined", required=False)

	args = parser.parse_args()
	
	genomeFiles = args.i
	genes = args.g
	phylovizinput=True
	if(args.p):
		phylovizinput=args.p
	
	
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

	fp.close()
	
	gene_fp = open( genes, 'r')

	genepath=''
	basepath=''
	lGenesFiles = []
	argumentsList = []
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		lGenesFiles.append( gene )
		genepath=os.path.dirname(gene)
		basepath=os.path.join(genepath, "temp")
		if not os.path.exists(basepath):
			os.makedirs(basepath)
		filepath=os.path.join(basepath,str(os.path.basename(gene))+"_argList.txt")
		with open(filepath, 'wb') as f:
			var = [gene, listOfGenomes]
			pickle.dump(var, f)
		argumentsList.append(filepath)

#		callAlleles([gene, listOfGenomes, listOfCDSDicts, listOfGenomesDict])
		

	gene_fp.close()
	
	# ------------------------------------------------- #
	#           RUN PRODIGAL OVER ALL GENOMES           #
	# ------------------------------------------------- #




	print ("Starting Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	#poolJobs = Pool()

	totgenomes= len(listOfGenomes)
	
	"""for genome in listOfGenomes:
		#print genome
		#listOfCDSDicts.append(runProdigal(genome))
		filepath=os.path.join(basepath,str(os.path.basename(genome))+"_ORF.txt")
		with open(filepath, 'wb') as f:
			var = runProdigal(genome)
			pickle.dump(var, f)"""
	joblist =[]
	with drmaa.Session() as s:
		for genome in listOfGenomes:
			#print('Creating job template')
			jt = s.createJobTemplate()
			#print os.path.join(os.getcwd(), 'callAlleles.py')
			jt.remoteCommand = os.path.join(os.getcwd(), 'runProdigal.py')
			#print argList
			jt.args = [str(genome),basepath]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(genome)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			#print('Cleaning up')
			s.deleteJobTemplate(jt)
		s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
		#for curjob in joblist:
		#	print 'Collecting job ' + curjob
		#	retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
		#	print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)
	
	
	
	print ("Finishing Prodigal at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	# ----------------------------- #
	# Each gene has a different job #
	# ----------------------------- #
	


	print
	print ("Starting Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	

	#output=callAlleles([gene, listOfGenomes, listOfCDSDicts, listOfGenomesDict])
	#print output
#	raise SystemExit

	totloci= len(argumentsList)
	joblist =[]
	with drmaa.Session() as s:
		for argList in argumentsList:
			#print('Creating job template')
			jt = s.createJobTemplate()
			#print os.path.join(os.getcwd(), 'callAlleles.py')
			jt.remoteCommand = os.path.join(os.getcwd(), 'callAlleles.py')
			#print argList
			jt.args = [str(argList),basepath]
			jt.joinFiles=True
			jt.nativeSpecification='-V'
			jobid = s.runJob(jt)
			joblist.append(jobid)
			with open("jobsid.txt","a") as f:
				f.write(str(argList)+"\n"+str(jobid))
			print('Your job has been submitted with ID %s' % jobid)

			#print('Cleaning up')
			s.deleteJobTemplate(jt)
		#s.synchronize(joblist, drmaa.Session.TIMEOUT_WAIT_FOREVER, True)
		for curjob in joblist:
			print 'Collecting job ' + curjob
			retval = s.wait(curjob, drmaa.Session.TIMEOUT_WAIT_FOREVER)
			print 'Job: ' + str(retval.jobId) + ' finished with status ' + str(retval.hasExited)
	

	print ("Finished Allele Calling at : "+time.strftime("%H:%M:%S-%d/%m/%Y"))
	
	output=[]
	for gene in lGenesFiles:
		filepath=os.path.join(basepath, os.path.basename(gene)+"_result.txt")
		with open(filepath,'rb') as f:
			var = pickle.load(f)
			output.append(var)
			
	shutil.rmtree(basepath)
	shutil.rmtree(os.path.join(os.path.dirname(gene), "blastdbs"))
	
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
					
					alleleschema.append(geneOut[1][gene])
					"""genename=(geneOut[1][gene]).split("_")

					if(len(genename)!=1):
						alleleschema.append(genename[1])
					else:
						alleleschema.append(genename[0])"""
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
	
