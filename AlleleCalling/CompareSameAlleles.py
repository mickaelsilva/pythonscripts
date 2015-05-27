import HTSeq
import argparse
import os.path
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
import shutil
import csv

def main():

	#run -i CompareSameAlleles.py -i ./listgenesCompare3.txt -g ./listgenesCompare1.txt

	parser = argparse.ArgumentParser(description="Given two list of genes of the same genes, checks if the alleles are equal")
	parser.add_argument('-i', nargs='?', type=str, help='1st list of genes files to compare', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='2nd list of genes files to compare', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='clean output from an allele calling', required=True)
	
	args = parser.parse_args()
	geneFiles1 = args.i
	geneFiles2 = args.g
	outputfile = args.o
	
	listgenesToCompare=[]
	
	gene_fp = open( geneFiles2, 'r')
	for gene in gene_fp:
		gene = gene.rstrip('\n')
		listgenesToCompare.append(gene)
	
	gene_fp.close()
	
	phylovinput=[]
	
	with open(outputfile, 'r') as f:
		reader = csv.reader(f, delimiter="\t")
		phylovinput = list(reader)
		
	gene_fp = open( geneFiles1, 'r')
	
	orderedlistgenes=[]
	finalResult={}
	totalalleles=0
	totalalleles2=0
	notfoundAlleles=0
	containedid=0
	containsid=0
	if not os.path.exists("./notEqualAlleles"):
		os.makedirs("./notEqualAlleles")
	for gene in gene_fp:
		orderedlistgenes.append(gene)
		geneDict={}
		geneDict2={}
		generesult=[0,0,0,0]
		gene = gene.rstrip('\n')
		
		genename=str(os.path.basename(gene))
		#print genename
		samegene=filter(lambda x: genename in x, listgenesToCompare)
		
		if samegene and "LNF" not in genename and "undefined" not in genename:
			#print gene
			gene_fp2 = HTSeq.FastaReader(gene)
			alleleI=0
			for allele in gene_fp2:
				geneDict[ allele.seq ] = allele.name
				alleleI += 1
			totalalleles+=alleleI
			#print alleleI
			#print samegene
			gene_fp2 = HTSeq.FastaReader(samegene[0])
			alleleI=0
			for allele in gene_fp2:
				geneDict2[ allele.seq ] = allele.name
				alleleI += 1
			totalalleles2+=alleleI
			a=0
			#print geneDict, geneDict2
			for key,alleleNam in geneDict.iteritems():
				notfound=True
				contained=False
				contains=False
				savedKey=''
				savedAlleleName=''
				mostsimilarseq=[0,'','']
				for key2,alleleName in geneDict2.iteritems():
					#print key
					if key == key2:
						generesult[0]+=1
						notfound=False
						contained=False
						contains=False
						break
					elif key in key2:						
						notfound=False
						contains=True
						savedKey=key2
						savedAlleleName=alleleName
									
					elif key2 in key and not contains:						
						notfound=False
						contained=True
						savedKey=key2
						savedAlleleName=alleleName
						
						
				if 	notfound:
					genomename=alleleNam.split("INF")
					genomename=genomename[1]
					genomename=genomename.split(".")
					genomename=genomename[0]
					genomename=genomename[1:]
					#print genomename
					
					geneindex=0
					
					for genename2 in phylovinput[0]:
						if genename2==genename:
							#print genename2,genename
							break
						geneindex+=1
					#print geneindex
							
					alleleid=0	
					try:	
						for item in phylovinput:
							#print item[0], genomename
							if item[0]==genomename:
								
								alleleid=item[geneindex]
								break
					except:
						alleleid="allele_not_present"	
					#print alleleid
					
					#geneDict.keys()[alleleid+1] allele indexes start at 0
					
					generesult[3]+=1
					notfoundAlleles+=1
					with open("./notEqualAlleles/notfound.fasta", 'a') as f:
						f.write(">"+genename+"_"+alleleNam+"\n"+key+"\n")
						try:
							f.write(">"+genename+"_Allele"+str(alleleid)+"\n"+geneDict.keys()[int(alleleid)-1]+"\n")
						except:
							f.write(">"+genename+"_Allele_"+str(alleleid)+"\n"+"\n")
					
				elif contained:
					generesult[2]+=1
					notfoundAlleles+=1
					containedid+=1
					with open("./notEqualAlleles/contained.fasta", 'a') as f:
						f.write(">"+genename+"_"+alleleNam+"\n"+key+"\n")
						f.write(">"+genename+"_"+savedAlleleName+"\n"+savedKey+"\n")
						
				elif contains:
					generesult[1]+=1
					containsid+=1
					with open("./notEqualAlleles/contains.fasta", 'a') as f:
						f.write(">"+genename+"_"+alleleNam+"\n"+key+"\n")
						f.write(">"+genename+"_"+savedAlleleName+"\n"+savedKey+"\n")
							
		#print gene, generesult
				
		finalResult[genename]=generesult
		#print 	generesult[3]
	gene_fp.close()
	#print finalResult
	
	finalResulttsv='Stats:\tEQUALS\tCONTAINS\tCONTAINED\tDIFFERENT'

	totals=[0,0,0,0]
	#print finalResult
	for key, val in finalResult.iteritems():
		
		if "LNF" not in key and "undefined" not in key: 
			finalResulttsv+= "\n" + key + "\t"
			#print type(key)
			value=0
			for item in val:
				totals[value]+=item
				finalResulttsv+= str(item) + "\t"
				value+=1
			
			
	print "%s alleles are equal, %s alleles contain other alleles from the second allele list, %s alleles are contained on alleles from the second list and %s are different" % (totals	[0],totals	[1],totals	[2],totals	[3])
	print "from a total of %s and %s alleles" % (totalalleles, totalalleles2)
	with open("sameAllelesResult.txt", 'w') as f:
				f.write(str(finalResulttsv))

if __name__ == "__main__":
	main()
