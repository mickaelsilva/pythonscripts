import HTSeq
import argparse
import os.path
from Bio.Seq import Seq
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbiblastpCommandline
import collections


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
				except Exception as e:
					print "translated error"
					print e
					protseq=""
					#raise
	return protseq,seq,reversedSeq

def main():

	parser = argparse.ArgumentParser(description="Given an ffn file, recovers the genes that are not paralogs and have a size bigger than the g parameter provided")
	parser.add_argument('-i', nargs='?', type=str, help='ffn file', required=True)
	parser.add_argument('-g', nargs='?', type=int, help='int minimum size', required=True)
	
	args = parser.parse_args()
	genes = args.i
	sizethresh = args.g
	passSteps = False

	
	#translate to protein and create new file
	abspath=os.path.abspath(genes)
	filename=os.path.basename(genes)
	abspath=abspath.replace(filename,'')
	proteinfile=os.path.join(abspath,'proteins.fasta') 
	
	geneDict = {}
	protDict={}
	orderedprotDict=collections.OrderedDict()
	alreadyIn=[]
	totalgenes=0
	repeatedgenes=0
	smallgenes=0
	
	if not passSteps:
		print "not passing steps"
		with open(proteinfile, "wb") as f:
			g_fp = HTSeq.FastaReader( genes )
			totalgenes+=1
			for gene in g_fp:
				dnaseq=	str(gene.seq)
				protseq,x,y=translateSeq(dnaseq)
				if len(protseq)>1:
					
					if str(protseq) in alreadyIn:
						repeatedgenes+=1
						#print gene.name + " already saved "
					
					elif len(str(protseq))<67:
						smallgenes+=1
						
					else:	
						alreadyIn.append(str(protseq))
						protname=">"+str(gene.name)+"\n"
							
						#print protseq
							
						f.write(protname+str(protseq)+"\n")
						protDict[protname] = str(protseq)
						geneDict[str(gene.name)] = gene.seq
				else:

					print gene.name
			
			orderedprotList=[]
			orderedprotList=sorted(protDict.items(), key=lambda x: len(x[1]), reverse=True)
			
			
			i=0
			while i < len(orderedprotList):
				elem=orderedprotList[i]
				orderedprotDict[elem[0]] = elem[1]
				i+=1
				
			#print orderedprotDict
		print str(repeatedgenes) + " repeated genes out of "+ str(totalgenes)
		print str(smallgenes) + " small genes out of "+ str(totalgenes)
		print "protein file created"
				
		# first step -  remove genes contained in other genes or 100% equal genes
		
		# list of results - the output of the function
		resultsList = []
		
		auxDict={}
		g_fp = HTSeq.FastaReader( proteinfile )
		g=0
		j=0
		
		print "Checking if proteins are equal or substring of others..."
		
		# for each gene from all the annotated genes - starting with an empty dictionary, only add a new gene if the "to be added gene" is not contained or equal to a gene already added to the dictionary
		auxprot=[]
		#auxnames=auxDict.keys()
		for elem in orderedprotDict.items():
			#print len(elem[1])
			contained=False
			
			prot=str(elem[1])
			if any(prot in x for x in auxprot):
				g+=1
				contained=True
			
			else:
				auxDict[elem[1]] = elem[0]
				auxprot.append(str(elem[1]))
			
			print str(j)+ " out of " + str(len(orderedprotDict)	)
				
			j+=1
			#print "____" +str(j)
		print "%s genes are contained in other genes" %  (g)
		
		#overwrite the original file, obtaining a new file with unique genes
		
		with open(proteinfile, "wb") as f:
			allsequences=''
			for k,v in auxDict.iteritems():
				allsequences+=v+k+"\n"
			f.write(allsequences)
	
	else:
		
		totalgenes=0
		smallgenes=0
		g_fp = HTSeq.FastaReader( genes )
		totalgenes+=1
		for gene in g_fp:
			dnaseq=	str(gene.seq)
			protseq,x,y=translateSeq(dnaseq)
			if len(protseq)>1:
				
				if str(protseq) in alreadyIn:
					repeatedgenes+=1
					#print gene.name + " already saved "
				
				elif len(str(protseq))<67:
					smallgenes+=1
					
				else:	
					alreadyIn.append(str(protseq))
					protname=">"+str(gene.name)+"\n"
							
						#print protseq
							
					protDict[protname] = str(protseq)
					geneDict[str(gene.name)] = gene.seq
			else:

				print gene.name
	
	
	
	geneFile = os.path.abspath( proteinfile )
	print proteinfile
	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1, True )
	
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'
					# ------------------------------ RUNNING BLAST ------------------------------ #

	cline = NcbiblastpCommandline(query=geneFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline, blast_out_file, geneFile)
	toRemove=[]
	genesToKeep=[]
	log=["removed\tcause\texplanation"]
	for blast_record in blast_records:
		
		allelename=blast_record.query
		#if not blast_record.query == allelename:
		#	print str(blast_record.query), allelename
		allelename=allelename.split(" ")
		allelename=allelename[0]
		alleleLength=len(geneDict[allelename])
		#print alleleLength
		try:
			
			#if gene A is not on the toRemove list yet, add to genesToKeep list
			
			if str(blast_record.query) not in toRemove:
				genesToKeep.append(blast_record.query)
				
				i=0
				#if first alignement is not against self, gene B is bigger than gene A and very simillar - remove gene A from genesToKeep and add gene B instead
				if  not str(blast_record.query) == str((blast_record.alignments[0]).hit_def):
					genesToKeep.remove(str(blast_record.query))
					toRemove.append(str(blast_record.query))
					log.append(str(blast_record.query)+"\t"+str((blast_record.alignments[0]).hit_def)+"\t"+"2 is first best match")
					
					#if gene B is not on the toRemove list, add to genesToKeep list
					if str((blast_record.alignments[0]).hit_def) not in toRemove:
						genesToKeep.append(str((blast_record.alignments[0]).hit_def))
					#print blast_record.query
					#print str((blast_record.alignments[0]).hit_def)
					#print "first allignment is not the self alignment"
					raise
				
				selfblastscore=(((blast_record.alignments[0]).hsps)[0]).score


				
				
				while i<len(blast_record.alignments):
					align=blast_record.alignments[i]
					
					match=(align.hsps)[0]
					scoreRatio=float(match.score)/float(selfblastscore)
					
					
					#scoreRatio = float(match.score) / float(alleleLength)
					alleleLength2=len(geneDict[str(align.hit_def)])
					
					#if good match and gene B not in toremove list
					if(scoreRatio>0.6 and not str(align.hit_def) == str(blast_record.query) and str(align.hit_def) not in toRemove):
						
						#if gene B is bigger than gene A, keep bigger gene B
						if alleleLength2>alleleLength :
							genesToKeep.append(str(align.hit_def))
							genesToKeep.remove(str(blast_record.query))
							toRemove.append(str(blast_record.query))
							log.append(str(blast_record.query)+"\t"+str(align.hit_def)+"\t"+"2 is bigger and bsr >0.6")
							
							#print "gene "+str(align.hit_def)+" is bigger than gene "+str(blast_record.query)
							#print str(alleleLength2), str(alleleLength), scoreRatio
							raise
						#else add gene B to toremove list
						elif str(align.hit_def) in genesToKeep:
							genesToKeep.remove(str(align.hit_def))
							toRemove.append(str(align.hit_def))
							log.append(str(align.hit_def)+"\t"+str(blast_record.query)+"\t"+"2 is bigger and bsr >0.6")
							
					i+=1
			
			#else gene A is on toRemove list, add all similar genes (not in genesToKeep) list to the toRemove list
			else:		
						
				i=0
				selfblastscore=0
				for align in blast_record.alignments:
					if not (str(align.hit_def) == str(blast_record.query)):
						selfblastscore=((align.hsps)[0]).score
						print "gene "+str(align.hit_def)+" is bigger than gene "+str(blast_record.query)
						raise
				
				while i<len(blast_record.alignments):
					align=blast_record.alignments[i]
					match=(align.hsps)[0]
					scoreRatio=float(match.score)/float(selfblastscore)
					
					if align.hit_def not in genesToKeep and not str(align.hit_def) == str(blast_record.query) and scoreRatio>0.6 :
						toRemove.append(align.hit_def)
						log.append(str(align.hit_def)+"\t"+str(blast_record.query)+"\t"+"2 was on the removed list and bsr >0.6")
							
					else:
						pass
						#notparalogs.remove(align.hit_def)
						#paralogs.append(align.hit_def)
						#paralogs.remove(alignment.hit_def)
					i+=1
			

		except Exception as e:
			#print e
			pass
	with open("logfile.txt", "wb") as f:
		for elem in log:
			
			f.write(str(elem)+"\n")
	
	"""with open("aasdasd.txt", "a") as f:
		for paralog in paralogs:
			name=paralog.split(" ")
			name=name[0].replace("|","_")
			f.write(str(name)+"\n")"""
	
	genesToKeep=list(set(genesToKeep))
	toRemove=list(set(toRemove))
	s = set(toRemove)
	notcommonToKeep= [x for x in genesToKeep if x not in s]
	print len(toRemove)
	print len(genesToKeep)
	print len(notcommonToKeep)
	#x=y
	pathfiles=os.path.dirname(geneFile)
	pathfiles=pathfiles+"/"
	#print pathfiles
	#print paralogs
	g_fp = HTSeq.FastaReader( genes )
	removedparalogs=0
	removedsize=0
	totalgenes=0
	rest=0
	#raise i
	for contig in g_fp:
		totalgenes+=1
		name = contig.name+" "+contig.descr
		name2= contig.name
		
		
		if name2 not in toRemove and name2 in genesToKeep:
			if int(len(contig.seq))>sizethresh:
				namefile=contig.name
				namefile=namefile.replace("|","_")
				with open(pathfiles+namefile+".fasta", "wb") as f:
					f.write(">1\n"+contig.seq+"\n")
				rest+=1	
			else:
				removedsize+=1
		else:
			#print name
			#toRemove.remove(name2)
			removedparalogs+=1
		
	#with open("aasdasd.txt", "a") as f:
		#for item in notparalogs:
		#for item in paralogs:
			#f.write(str(item)+"\n")
	print "%s genes are contained in other genes" %  (g)
	print "Removed %s same Locus genes" % str(removedparalogs)
	print "Removed %s because of size " % str(removedsize)
	print "%s Scheme genes " % str(rest)
	print "total genes:" + str(totalgenes)
	
	
if __name__ == "__main__":
	main()
