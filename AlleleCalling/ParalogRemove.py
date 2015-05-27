import HTSeq
import argparse
import os.path
from CommonFastaFunctions import Create_Blastdb
from CommonFastaFunctions import LoadAlelleFasta
from CommonFastaFunctions import LoadAlellicProfileGeneric
from CommonFastaFunctions import WriteFasta
from CommonFastaFunctions import runBlast
from CommonFastaFunctions import runBlastParser
from Bio.Blast.Applications import NcbiblastnCommandline


def main():

	parser = argparse.ArgumentParser(description="Given an ffn file, recovers the genes that are not paralogs and have a size bigger than the g parameter provided")
	parser.add_argument('-i', nargs='?', type=str, help='ffn file', required=True)
	parser.add_argument('-g', nargs='?', type=int, help='int minimum size', required=True)
	
	args = parser.parse_args()
	genes = args.i
	sizethresh = args.g
	
	gene_fp = HTSeq.FastaReader(genes)
	geneFile = os.path.abspath( genes )
	Gene_Blast_DB_name = Create_Blastdb( geneFile, 1 )
	
	geneF = os.path.splitext( geneFile )[0]
	blast_out_file = geneF + '.xml'

	# list of results - the output of the function
	resultsList = []

					# ------------------------------ RUNNING BLAST ------------------------------ #

	cline = NcbiblastnCommandline(query=geneFile, db=Gene_Blast_DB_name, evalue=0.001, out=blast_out_file, outfmt=5)
	blast_records = runBlastParser(cline, blast_out_file, geneFile)
	paralogs=[]
	for blast_record in blast_records:
		try:
			alignment=blast_record.alignments[1]
			paralogs.append( alignment.hit_def)

		except:
			continue
	
	pathfiles=os.path.dirname(geneFile)
	pathfiles=pathfiles+"/"
	print pathfiles
	
	g_fp = HTSeq.FastaReader( genes )
	removedparalogs=0
	removedsize=0
	for contig in g_fp:
		name = contig.name+" "+contig.descr
		if name not in paralogs:
			if int(len(contig.seq))>sizethresh:
				namefile=contig.name
				namefile=namefile.replace("|","_")
				with open(pathfiles+namefile+".fasta", "wb") as f:
					f.write(">1\n"+contig.seq+"\n")
			else:
				removedsize+=1
		else:
			print name
			removedparalogs+=1
	print "Removed %s paralog genes" % str(removedparalogs)
	print "Removed %s because of size :" % str(removedsize)
	
	
	
if __name__ == "__main__":
	main()
