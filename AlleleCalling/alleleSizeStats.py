import sys
import os
import csv
import numpy as np
from numpy import array
import argparse
import matplotlib.pyplot as plt
from matplotlib import colors
import HTSeq
from operator import itemgetter
from collections import OrderedDict

def main():

	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genome files (list of fasta files)', required=True)
	
	args = parser.parse_args()
	genes = args.i
	
	gene_fp = open( genes, 'r')

	statistics=[]
	
	genenumber=0
	conservedlengthgenes=[]
	notconservedlengthgenes=[]
	for gene in gene_fp:
		
		gene = gene.rstrip('\n')

		gene_fp2 = HTSeq.FastaReader(gene)
		maxsize=0
		minsize=99999999999
		sizesum=0
		allelenumber=0
		aux=[0,0]
		sizes=[]
		
		for allele in gene_fp2: 
			#print allele
			allelesize=len(allele.seq)
			sizes.append(allelesize)
			#print allelesize
			if 	allelesize>=maxsize or allelesize<=minsize:
				if allelesize>maxsize:
					maxsize=allelesize
				if allelesize<minsize:
					minsize=allelesize
			else:
				aux.append(allelesize)
			sizesum+=allelesize
			allelenumber+=1
		aux[0]=	minsize
		mean=float(sizesum)/allelenumber
		#print mean
		aux[1]=	maxsize
		#aux[2]=	maxsize
		
		i=0
		
		for size in sizes:
			
			if (not float(size)> mean*1.2 and not float(size)< mean*0.8):
				
				i+=1
		rate = i/float(allelenumber)
		if rate>=0.8 :
			#print i,allelenumber
			conservedlengthgenes.append(os.path.basename(gene))
		else:
			notconservedlengthgenes.append(os.path.basename(gene))
		#aux=[minsize,sizesum/allelenumber,maxsize]
		#print aux
		statistics.append(aux)	
		
		#print statistics
		#adasd
	sortbymedia=sorted(statistics, key=itemgetter(1))
	sortbymedia.reverse()
	
	print "there are "+ str(len(conservedlengthgenes))+ " genes that have their size conserved"
	
	gene_fp = open( "conservedgenes.txt", 'r')
	
	conservedgenes=[]
	for gene in gene_fp:
		
		gene = gene.rstrip('\n')
		conservedgenes.append(gene)
	#print 	conservedgenes
	#print conservedlengthgenes
	asd=list(set(conservedgenes) - set(conservedlengthgenes))
	print (asd)
	
	s = set(conservedgenes)
	temp3 = [x for x in conservedlengthgenes if x in s]
	
	#print temp3
	#print len(temp3)
	print "there are "+ str(len(temp3))+ " genes that have their size conserved and present on the phyloviz input out of " + str(len(conservedgenes))+" genes"
	dddsaa
	"""for item in sortbymedia:
		
		plt.plot(genenumber,item[0],'o',label="minimum",color='green')
		plt.plot(genenumber,item[1],'o',label="mean",color='red')
		plt.plot(genenumber,item[2],'o',label="maximum",color='black')
		del item[0:3]
		#print sortbymedia
		for elem in item:
		
			plt.plot(genenumber,elem,'o',color='yellow')"""
	plt.figure()	
	bp=plt.boxplot(sortbymedia,patch_artist=True)
	
	"""for item in sortbymedia:
		genenumber+=1
		plt.boxplot(item,genenumber)"""
		
	#n, bins, patches = plt.hist(statistics, histtype='bar',stacked=True)
	
	#n, bins, patches = plt.hist(statistics, 10, normed=1, histtype='bar',color=['crimson', 'burlywood', 'chartreuse'],label=['Crimson', 'Burlywood', 'Chartreuse'])			
	handles, labels = plt.gca().get_legend_handles_labels()

	#print labels
	handle_list, label_list = [], []
	by_label = OrderedDict(zip(labels, handles))
	
	frame1 = plt.gca()			
	frame1.axes.get_xaxis().set_visible(False)
	plt.ylabel('nucleotide length of the gene')
	#plt.legend(by_label.values(), by_label.keys(),numpoints=1 , loc='lower right')
	plt.title("allele size comparison per gene for "+str(len(sortbymedia))+"genes")
	plt.setp(bp['boxes'], color='black',facecolor="red")
	plt.setp(bp['whiskers'], color='blue')
	plt.setp(bp['fliers'], color='yellow', marker='+')
	#plt.xlim([-1,5])
	plt.show()			
	plt.close()

if __name__ == "__main__":
    main()
