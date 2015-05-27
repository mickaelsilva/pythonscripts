import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import colors
import HTSeq
from operator import itemgetter
from collections import OrderedDict
from copy import deepcopy
from matplotlib.ticker import ScalarFormatter

def main():

	parser = argparse.ArgumentParser(description="This program screens a set of genes in a fasta file.")
	parser.add_argument('-i', nargs='?', type=str, help='List of genes files (list of fasta files)', required=True)
	parser.add_argument('-t', nargs='?', type=float, help='threshold', required=False)
	parser.add_argument('-c', nargs='?', type=bool, help='All alleles must be conserved', required=False)
	parser.add_argument('-r', nargs='?', type=bool, help='Return values', required=False)
	
	args = parser.parse_args()
	genes = args.i
	

	try:
		threshold=float(args.t)
	except:
		threshold=0.05
		pass
	try:
		OneNotConserved=bool(args.c)
	except:
		OneNotConserved=False
		pass
	
	try:
		ReturnValues=bool(args.r)
	except:
		ReturnValues=False
		pass
	
	return getStats(genes,threshold,OneNotConserved,ReturnValues)

def buildPlot(nparray,ReturnValues):
	
	try:
		plt.close('all')
	except:
		pass
	if not ReturnValues	:
		plt.figure()
	else:
		fig,ax = plt.subplots(figsize=(25.5,15.0))
		#bp=plt.boxplot(nparray,patch_artist=True)
		bp=plt.boxplot(nparray,patch_artist=True, whis=[5,95])

		
	handles, labels = plt.gca().get_legend_handles_labels()

	#print labels
	handle_list, label_list = [], []
	by_label = OrderedDict(zip(labels, handles))
	
	frame1 = plt.gca()			
	frame1.axes.get_xaxis().set_visible(False)
	plt.ylabel('nucleotide length of the gene')
	plt.title("allele size comparison per gene for "+str(len(nparray))+"genes")
	plt.setp(bp['boxes'], color='black',facecolor="white")
	plt.setp(bp['whiskers'], color='blue')
	plt.setp(bp['fliers'],marker='o', markerfacecolor='green',linestyle='none')
	
	return plt,ax
	

def getStats(genes,threshold,OneNotConserved,ReturnValues):
	
	gene_fp = open( genes, 'r')

	statistics=[]
	
	genenumber=0
	conservedlengthgenes=[]
	notconservedlengthgenes=[]
	z=0
	conservedgenes=[]
	print "Genes with only 1 allele:"
	modaStats=[]
	for gene in gene_fp:
		
		gene = gene.rstrip('\n')
		conservedgenes.append(gene)
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
		#if maxsize>2200 and minsize<1600:
			#print gene
		aux[0]=	minsize
		mean=float(sizesum)/allelenumber
		#print mean
		aux[1]=	maxsize
		#aux[2]=	str(gene)
		
		i=0
		
		moda=max(set(sizes), key=sizes.count)
		modaStats.append(moda)
		#if "Unique_Acinetobacter_baumannii_AB0057.1.peg.gi_213157425_ref_YP_002319470.1_.fasta" in gene:
			#print moda
			#print sizes
		for size in sizes:
			#if "Unique_Acinetobacter_baumannii_AB0057.1.peg.gi_213157425_ref_YP_002319470.1_.fasta" in gene:
				#print size, moda*1.2, moda*0.8
			#if (not float(size)> mean*1.1 and not float(size)< mean*0.9):
			if (not float(size)> moda*(1+threshold) and not float(size)< moda*(1-threshold)):
				i+=1
		rate = i/float(allelenumber)
		#if "Unique_Acinetobacter_baumannii_AB0057.1.peg.gi_213157425_ref_YP_002319470.1_.fasta" in gene:
			#print rate,i
		if not OneNotConserved and (rate>=1 or (len(sizes)-i)<2) :
			if allelenumber==1:
				z+=1
				
				print os.path.basename(gene)
			conservedlengthgenes.append(os.path.basename(gene))
		elif OneNotConserved and rate>=1  :
			#print i,allelenumber
			if allelenumber==1:
				z+=1
				
				print os.path.basename(gene)
			conservedlengthgenes.append(os.path.basename(gene))
		else:
			notconservedlengthgenes.append(os.path.basename(gene))
		#aux=[minsize,sizesum/allelenumber,maxsize]
		#print aux
		statistics.append(aux)	
		
		#print statistics
		#adasd
	print "\n"+str(z)+ " genes with only one allele\n"
	sortbymedia=sorted(statistics, key=itemgetter(1))
	sortbymedia.reverse()
	"""j=0
	while j<100:
		print (sortbymedia[j])
		j+=1"""
	
	print str(len(conservedlengthgenes)) +" conserved genes"
	print str(len(notconservedlengthgenes)) +" not conserved genes"
	
	print "\nBuilding box plot..."
	
	if not ReturnValues	:
		with open ("notconserved.txt","wb") as f:
			for gene in notconservedlengthgenes:
				f.write(str(gene)+"\n")
	
	
	asd=list(set(conservedgenes) - set(conservedlengthgenes))
	#print (asd)
	
	s = set(conservedgenes)
	temp3 = [x for x in conservedlengthgenes if x in s]
	
	#print temp3
	#print len(temp3)
	#print "there are "+ str(len(temp3))+ " genes that have their size conserved and present on the phyloviz input out of " + str(len(conservedgenes))+" genes"
	#dddsaa
	"""for item in sortbymedia:
		
		plt.plot(genenumber,item[0],'o',label="minimum",color='green')
		plt.plot(genenumber,item[1],'o',label="mean",color='red')
		plt.plot(genenumber,item[2],'o',label="maximum",color='black')
		del item[0:3]
		#print sortbymedia
		for elem in item:
		
			plt.plot(genenumber,elem,'o',color='yellow')"""
	
	#plt.xlim([-1,5])	
		
	if not ReturnValues	:
		plt=buildPlot(sortbymedia,ReturnValues)
		
		plt.show()	
		plt.close('all')
	else:
		print "Creating the plot image files"
		
		genebasename=str(os.path.basename(genes))
		genebasename=genebasename.split(".")
		genebasename=genebasename[0]
		imagesDir="./resultsHTML/"+genebasename+"/"
		
		
		
		copySortbymedia = deepcopy(sortbymedia)
		j=0
		while len(copySortbymedia)>0:
			i=0
			
			batch=[]
			while i<100:
				try:
					batch.append(copySortbymedia.pop(0))
					i+=1
				except:
					i+=200
			plt,ax=buildPlot(batch,ReturnValues)
			j+=1
			"""plt.yscale('log', basey=10)
			for axis in [ax.xaxis, ax.yaxis]:
				axis.set_major_formatter(ScalarFormatter())"""
			plt.savefig(imagesDir+"plot"+str(j)+".png", bbox_inches='tight')
			
			
		plt,ax=buildPlot(sortbymedia,ReturnValues)
		plt.yscale('log', basey=10)
		for axis in [ax.xaxis, ax.yaxis]:
			axis.set_major_formatter(ScalarFormatter())
		plt.savefig(imagesDir+"plot.png", bbox_inches='tight')
		
		plt.close('all')
		
		#fig,ax = plt.subplots(figsize=(25.5,15.0))
		"""print modaStats
		print len (modaStats)
		print type(modaStats[0])
		print type(modaStats)
		
		asd=set(modaStats)
		
		bp=plt.hist(modaStats,sorted(list(asd)),histtype='stepfilled')"""
		fig, ax = plt.subplots(figsize=(25.5,15.0))
		bp=plt.hist(modaStats,100,rwidth=0.8)
		plt.ylabel('Number of times')
		plt.xlabel('Allele Size')
		start, end = ax.get_xlim()
		print start, end
		ticks=range(int(start),int(end),500)
		plt.xticks(ticks)
		plt.grid(True)
		#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
		#plt.show()
		plt.savefig(imagesDir+"histogram_mode.png", bbox_inches='tight',histtype='step')
		return notconservedlengthgenes,len(conservedgenes),z

if __name__ == "__main__":
    main()
