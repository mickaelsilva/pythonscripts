import os
import argparse
import matplotlib.pyplot as plt
from matplotlib import colors
import HTSeq
from operator import itemgetter
from collections import OrderedDict
from copy import deepcopy
from matplotlib.ticker import ScalarFormatter
import mpld3
from mpld3 import utils, plugins
import numpy

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

class ClickInfo(plugins.PluginBase):
    """Plugin for getting info on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clickinfo", ClickInfo);
    ClickInfo.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo.prototype.constructor = ClickInfo;
    ClickInfo.prototype.requiredProps = ["id"];
    ClickInfo.prototype.defaultProps = {labels:null}
    function ClickInfo(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };
    
    ClickInfo.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        var labels = this.props.labels;
        obj.elements().on("mousedown",function(d, i){ 
                            window.open(labels, '_blank')});
    }
    """
    def __init__(self, points, labels):
		self.points = points
		self.labels = labels
		self.dict_ = {"type": "clickinfo","id": utils.get_id(points),"labels": labels}


def buildPlot(nparray,ReturnValues):
	
	try:
		plt.close('all')
	except:
		pass
	if not ReturnValues	:
		plt.figure()
	else:
		fig,ax = plt.subplots(figsize=(20,10))
		#bp=plt.boxplot(nparray,patch_artist=True)
		bp=plt.boxplot(nparray, whis=[5,95])
	
	
		
	handles, labels = plt.gca().get_legend_handles_labels()

	#print labels
	handle_list, label_list = [], []
	by_label = OrderedDict(zip(labels, handles))
	
	frame1 = plt.gca()			
	frame1.axes.get_xaxis().set_visible(False)
	frame1.axes.get_xaxis().set_ticks([])
	plt.ylabel('nucleotide length of the gene')
	plt.title("allele size comparison per gene for "+str(len(nparray))+"genes")
	
	#for box in bp['boxes']:
		# change outline color
	#	box.set( color='white', linewidth=2)
    # change fill color
	#	box.set( facecolor = 'black' )
	
	#plt.setp(bp['boxes'], color='blue')
	plt.setp(bp['whiskers'], color='blue')
	plt.setp(bp['fliers'],marker='o', markerfacecolor='green',linestyle='none')
	
	return plt,ax,fig,bp
	

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
	allsizes=[]
	
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
		#firstallelesize=len((gene_fp2[0]).seq)
		
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
		aux[1]=	maxsize

		
		i=0
		#print sizes
		moda=max(set(sizes), key=sizes.count)
		median=numpy.median(numpy.array(sizes))
		#print median
		modaStats.append(moda)

		for size in sizes:

			if (not float(size)> moda*(1+threshold) and not float(size)< moda*(1-threshold)):
				i+=1
		rate = i/float(allelenumber)

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


		#aux.append(gene)
		#aux.append(median)
		#aux.append(os.path.basename(gene))
		sizes.append(gene)
		sizes.append(sizes[0])
		sizes.append(median)
		
		#statistics.append(aux)	
		allsizes.append(sizes)


	print "\n"+str(z)+ " genes with only one allele\n"
	
	sortbymedia=sorted(allsizes, key=itemgetter(-1))
	sortbymedia.reverse()
	
	
	allsizes2=deepcopy(allsizes)
	
	sortbyfirstallele=[]
	sortbyfirstallele=sorted(allsizes2, key=itemgetter(-2))
	sortbyfirstallele.reverse()
	
	
	orderedlistgene=[]
	for elem in sortbymedia:
		elem.pop(-1)
		elem.pop(-1)
		orderedlistgene.append(elem.pop(-1))
		
	orderedlistgeneFirst=[]

	for elem in sortbyfirstallele:
		elem.pop(-1)
		elem.pop(-1)
		orderedlistgeneFirst.append(elem.pop(-1))	
			
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
	

		
	if not ReturnValues	:
		plt=buildPlot(sortbymedia,ReturnValues)
		
		plt.show()	
		plt.close('all')
	else:
		print "Creating the plot image files"
		
		genebasename=str(os.path.basename(genes))
		genebasename=genebasename.split(".")
		genebasename=genebasename[0]
		
		
		plt,ax,fig,bp=buildPlot(sortbymedia,ReturnValues)
		
		allboxes=bp.get('boxes')
		i=0
		for box in allboxes:
			mpld3.plugins.connect(fig, plugins.LineLabelTooltip(box,label=os.path.basename(orderedlistgene[i]),voffset=50, hoffset=10))
			mpld3.plugins.connect(fig, ClickInfo(box,(orderedlistgene[i])))
			i+=1
		allmedians=bp.get('medians')
		i=0
		for median in allmedians:
			mpld3.plugins.connect(fig, plugins.LineLabelTooltip(median,label=os.path.basename(orderedlistgene[i]),voffset=50, hoffset=10))
			mpld3.plugins.connect(fig, ClickInfo(median,(orderedlistgene[i])))
			i+=1
		
			
		#mpld3.show()

		"""boxplothtml=mpld3.fig_to_dict(fig)
		
		plt.close('all')
		
		plt,ax,fig,bp=buildPlot(sortbyfirstallele,ReturnValues)
		
		allboxes=bp.get('boxes')
		i=0
		for box in allboxes:
			mpld3.plugins.connect(fig, plugins.LineLabelTooltip(box,label=os.path.basename(orderedlistgeneFirst[i]),voffset=50, hoffset=10))
			mpld3.plugins.connect(fig, ClickInfo(box,(orderedlistgeneFirst[i])))
			i+=1
		allmedians=bp.get('medians')
		i=0
		for median in allmedians:
			mpld3.plugins.connect(fig, plugins.LineLabelTooltip(median,label=os.path.basename(orderedlistgeneFirst[i]),voffset=50, hoffset=10))
			mpld3.plugins.connect(fig, ClickInfo(median,(orderedlistgeneFirst[i])))
			i+=1
		
			
		mpld3.show()"""
		
		#plt.savefig(imagesDir+"plot.png", bbox_inches='tight')
		
		plt.close('all')
		

		fig, ax = plt.subplots(figsize=(25.5,15.0))
		bp=plt.hist(modaStats,100,rwidth=0.8)
		plt.ylabel('Number of occurrences')
		plt.xlabel('Allele Size')
		start, end = ax.get_xlim()
		ticks=range(int(start),int(end),500)
		plt.xticks(ticks)
		plt.grid(True)
		#ax.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.1f'))
		#plt.show()
		
		#histplothtml=mpld3.fig_to_dict(fig,template_type="simple")
		histplothtml=mpld3.fig_to_dict(fig)
		#mpld3.show()
		
		#plt.savefig(imagesDir+"histogram_mode.png", bbox_inches='tight',histtype='step')
		return notconservedlengthgenes,len(conservedgenes),z,boxplothtml,histplothtml

if __name__ == "__main__":
    main()
