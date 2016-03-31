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

class ClickInfo2(plugins.PluginBase):
    """Plugin for getting info on click"""
    
    JAVASCRIPT = """
    mpld3.register_plugin("clickinfo", ClickInfo2);
    ClickInfo2.prototype = Object.create(mpld3.Plugin.prototype);
    ClickInfo2.prototype.constructor = ClickInfo2;
    ClickInfo2.prototype.requiredProps = ["id"];
    ClickInfo2.prototype.defaultProps = {labels:null}
    function ClickInfo2(fig, props){
        mpld3.Plugin.call(this, fig, props);
    };

    ClickInfo2.prototype.draw = function(){
        var obj = mpld3.get_element(this.props.id);
        labels = this.props.labels;
        obj.elements().on("mousedown",
                          function(d, i){ 
                            window.open(labels[i], '_blank')});
    }
    """
    def __init__(self, points, labels):
        self.points = points
        self.labels = labels
        """if isinstance(points, matplotlib.lines.Line2D):
            suffix = "pts"
        else:
            suffix = None"""
        self.dict_ = {"type": "clickinfo","id": utils.get_id(points, "pts"),"labels": labels}


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
	tabStats=[]
	allNumberAlleles=[]

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
		


		#per gene get all sizes, minimin size, maximum size, media and mode
		for allele in gene_fp2: 
			
			allelesize=len(allele.seq)
			sizes.append(allelesize)
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
		
		sizesnpList=numpy.array(sizes)
		tabStats.append([gene,numpy.amin(sizesnpList),numpy.amax(sizesnpList),numpy.mean(sizesnpList),numpy.std(sizesnpList)])
		
		i=0

		moda=max(set(sizes), key=sizes.count)
		median=numpy.median(numpy.array(sizes))

		modaStats.append(moda)
		
		
		#check if the gene is conserved considering the threshold and the -p parameter
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



		sizes.append(gene)
		sizes.append(sizes[0])
		sizes.append(median)
		
	
	
		allsizes.append(sizes)
		allNumberAlleles.append([gene,allelenumber])

	with open("tabStats.txt", "wb") as t:
		t.write("gene\tmin\tmax\tmean\tstandard deviation")
		for elem in tabStats:
			for elem2 in elem:
				t.write(str(elem2)+"\t")
			t.write("\n")
	print "\n"+str(z)+ " genes with only one allele\n"
	
	#order genes by median
	sortbymedia=sorted(allsizes, key=itemgetter(-1))
	sortbymedia.reverse()
	
	#order genes by number of alleles
	sortbyNumberAlleles=sorted(allNumberAlleles, key=itemgetter(-1))
	sortbyNumberAlleles.reverse()
	
	
	"""allsizes2=deepcopy(allsizes)
	
	#order genes by the first allele
	sortbyfirstallele=[]
	sortbyfirstallele=sorted(allsizes2, key=itemgetter(-2))
	sortbyfirstallele.reverse()"""
	
	
	orderedlistgene=[]

	for elem in sortbymedia:
		elem.pop(-1)
		elem.pop(-1)
		orderedlistgene.append(elem.pop(-1))
		
	
	orderedlistgene2=[]
	for elem in sortbyNumberAlleles:
		orderedlistgene2.append(elem.pop(0))
			
		
	"""orderedlistgeneFirst=[]

	for elem in sortbyfirstallele:
		elem.pop(-1)
		elem.pop(-1)
		orderedlistgeneFirst.append(elem.pop(-1))"""
			
	print str(len(conservedlengthgenes)) +" conserved genes"
	print str(len(notconservedlengthgenes)) +" not conserved genes"

	
	print "\nBuilding box plot..."
	
	if not ReturnValues	:
		with open ("notconserved.txt","wb") as f:
			for gene in notconservedlengthgenes:
				f.write(str(gene)+"\n")


	

	#create the boxplot and build the html representation	
	if not ReturnValues	:
		aplot=buildPlot(sortbymedia,ReturnValues)
		
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
		
		ax.yaxis.labelpad = 40	
		boxplothtml=mpld3.fig_to_dict(fig)
		
		"""
		#code to build the boxplot ordered by first allele of the gene file
		
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
				
		plt.close('all')
		orderedlistgene2_basename=[]
		for elem in orderedlistgene2:
			orderedlistgene2_basename.append(os.path.basename(elem))
		
		
		fig, ax= plt.subplots(figsize=(10.5,5.0))
		points = ax.plot(sortbyNumberAlleles, 'o')
		
		plt.ylabel('Number of alleles')
		frame1 = plt.gca()			
		frame1.axes.get_xaxis().set_visible(False)
		frame1.axes.get_xaxis().set_ticks([])
		ax.yaxis.labelpad = 40
		
		mpld3.plugins.connect(fig, plugins.PointLabelTooltip(points[0],labels=orderedlistgene2_basename))

		mpld3.plugins.connect(fig, ClickInfo2(points[0], orderedlistgene2))
		
		
		
		numberallelesplothtml=mpld3.fig_to_dict(fig)
		plt.close('all')
	
		
		fig, ax = plt.subplots(figsize=(10.5,5.0))
		bp=plt.hist(modaStats,100,rwidth=0.8)
		plt.ylabel('Number of occurrences')
		plt.xlabel('Allele Size')
		ax.yaxis.labelpad = 40
		start, end = ax.get_xlim()
		ticks=range(int(start),int(end),250)
		plt.xticks(ticks)
		plt.grid(True)
		
		
		histplothtml=mpld3.fig_to_dict(fig)


		return notconservedlengthgenes,len(conservedgenes),z,boxplothtml,histplothtml,numberallelesplothtml

if __name__ == "__main__":
    main()
