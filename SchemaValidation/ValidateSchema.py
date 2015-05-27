#!/usr/bin/python
import CheckCDS
import alleleSizeStats
import os
import argparse


def main():
	
	parser = argparse.ArgumentParser(description="This program analyses a set of gene files, analyzing the alleles CDS and the length of the alleles per gene")
	parser.add_argument('-i', nargs='?', type=str, help='list genes, directory or .txt file with the full path', required=True)
	parser.add_argument('-p', nargs='?', type=bool, help='One bad allele still makes gene conserved', required=False)
	parser.add_argument('-t', nargs='?', type=float, help='Threshold', required=False)
	
	args=parser.parse_args()
	genes = args.i
	
	try:
		threshold=float(args.t)
	except:
		threshold=0.05
		pass
	try:
		OneBadGeneNotConserved=bool(args.p)
	except:
		OneBadGeneNotConserved=False
		pass
	
	try:
		f=open( genes, 'r')
		f.close()
	except IOError:
		listbasename=os.path.basename(os.path.normpath(genes))
		
		with open("listGenes"+listbasename+".txt", "wb") as f:
			for file in os.listdir(genes):
				f.write( str(genes)+str(file)+"\n")
		genes="listGenes"+listbasename+".txt"
		
	
	genebasename=str(os.path.basename(genes))
	genebasename=genebasename.split(".")
	genebasename=genebasename[0]
	imagesDir="./resultsHTML/"+genebasename+"/"
	
	
	if not os.path.exists(imagesDir):
		os.makedirs(imagesDir)
		
	notConservedgenes,totalgenes,genesWOneAllele=alleleSizeStats.getStats(genes,threshold,OneBadGeneNotConserved,True)
	print notConservedgenes
	listStopc,listnotStart,listnotMultiple=CheckCDS.analyzeCDS(genes,True)
	
	with open("./resultsHTML/"+genebasename+"_results.html", "wb") as f:
		f.write("<!DOCTYPE html>\n<html>\n<head>\n<title>Schema Validation Results</title>\n</head>\n<body>\n<h1>Allele CDS analysis results</h1>\n")
		f.write("<h2>Alleles with stop codon inside</h2>\n<ul style='list-style-type:circle'>\n")		
		for item in listStopc:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		f.write("<h2>Alleles without start codon</h2>\n<ul style='list-style-type:circle'>\n")		
		for item in listnotStart:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		
		f.write("<h2>Not multiple of 3 alleles</h2>\n<ul style='list-style-type:circle'>\n")		
		for item in listnotMultiple:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		
		f.write("<h1>Allele size analysis using a mode +/- "+str(threshold)+"</h1>\n<h3>"+str(totalgenes)+" total genes</h3>\n<h3>"+str(len(notConservedgenes))+" genes not conserved</h3>\n<h3>"+str(genesWOneAllele)+" genes with only one allele</h3>\n")
		f.write("<h2>Distribution of allele mode sizes per gene</h2>\n<img src='./"+genebasename+"/histogram_mode.png'>\n<br>")
		f.write("<h2>Boxplot for all genes in log10 scale</h2><p>Box plot for each gene on a descending order of the maximum allele size</p><p>-->Box represent the 25 and 75 percentiles (1st and 3rd quartile)</p><p>-->Box plot whiskers representing the 5 and 95 percentile</p><p>-->Red line represent the median (2nd quartile)</p><p>-->Green dots are outliers </p>\n<img src='./"+genebasename+"/plot.png'>\n")
		
		
		
		i=1
		
		while i<100:
			fname=imagesDir+"plot"+str(i)+".png"
			if os.path.isfile(fname):
				i+=1
			else:
				break
		j=1
		f.write(" \n<br><h2>Get a closer view of the results?</h2><button onclick='showImages()'>Show more detailed images</button>\n<div id='detailedImages' style='display: none;'>")
		f.write("\n<h2>batch of 100 genes box plot in linear scale</h2>")
		while j<i:
			fname="./"+genebasename+"/plot"+str(j)+".png"
			f.write("<img src='"+fname+"'>\n")
			j+=1
		f.write("</div>\n<script type='text/javascript'>function showImages(){var t=document.getElementById('detailedImages').removeAttribute('style');}</script></body>\n</html>")
	
if __name__ == "__main__":
    main()

