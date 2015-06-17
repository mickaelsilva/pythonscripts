#!/usr/bin/python
import CheckCDS
import alleleSizeStats
import os
import argparse
import json

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
	
		
	notConservedgenes,totalgenes,genesWOneAllele,boxplot,histplot=alleleSizeStats.getStats(genes,threshold,OneBadGeneNotConserved,True)
	
	boxplot=str(json.dumps(boxplot))
	histplot=str(json.dumps(histplot))

	listStopc,listnotStart,listnotMultiple=CheckCDS.analyzeCDS(genes,True)
	
	with open("./resultsHTML/"+genebasename+"_results.html", "wb") as f:
		f.write("<!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script><script type='text/javascript' src='https://mpld3.github.io/js/mpld3.v0.2.js'></script>\n")
		f.write("""<style type="text/css">
		ul {
    /*min-height: 300px;*/
    -webkit-column-count: 4;
       -moz-column-count: 4;
            column-count: 4; /*4 is just placeholder -- can be anything*/
}
li {
    display: table;
    padding-bottom: 20px; 
    margin-right: 30px;
}
li a {
    color: rgb(0, 162, 232);
}
</style>""")

		
		f.write("<title>Schema Validation Results</title>\n</head>\n<body>\n<h1>Allele CDS analysis results</h1>\n")
		f.write("<h2>Alleles with stop codon inside</h2>\n<ul>\n")		
		for item in listStopc:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		f.write("<h2>Alleles without start codon</h2>\n<ul>\n")		
		for item in listnotStart:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		
		f.write("<h2>Not multiple of 3 alleles</h2>\n<ul>\n")		
		for item in listnotMultiple:
			f.write("<li>"+str(item)+"</li>")
			f.write("\n")
		f.write("</ul>\n")
		
		f.write("<h1>Allele size analysis using a mode +/- "+str(threshold)+"</h1>\n<h3>"+str(totalgenes)+" total genes</h3>\n<h3>"+str(len(notConservedgenes))+" genes not conserved</h3>\n<h3>"+str(genesWOneAllele)+" genes with only one allele</h3>\n")
		f.write("<h2>Distribution of allele mode sizes per gene</h2>\n<div id='fig02'></div>\n<br>")
		f.write("<h2>Boxplot for all genes</h2><p>Box plot for each gene on a descending order of the maximum allele size</p><p>Use the zoom button and hover the mouse over a box/median to see the gene name</p><p>-->Box represent the 25 and 75 percentiles (1st and 3rd quartile)</p><p>-->Box plot whiskers representing the 5 and 95 percentile</p><p>-->Red line represent the median (2nd quartile)</p><p>-->Green dots are outliers </p>\n<div id='fig01'></div>\n")
		

		f.write("\n<script type='text/javascript'>var json01 ="+str(boxplot)+";\nvar json02 ="+str(histplot)+";\nmpld3.draw_figure('fig01', json01);mpld3.draw_figure('fig02', json02);</script>")

		f.write("</body>\n</html>")
	
if __name__ == "__main__":
    main()

