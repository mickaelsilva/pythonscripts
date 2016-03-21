#!/usr/bin/python
import CheckCDS
import alleleSizeStats
import os
import argparse
import json
from operator import itemgetter

def main():
	
	parser = argparse.ArgumentParser(description="This program analyses a set of gene files, analyzing the alleles CDS and the length of the alleles per gene")
	parser.add_argument('-i', nargs='?', type=str, help='list genes, directory or .txt file with the full path', required=True)
	parser.add_argument('-p', nargs='?', type=bool, help='One bad allele still makes gene conserved', required=False)
	parser.add_argument('-ta', nargs='?', type=int, help='Threshold', required=True)
	parser.add_argument('-t', nargs='?', type=float, help='Threshold', required=False)
	
	args=parser.parse_args()
	genes = args.i
	transTable = args.ta
	
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

	statsPerGene=CheckCDS.analyzeCDS(genes,transTable,True)
	
	# stats values are ordered in a list allelesNotMultiple3,listStopcodonsInside,listnotStartCodon,numberOfAlleles
	
	
	with open("./resultsHTML/"+genebasename+"_results.html", "wb") as f:
		f.write("<!DOCTYPE html>\n<html>\n<head><script type='text/javascript' src='https://cdnjs.cloudflare.com/ajax/libs/d3/3.5.5/d3.min.js'></script><script type='text/javascript' src='https://mpld3.github.io/js/mpld3.v0.2.js'></script>\n")
		f.write("<script type='text/javascript' src='https://ajax.googleapis.com/ajax/libs/jquery/2.1.3/jquery.min.js'></script>")
		f.write("""<script type='text/javascript'>
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
    </script>""")
		
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

.tg  {border-collapse:collapse;border-spacing:0;border-color:#999;}
.tg td{font-family:Arial, sans-serif;font-size:14px;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#444;background-color:#F7FDFA;}
.tg th{font-family:Arial, sans-serif;font-size:14px;font-weight:normal;padding:10px 5px;border-style:solid;border-width:1px;overflow:hidden;word-break:normal;border-color:#999;color:#fff;background-color:#26ADE4;}
.tg .tg-vn4c{background-color:#D2E4FC;text-align:center}
</style>""")
		
		
		f.write("<h1>Allele size analysis using a mode +/- "+str(threshold)+"</h1><p> Genes are considered not conserved if >1 allele are outside the mode +/-0.05 size. Genes with only 1 allele outside the threshold are considered conserved\n<h3>"+str(totalgenes)+" total genes</h3>\n<h3>"+str(len(notConservedgenes))+" genes not conserved</h3>\n<h3>"+str(genesWOneAllele)+" genes with only one allele</h3>\n")
		f.write("<h2>Distribution of allele mode sizes per gene</h2>\n<div id='fig02'></div>\n<br>")
		f.write("<h2>Boxplot for all genes</h2><p>Box plot for each gene on a descending order of the median allele sizes</p><p>Use the zoom button and hover the mouse over a box/median to see the gene name</p><p>-->Box represent the 25 and 75 percentiles (1st and 3rd quartile)</p><p>-->Box plot whiskers representing the 5 and 95 percentile</p><p>-->Red line represent the median (2nd quartile)</p><p>-->Green dots are outliers </p>\n<div id='fig01'></div>\n")
		
		f.write("<title>Schema Validation Results</title>\n</head>\n<body>\n<h1>Allele CDS analysis results</h1>\n<p>Summary table of the alleles with issues per gene using the <a href='http://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi#SG11'>NCBI translation table 11</a></p><p>Click on the gene name to open the fasta file</p><p>click on the boxes with the % to get the index of the alleles with issues</p> ")
		
		
		f.write("""<table class="tg">
  <tr>
    <th class="tg-031e">Gene</th>
    <th class="tg-031e">Number alleles not multiple of 3</th>
    <th class="tg-031e">Number alleles w/ >1 stop codons</th>
    <th class="tg-031e">Number alleles wo/ Start Codon</th>
    <th class="tg-031e">Number of alleles (% alleles w/ issues) </th>
  </tr>""")
		ordered=[]
		for key, value in statsPerGene.iteritems():
			aux=[]
			numberMultip=float(len(value[0]))
			numberStop=float(len(value[1]))
			numberStart=float(len(value[2]))
			total=float(value[3])
			totalpercent=((numberMultip+numberStart+numberStop)/total)*100
			aux.append(key)
			aux.append(totalpercent)
			ordered.append(aux)
		ordered=sorted(ordered, key=itemgetter(-1))
		ordered.reverse()
		
		newlist=[]	
		for item in ordered:	
			
			aux=[]
			aux.append(item[0])
			
			
			
			value=statsPerGene[item[0]]
			numberMultip=float(len(value[0]))
			numberStop=float(len(value[1]))
			numberStart=float(len(value[2]))
			total=float(value[3])
			
			aux.append(value)
			newlist.append(aux)
			name=os.path.basename(str(item[0]))
			name=name.split(".")
			name=name[0]
			if (numberMultip>0 or numberStop>0 or numberStart>0):
				f.write("<tr id="+str(item[0])+""">\n<td class='tg-vn4c' onclick="window.open('"""+str(item[0])+"""')" style='cursor:pointer'>"""+name+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberMultip))+" ("+str('{0:.2f}'.format((numberMultip/total)*100))+"%)"+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberStop))+" ("+str('{0:.2f}'.format((numberStop/total)*100))+"%)"+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(numberStart))+" ("+str('{0:.2f}'.format((numberStart/total)*100))+"%)"+"</td>\n<td class='tg-vn4c' onclick='a(this);'>"+str(int(total))+" ("+str('{0:.2f}'.format(((numberMultip+numberStart+numberStop)/total)*100))+"%)"+"</td>\n</tr>")
		
		
		f.write("</table>")
		f.write("""<div id='AllelesWissues'></div><button onclick="$('#AllelesWissues').empty();">clean</button>""")
		f.write("""\n<script type='text/javascript'>function a(element) {
	var id = $(element).closest("tr").attr("id");
	var badalleles=[];
	for (i = 0; i < alleles.length; i++) { 
		if((alleles[i])[0]==id){
			badalleles=(alleles[i])[1];
			break;
			}
		}
	var notmulti=(badalleles[0]).join('; ');
	var stopcodon=(badalleles[1]).join('; ');
	var startcodon=(badalleles[2]).join('; ');
	var name=(id.split("/")).slice(-1)[0]
	name=(name.split("."))[0]
	$('#AllelesWissues').append('<h2> Gene: '+name+'</h2>');
	$('#AllelesWissues').append('<p> Alleles not multiple of 3: '+notmulti+'</p><p> Alleles with >1 stop codon: '+stopcodon+'</p><p>Alleles without start codon: '+startcodon+'</p>');

	$('html,body').animate({
        scrollTop: $('#AllelesWissues').offset().top},'slow');
	
	}
	</script>""")
		
		f.write("\n<script type='text/javascript'>var alleles="+json.dumps(newlist)+"</script>")
		f.write("\n<script type='text/javascript'>var json01 ="+str(boxplot)+";\nvar json02 ="+str(histplot)+";\nmpld3.draw_figure('fig01', json01);mpld3.draw_figure('fig02', json02);</script>")

		f.write("</body>\n</html>")
	
if __name__ == "__main__":
    main()

