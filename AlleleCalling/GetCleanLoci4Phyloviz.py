import sys
import csv
import numpy as np
from numpy import array
import argparse
from collections import OrderedDict
import HTSeq


def presence (d2):	

	d2c=np.copy(d2)
	
	#genelist= d2[:1,1:]
	for item in d2c:
		print len(item)
	genomeslist= d2c[1:,:1]
	
	d2c = d2c[1:,1:-1]
	
	row=0
	while row<d2c.shape[0]:
		column=0
		while column<d2c.shape[1]:
			if d2c[row,column] == "LNF":
				d2c[row,column]=0
			else:
				d2c[row,column]=1
			column+=1
		row+=1
	
	#print d2
	print len(d2c)
	d2c=d2c.T
	print len(d2c)
	genomeslist=(genomeslist.tolist())

	return d2c,genomeslist
	
def lostGenesPontuation(matrix,genomeslist):
	
	numbergenomes=len(genomeslist)
	pontuationmatrix=[0]*numbergenomes
	
	for gene in matrix:
		aux=[0]*numbergenomes
		presentgenes=0
		genomeindex=0
		for value in gene:
			if int(value)==0:
				aux[genomeindex]=aux[genomeindex]-1
				#print aux
				
			else:
				presentgenes+=1
			genomeindex+=1
		#if presentgenes>(numbergenomes-(numbergenomes/2)):
		if presentgenes>(numbergenomes-3):		
			pontuationmatrix=[x + y for x, y in zip(pontuationmatrix, aux)]
				
			
	
	cleangenomelist=[]
	for genome in genomeslist:
		cleangenomelist.append(genome[0])
		
	ordered=OrderedDict(zip(cleangenomelist,pontuationmatrix ))
	ordered = OrderedDict(sorted(ordered.items(), key=lambda(k,v):(v,k)))
	ordered=OrderedDict(ordered.items()[::-1])
	#print ordered

	for genome in ordered.items():
		print genome
	return False

def notGoodAlleles(genesDirect,geneName, rangeFloat):
	
	sizes=[]
	notgood=[]
	
	filename=str(genesDirect)+str(geneName)
	print filename
	
	gene_fp2 = HTSeq.FastaReader(filename)
	allelenumber=0
	sizesum=0
	for allele in gene_fp2: 
		#print allele
		allelesize=len(allele.seq)
		sizes.append(allelesize)
		sizesum+=allelesize
		allelenumber+=1
	mean=float(sizesum)/allelenumber
		
		
	mediana=np.percentile(sizes, 50, axis=0)	
	#quart3=np.percentile(sizes, 75, axis=0)	
	
	i=0
	#print sizes
	#print quart1,quart3
	
	#print mean
		
	for size in sizes:
			
		if ( float(size)> mediana*(1+rangeFloat) or float(size)< mediana*(1-rangeFloat)):
			notgood.append(i+1)
		i+=1
	
	#print notgood		
	return notgood
	
def removegenomes(d2a,shortGenomeList):
	
	rowid=1
	deleted=0
	
	while rowid< d2a.shape[0]:
		inside=False
		for genome in shortGenomeList:
			#print genome
			
			if genome == d2a[rowid][0] :	
				inside=True
				break

			#columnid+=1
		if not inside:
			d2a=np.delete(d2a, rowid, axis=0)
			deleted+=1
		rowid+=1
	#print d2a
	#print len (d2a)
	
	if deleted>0:
		return removegenomes(d2a,shortGenomeList)
	else:
		#print d2a
		return d2a

def clean (inputfile,outputfile,totaldeletedgenes,shortGenomeList,genesDirect,rangeFloat):
	
	with open(inputfile) as f:
		reader = csv.reader(f, delimiter="\t")
		d = list(reader)
	
	d2 = array(d)
	if shortGenomeList:
		d2=removegenomes(d2,shortGenomeList)
	
	matrix2,genomeslist2=presence (d2)

	lostGenesPontuation(matrix2,genomeslist2)
	#print d2
	d2=d2.T
	rowid=1
	deleted=0
	badAlleles=False
	lnfdel=0
	balldel=0
	while rowid< d2.shape[0]:
		columnid=0
		
		if genesDirect and columnid==0:
			badAlleles=notGoodAlleles(genesDirect,d2[rowid][columnid], rangeFloat)
			columnid=+1
		#print badAlleles
		while columnid< d2.shape[1]:
			
			if ( "undefined" in d2[rowid][columnid] or "incomplete" in d2[rowid][columnid] or "LNF" in d2[rowid][columnid] or "small" in d2[rowid][columnid] or "LOT" in d2[rowid][columnid] ):	
				d2=np.delete(d2, rowid, 0)
				totaldeletedgenes+=1
				deleted+=1
				rowid-=1
				columnid=1
				lnfdel+=1
				break
			
			if badAlleles :
				allele=d2[rowid][columnid]
				allele=allele.replace("INF-","",)
				
				try:
					if int(allele) in badAlleles:
						#print int(allele),badAlleles
						d2=np.delete(d2, rowid, 0)
						totaldeletedgenes+=1
						deleted+=1
						rowid-=1
						columnid=1
						balldel+=1
						break
				except:
					pass
			
			columnid+=1
		rowid+=1
	
	d2=d2.T
	d2=d2.tolist()
	
	
	with open(outputfile, "wb") as f:
		writer = csv.writer(f,delimiter='	')
		writer.writerows(d2)
	
	

	file = open(outputfile)
	contents = file.read()
	contents = contents.replace('INF-', '')
	contents = contents.replace('INF1:-', '')
	contents = contents.replace('INF2:-', '')
	contents = contents.replace('INF3:-', '')
	contents = contents.replace('INF4:-', '')
	contents = contents.replace('INF5:-', '')
	contents = contents.replace('INF6:-', '')
	contents = contents.replace('INF7:-', '')
	contents = contents.replace('NA-', '')
	contents = contents.replace('NA1-', '')
	contents = contents.replace('NA2-', '')
	contents = contents.replace('NA3-', '')
	contents = contents.replace('NA4-', '')
	contents = contents.replace('NA5-', '')
	contents = contents.replace('NA6-', '')
	
	with open(outputfile, 'w') as f:
			f.write(contents)
	
	#if deleted==0:
	print str(lnfdel),str(balldel)
	print "deleted : %s genes" % totaldeletedgenes
	print "total genes remaining : "+ str(rowid)
	#else:
		#clean(outputfile,outputfile,totaldeletedgenes,shortGenomeList,genesDirect,rangeFloat)
	

def main():

	parser = argparse.ArgumentParser(description="This program cleans an output file for phyloviz")
	parser.add_argument('-i', nargs='?', type=str, help='output to clean', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='name of the clean file', required=True)
	parser.add_argument('-o', nargs='?', type=str, help='info file', required=False)
	parser.add_argument('-p', nargs='?', type=int, help='group by property', required=False)
	parser.add_argument('-s', nargs='?', type=str, help='specify group', required=False)
	parser.add_argument('-d', nargs='?', type=str, help='genes directory', required=False)
	
	args = parser.parse_args()

	
	pathOutputfile = args.i
	newfile = args.g
	
	
	groupselected=False
	try:
		info = args.o
		property = (args.p)-1
		groupselected = args.s
		groupselected=groupselected.split(',')
	except:
		pass
	genesdirect=False
	try:
		genesdirect = (args.d)
	except:
		pass
	print groupselected
	
	if groupselected:
		with open(info) as f:
			reader = csv.reader(f, delimiter="\t")
			oldlist = list(reader)
			
		print "sorted by : " +str((oldlist[0])[property])
		oldlist.remove(oldlist[0])

		
		values = set(map(lambda x:x[property], oldlist))
		#print values
		newlist = [[y[0] for y in oldlist if y[property]==x] for x in values]
	
		shortGenomeList=[]
		values=list(values)
		for value in groupselected:
			index=values.index(value)
			listvalue=newlist[index]
			for genome in listvalue:
				shortGenomeList.append(genome)
	else:
		shortGenomeList=False
		
	print shortGenomeList
	
	
	
	#print matrix
	
	clean(pathOutputfile,newfile,0,shortGenomeList,genesdirect,0.2)

	

	
	
	
	
		
if __name__ == "__main__":
	main()
