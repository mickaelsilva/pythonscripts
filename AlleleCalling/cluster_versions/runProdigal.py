#!/usr/bin/python
import sys
import os
import os.path
import string
import subprocess
import pickle

def main():
	
	try:
		input_file = sys.argv[1]
		tempPath = sys.argv[2]
	except IndexError:
		print "usage: list_pickle_obj"
	
	contigsFasta=input_file
	
	basepath=tempPath
	print basepath
	
        # ------------ #
        # RUN PRODIGAL #
        # ------------ #  #/scratch/NGStools/prodigal-2.50/prodigal -i /scratch/spneumoniae_600/ERR067978/velvet_assembly/contigs.fa -c -m -g 11 -p single -f sco -q > test_all.txt
	prodigal_path='prodigal'
	proc = subprocess.Popen([prodigal_path, '-i', contigsFasta, '-c', '-m', '-g', '11', '-p', 'single', '-f', 'sco', '-q'], stdout=subprocess.PIPE)
	
	cdsDict = {}
	tempList = []
	line = ' '
	while line != '':
		
		# when it finds a contig tag
		if "seqhdr" in line:
		# add contig to cdsDict and start new entry
			
			if len(tempList) > 0:

				 # --- brute force parsing of the contig tag - better solution is advisable --- #
				
				i=0
				for l in contigTag:
					if l == ' ':
						break
					i+=1
				contigTag=contigTag[:i]
				
				cdsDict[contigTag.replace("\r","")] = tempList
				tempList = []
			
			contigTag = line.split('"')[-2]
			
                # when it finds a line with cds indexes
		elif line[0] == '>':
			
                        # parsing
			cdsL = line.split('_')

                        # --- each element of this list is a pair of indices - the start and the end of a CDS --- #

			tempList.append([ int(cdsL[1]) - 1 , int(cdsL[2]) ])		# start index correction needed because prodigal indexes start in 1 instead of 0

                # reads the stdout from 'prodigal'
		line = proc.stdout.readline()

	# ADD LAST
	if len(tempList) > 0:
		
		 # --- brute force parsing of the contig tag - better solution is advisable --- #

		i=0
		for l in contigTag:
			if l == ' ':
				break
			i+=1
		contigTag=contigTag[:i]

		cdsDict[contigTag.replace("\r","")] =tempList
	#print cdsDict.keys()
	
	filepath=os.path.join(basepath,str(os.path.basename(contigsFasta))+"_ORF.txt")
	with open(filepath, 'wb') as f:
		var = cdsDict
		pickle.dump(var, f)
	
	
	return True
	
if __name__ == "__main__":
    main()
