#! /usr/bin/env python

import sys
import subprocess
import ftplib
import argparse
import os.path

def download(dirs2,target_dir2,ref2,success2,f2,link2):
	#new folder for each reference with reference id name
	subprocess.call(['mkdir', target_dir2+"/"+ref2])
				
	#get fasta file for each read file name
	for item in dirs2:
					
		f2.cwd(link2)
		final_target_dir=target_dir2+"/"+ref2 +"/"+ item
		file = open(final_target_dir, 'wb')
		print "Downloading: %s" % item

		f2.retrbinary('RETR %s' % item, file.write)
		file.close()
		print "Downloaded %s" % item
		success2+=1		

	return success2

def main():
	
	#run -i seqFromWeb.py -i allfastq.txt -g ./downloadseq -m "Illumina"
	parser = argparse.ArgumentParser(description="This program downloads sequencing runs given the sra RUN ID in a list to a selected directory")
	parser.add_argument('-i', nargs='?', type=str, help='RUN ID list', required=True)
	parser.add_argument('-g', nargs='?', type=str, help='target directory', required=True)
	parser.add_argument('-m', nargs='?', type=str, help='sequencing machine name', required=False)
	
	args=parser.parse_args()

	
	target_id_file = args.i
	target_dir = args.g
	
	if not os.path.exists(target_dir):
		os.makedirs(target_dir)
	
	try:
		machinetype = args.m
		machinetype=machinetype.lower()
	except :
		pass	
	
	with open(target_id_file) as f:
		references = f.read().splitlines()
		
	failed=0
	success=0
	f=ftplib.FTP('ftp.sra.ebi.ac.uk')
	f.login()
	#for each reference id get the read files
	for ref in references:
		toDown=True
		try:
			ref=ref.split("\t")
			model=ref[1]
			ref=ref[0]
			
		except:
			
			continue
		else:
			if machinetype and machinetype not in model.lower():

				toDown=False
		
		if toDown:
			
			try:
				
				firstid=ref[0:6]
				#get the read files name from the reference id
				link='/vol1/fastq/'+firstid+"/"+ref
				f.cwd(link)
				dirs=f.nlst();
				
			except:
				try:
					firstid=ref[0:6]
					#get the read files name from the reference id
					link='/vol1/fastq/'+firstid+"/00"+ref[-1]+"/"+ref
					f.cwd(link)
					dirs=f.nlst();
					
					
						
				except Exception, e:
					failed +=1
					print "Bad ID: " + ref
				else:
					success=download(dirs,target_dir,ref,success,f,link)	
			
			else:
				
				success=download(dirs,target_dir,ref,success,f,link)
						
				
	f.quit()	
	print "Successfully downloaded %s files and %s ID references were wrong" % (success,failed)	
			
if __name__ == "__main__":
    main()

