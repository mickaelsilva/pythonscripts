SeqOfWeb
=============

#seqFromWebTaxon.py

Given a taxon name (recognized by SRA webservice as a taxon), provides a list (-o parameter is the name of the output file) of the identifiers of sequencing runs :

Running the scrip:

	% SeqFromWebTaxon.py -i "acinetobacter baumannii" -o allfastq.txt -g True
	
output file example:
```
DRR006296	Illumina Genome Analyzer IIx
DRR006297	Illumina Genome Analyzer IIx
ERR110008	Illumina HiSeq 2000
ERR110009	Illumina HiSeq 2000
ERR110010	Illumina HiSeq 2000
```

Running the scrip:

	% SeqFromWebTaxon.py -i "acinetobacter baumannii" -o allfastq.txt

```
DRR006296
DRR006297
ERR110008
ERR110009
ERR110010
```
the output file may be used as an input file to the seqFromWeb.py script

Parameters

`-i` Taxon name
`-o` output file
`-g` (optional) True to include the sequencer information on the output

=============

#seqFromWeb.py

Given a list of ERR id's (listOfIDs.txt), downloads all the files represented by the ERR to target dir (dirToDown).
example of a working ERR:

ERR051577

which will be downloaded as ftp://ftp.sra.ebi.ac.uk/vol1/fastq/ERR051/ERR051577/

example listOfIDs.txt:
```
ASDASDASD
ERR047928
ASD
ERR048171
```
Running the scrip:

	% seqFromWeb.py -i listOfIDs.txt -g dirToDown
	
ERR are tested to check their existance and will promp a message such as "BAD ID:" if not found

message promp example:
```
Bad ID: ASDASDASD
Bad ID: ERR047928
Bad ID: ASD
Downloading: ERR048171_1.fastq.gz
```
Using the optional parameter -m

Running the scrip:

	% seqFromWeb.py -i listOfIDs.txt -g dirToDown -m "Illumina"
	
with the input:

```
ERR586979	Illumina HiSeq 2000
ERR586980	Illumina HiSeq 2000
ERR662807	Illumina HiSeq 2000
SRR006006	454 GS FLX
SRR006007	454 GS FLX
SRR006008	454 GS FLX
```

only the first 3 runs will be downloaded

Parameters

`-i` list of run ID's to download
`-g` directory where the files will be downloaded
`-m` (optional) Jumps the runs that are not sequenced by this type of sequencer


