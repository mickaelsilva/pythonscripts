Allele Calling

=============
CommonFastaFunctions.py and runProdigal.py scripts are needed to run the main scripts
=============
#alleleCalling_probe_based_main.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles as probes

	% alleleCalling_probe_based_main.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` (optional) parameter to return a phyloviz output file type and a statistics.txt file

short example statistics file:

* EXC - allele has exact match (100% identity)
* NA - new allele found
* undefined - allele found is contained in an already defined allele but match size is more than 2 bases different from the defined allele
* LNF - locus not found
* LOT - locus is on the tip of the contig of the genome
* incomplete - match size is less than 80% of the allele size and identity% is smaller than 50%
* small - match size is less than 50% of the allele size
```
Stats:	EXC	NA	undefined	LNF	LOT	incomplete	small
NC_017162.fna	1026	4	0	0	0	0	0	
NC_009085.fna	248	1	0	0	0	0	0	

```
short example phyloviz file output:
```
FILE	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384132717_ref_YP_005515329.1_.fasta	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133246_ref_YP_005515858.1_.fasta
10_S10_L001.fasta	7	1
2_S12_L001.fasta	7	1
7_S7_L001.fasta	NA6:-8	NA5:-7
```
=============
#alleleCalling_ORFbased_main.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles as probes and using Prodigal corrects the BLAST hit of the alleles against the genome

	% alleleCalling_ORFbased_main.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` (optional) parameter to return a phyloviz output file type and a statistics.txt file

short example statistics file:

* EXC - allele has exact match (100% identity)
* INF - infered allele with prodigal
* LNF - locus not found
* LOT - locus is on the tip of the contig of the genome
* incomplete - match size is less than 80% of the allele size and identity% is smaller than 50%
* small - match size is less than 50% of the allele size
```
Stats:	EXC	INF	LNF	LOT	incomplete	small
NC_017162.fna	1026	193	0	0	0	0	
NC_009085.fna	248	974	0	0	0	0	

```
short example phyloviz file output:
```
FILE	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384132717_ref_YP_005515329.1_.fasta	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133246_ref_YP_005515858.1_.fasta
10_S10_L001.fasta	7	1
2_S12_L001.fasta	7	1
7_S7_L001.fasta	INF1:-8	INF1:-7
```
=============
#alleleCalling_ORFbased_protein_main.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles translated aminoacid sequences and the ORF translated sequences of a genome obtained using Prodigal

	% alleleCalling_ORFbased_protein_main.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
`-i` path to the list of genomes file

`-g` path to the list of alleles file

`-o` output file name

`-p` (optional) parameter to return a phyloviz output file type and a statistics.txt file

short example statistics file:

* EXC - allele has exact match (100% identity)
* INF - infered allele with prodigal
* LNF - locus not found
```
Stats:	EXC	INF	LNF
NC_017162.fna	1026	193	
NC_009085.fna	248	974	

```
short example phyloviz file output:
```
FILE	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384132717_ref_YP_005515329.1_.fasta	Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133246_ref_YP_005515858.1_.fasta
10_S10_L001.fasta	7	1
2_S12_L001.fasta	7	1
7_S7_L001.fasta	INF1:-8	INF1:-7
```
