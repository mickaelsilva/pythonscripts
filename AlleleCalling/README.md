Allele Calling

=============
#alleleCalling_baseprobe.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles as probes

	% alleleCalling_baseprobe.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
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
#alleleCalling_ORFbased.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles as probes and using Prodigal corrects the BLAST hit of the alleles against the genome

	% alleleCalling_ORFbased.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
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
#alleleCalling_ORFbased_protein.py

Given a list of genomes and a list of alleles, the program will perform an allele call using the defined alleles translated aminoacid sequences and the ORF translated sequences of a genome obtained using Prodigal

	% alleleCalling_ORFbased_protein.py -i listgenomes.txt -g listgenes.txt -o output_file_name.txt -p True
	
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

=============
#ParalogRemove.py

Given an ffn file, recovers the genes that are not paralogs and have a size bigger than the g parameter provided

	% ParalogRemove.py -i ../ffnallele/NC_010611.ffn -g 10000
	
`-i` .ffn file

`-g` minimum size of the gene

=============
#CompareSameLocus.py

Given two list of locus from two different schemas and creates a folder with paired files for all the locus that are determined to be the same based on BLAST search with identity greater than 80% and locus length variation up to 20%
	
	% CompareSameLocus.py -i ./listgenes6.txt -g listgenesffn.txt

short example of a list genes:
```
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384132717_ref_YP_005515329.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133246_ref_YP_005515858.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_AB0057.1.peg.gi_213158168_ref_YP_002320219.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_AB0057.1.peg.gi_213156045_ref_YP_002318090.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384131122_ref_YP_005513734.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133544_ref_YP_005516156.1_.fasta
/home/msilva/Desktop//Unique_Acinetobacter_baumannii_1656-2.1.peg.gi_384133133_ref_YP_005515745.1_.fasta
```

=============
#CreateSchema.py

dependencies:
CommonFastaFunctions.py
biopython
HTSeq

Given a concatenated ffn file, removes genes that are substring of bigger genes and genes smaller than choosen in the -g parameter. Blasts all the genes against each other and saves the bigger genes, removing the smaller genes with a 0.6>BSR

	% CreateSchema.py -i allffnfile.fasta -g 200
	
Output:
proteins.fasta containing the transaltion of all the genes from the given ffn file, without substring genes
*.fasta large set of .fasta files, 1 per gene

