README.txt

-----
Purpose of this repository is to provide the tools to reproduce the analysis performed for the paper "Genome-wide mutant fitness profiling predicts the mechanism of action of a Lipid II binding antibiotic" 

CITATION: []
-----

FASTQ TO TABULAR FILE: HOW TO USE GALAXY TO TRIM AND FILTER SEQUENCES AND ALIGN TO GENOME
FastQsanger files can be found at NCBI SRA [reference code]
Galaxy is a an open source, web-based platform for data intensive biomedical research. The following analyses were performed on the Galaxy Service hosted by Tufts Genomics. 

1) GET DATA: Get Data > Upload File
	- File Format: fastqsanger
	- Choose fastq file of interest
2) SPLIT BY TRANSPOSON CONSTRUCT. Our transposon libraries have six different transposon constructs. SEE: Santiago, M; Matano, LM; Moussa SH; Gilmore, MS; Walker, S; Meredith, TC. A new platform for ultra-high density
Staphylococcus aureus transposon libraries. BMC Genomics (2015): 
	- UPLOAD BARCODES: Get Data > Upload File
		- File Format: txt
		- Choose Galaxy-TransposonBarcodes.txt
	- SPLIT BY BARCODE: NGS: QC and manipulation > Barcode Splitter
		- Barcodes to use: Galaxy-TransposonBarcodes.txt
		- Library to split: fastqsanger file
		- Barcdes found at End of Sequence (3' end)
		- Number of allowed mismatches: 2
		- Number of allowed barcodes nucleotide deletions: 0
	- UPLOAD UNMATCHED: Click eye icon on Barcode Splitter output. Will display a series of hyperlinks to the split data. Copy the link for the one labeled 'unmatched', go to Get Data > Upload File, and upload via the link as a fastqsanger file
	- TRIM UNMATCHED: the MmeI enzyme used in our sample prep cuts either 15 or 16 bases from its recognition site. The first barcode split recognizes sequences that are 16 bases long. We must now get the sequences that were 15 bases long. 
		- NGS: QC and manipulation > Trim sequences
			- Library to clip: uploaded unmatched fastqsanger
			- First base to keep: 1
			- Last base to keep: 75
	- REPEAT BARCODE SPLIT using the output from the Trim sequences. 
	- UPLOAD MATCHING TRANSPOSON FILES: construct 1 from first barcode split and from second barcode split. This can be done by clicking on the eye icon for the barcode splitting, copying the appropriate data hyperlink, and uploading via Get Data > Upload File. 
	- CONCATENATE MATCHING TRANSPOSON FILES: Text Manipulation > Concatenate datasets tail-to-head. 
		- Concatenate Dataset: Uploaded file from first barcode split 
		- Click Insert Dataset
		- 1: Dataset: Uploaded file from second barcode split
3) TRIM DOWN TO JUST GENOMIC SEQUENCE: 
	- NGS: QC and manipulation > Trim sequences
		- Library to clip: concatenated dataset 
		- First base to keep: 5
		- Last base to keep: 20
4) FILTER BY QUALITY: 
	- NGS: QC and manipulation > Filter by quality
		- Library to filter: Trimmed sequence
		- Quality cut-off value: 20 
		- Percent of bases in sequence that must have quality equal to / higher than cut-off value: 90
5) ALIGN READS TO GENOME: 
	- Upload NCTC8325.fasta using Get Data > Upload File
	- NGS: Mapping > Map with Bowtie for Illumina
		- Will you select a reference genome from your history or use a built-in index?: Use one from the history
		- Select the reference genome: NCTC8325.fasta
		- Choose whether to use Default options for building indices or to Set your own: Default
		- Is this library mate-paired?: Single-end
		- FASTQ file: Filtering output file 
		- Bowtie settings to use: Full parameter list
		- Skip the first n reads (-s): 0
		- Only align the first n reads (-u): -1
		- Trim n bases from high-quality (left) end of each read before alignment (-5): 0
		- Trim n bases from low-quality  (right) end of each read before alignment (-3): 0 
		- Maximum number of mismatches permitted in the seed (-n): 2
		- Maximum permitted total of quality values at mismatched read positions: 70
		- Seed length (-l) 8
		- Whether or not to round to the nearest 10 and saturating at 30 (--nomaqround): Round to nearest 10 
		- Number of mismatches for SOAP-like alignment policy (-v): -1
		- Whether or not to try as hard as possible to find valid alignments when they exist (-y): Try hard
		- Report up to n valid alignments per read (-k): 1
		- Whether or not to report all valid alignments per read: Do not report all valid alignments
		- Suppress all alignments for a read if more than n reportable alignments exist (-m): -1
		- write all reads with a number of valid alignments exceeding the limit set with the -m option to a file (--max): No
		- Write all reads that could not be aligned to a file (--un): No
		- Whether or not to make Bowtie guarantee that reported singleton alignments are 'best' in terms of stratum and in terms of the quality of values at the mismatched positions (--best): Do not use best
		- maximum number of backtracks permitted when aligning a read (--maxbts): 125
		- Override the offrate of the index to n (-o): -1
		- Seed for pseudo-random number generator (--seed): -1
		- Suppress the header in the output SAM file: No 
		- Job resource Parameters: Use default job resource parameters
6) TABULATE ALIGNED READS: 
	- If available: Tufts: TnSeq > Hopcount
		- Upload NCTC8325.genbank file using Get Data > Upload File
		- Tufts: TnSeq > Hopcount
			- SAM file (alignment): Output from Bowtie
			- Input is from fastx_collapsed reads (count is in header): No
			- Genbank file containing annotations (genome): Select uploaded Genbank
			- Cutoff for number of reads per insertion to consider: 0
			- Exclude first X positions of each gene: 0
			- Exclude last X positions of gene: 0
			- Exclude first X% of each gene: 0
			- Exclude last X% of each gene: 0
			- Report genome sequence for each position: No
			- Report every position in the genome, even those with zero counts (not recommended): No
	- If Hopcount tool not available: 
		- Download SAM file(s) into directory containing sam_to_tabular.py. 
		- In Command Prompt/Terminal, navigate to directory with SAM file(s). 
		- >> python sam_to_tabular.py

-----

PERFORM PROMOTER ANALYSIS TO IDENTIFY GENOMIC WINDOWS WITH DIRECTIONALLY-BIASED INSERTIONS

1) CREATE DIRECTORY WITH: 
	- Tabular files from Galaxy analysis. Files must be named <treatment>_<transposon construct identifier>.tabular. Recognized Tn construct identifiers: blunt, cap, dual, erm, pen, tuf. TREATMENT FOR CONTROL FILES MUST BE untr
	- HG003_GeneList.txt
	- HG003_TAsites.txt
	- promoter_analysis_driver.py
2) RUN: 
	- In command prompt/terminal, navigate to directory with tabular files. 
	- >> python promoter_analysis_driver.py HG003_GeneList.txt HG003_TAsites.txt

3) OUTPUTS: 
	1. IGV files for each tabular file and a cpt tabular file that has the reads from all of the transposons with a single outward-facing promoter. IGV is a tab-delimited text file with no header and the following fields: 
		0. NCBI genome reference number
		1. Start of TA site
		2. End of TA site
		3. # Total Reads
		4. Locus Tag (if in gene)
		5. # Plus direction reads 
		6. # Minus direction reads 

	2. A promoter analysis text file (starts with PRO_) for each of the transposon files. This is a tab-delimited, headerless text file with the following fields: 
		0. Genomic window number 
		1. TA start
		2. Direction of hit 
		3. Log of the ratio of plus-direction control and experimental reads
		4. Log of the ratio of minus-direction control and experimental reads 
		5. Log of the ratio of plus direction and minus direction experimental reads


-----

PERFORM KNN ANALYSIS TO PREDICT MEChANISM OF ACTION

1) Convert tabular files from Galaxy to IGV format
	- Make directory with tabular files, HG003_TAsites.txt, HG003_GeneList.txt, and tabular_to_igv.py
	- In command prompt/terminal, navigate to directory with tabular files. 
	- >> python tabular_to_igv.py
	- OUTPUT: for each tabular file, you will get a tab-delimited, headerless IGV file with the following columns: 
		0. Genome reference ID
		1. Start of TA site
		2. End of TA site
		3. Total number of reads at that site
		4. Locus Tag if the TA site is in a gene
		5. Number of reads in the plus direction
		6. Number of reads in the minus direction

2) Use IGV files to create clust_*.txt files that are the input for the KNN algorithm. 
	- Make a directory with:
		- igv_to_clust.py 
		- normalization.py
		- mannwhitneyu.py 
		- fdr.py
		- make_clust_file.py
		- allctrls.igv - a control IGV file to normalize to. 
		- As many experimental igvs as desired, labeled as <treatment>.igv, ideally containing the summed read count from all 6 Tn constructs
		- saouhsc_orfs.txt
	- In command prompt/terminal, navigate to directory 
	- >> python igv_to_clust.py
	- OUTPUTS: 
		- A normalized control IGV for each experimental igv. 
		- A Mann-Whitney U output file for each experimental igv. Tab-delimited with no heading and the following columns: 
			0: Locus tag 
	
			1: Number of TA sites in the gene 
	
			2: Control sequencing reads in the gene 
	
			3: Experiment sequencing reads in the gene
	
			4: Mann Whitney U U value
	
			5: Mann Whitney U p-value 
	
			6: Experiment:Control read ratio 
	
			7: Length of the gene 
	
			8: Control reads per TA site 
	
			9: Experiment reads per TA site 
	
			10: Gene start 
	
			11: Gene stop 
		- A FDR (Benjamini-Hochberg) q-value list for each experimental igv
		- a 'clust' file for each experimental igv file. No heading, comma-delmited, and the following columns: 
			0: Locus Tag
	
			1: Number of TA sites in gene 
	
			2: Length of gene 
	
			3: Control sequencing reads 
	
			4: Experiment sequencing reads 
	
			5: Experiment:Control reads ratio 
	
			6: FDR q-value 
	
			7: Modified read ratio - if either the control reads or experiment reads were less than 1*10-5 times the total read count for the file, the read count was increased to this minimum value before the ratio was calculate. 
	
			8: Index - reads per length.
	
			9: Fitness value - 	gene rank in range (0,1] with interval 1/# of genes based on the product of the index and the modified read ratio.
			
3) Create a list of possible files to use in the KNN training set. This should be a tab-delimited txt file with no heading. Column 1 is a list of file names for the data files to be used in the training set. Column 2 is the numeric category the data file belongs in. Column 3 is empty. Column 4 is a list of categories. Column 5 is the numeric code (0 to [number of categories-1]) for the category.

3b) OPTIONAL: Create a list of files for unknown compounds you want to predict the mechanism of action for. This should be a text file with no heading and one file name per line. 

4) Perform leave-one-out validation and predict mechanism of action of unknowns
	- Create a directory with: 
		- The list created in step 3 (and 3b)
		- The clust files listed in the lists
		- KNN.py
	- In command prompt/terminal, navigate to directory. 
	- >> python KNN.py <training set list> <optional unknown list>
	- OUTPUTS: 
		- predictions.txt: a file listing the predictions for the leave-one-out validation for the training set and their respective probabilities. 
		- filesused.txt: a list of the files that were used for the analysis. Note that some input files are not used if they do not have enough genes that are sufficiently changed compared to the control or if too many of the genes are depleted of sequencing reads compared to the control. 
		- genematrix.txt: a binary matrix showing which files in the training set contributed which genes to the list of gene fitness values fed into the KNN algorithm. 
		- OPTIONAL: PRED_<name of unknowns list>.txt: a list of the predictions for the files with unknown mechanisms of action and their respective probabilities. 

-----

DESCRIPTIONS OF FILES

allctrls.igv: a sum of the reads in all of our control files. Tab-delimited. No header. Columns: 
	0. Genome identifier
	1. Start of TA site
	2. End of TA site
	3. Total number of reads at that site
	4. Locus tag if the TA site is within a gene
	5. Number of reads in the plus strand
	6. Number of reads in the minus strand

fdr.py: A python script that perform Benjamini-Hochberg FDR analysis on the output from the mannwhitneyu.py script.

	REQUIRED INPUTS: mannwhitneyu.py output file - A tab-delimited text file with no heading and the following columns: 
		0: Locus tag 
	
		1: Number of TA sites in the gene 
	
		2: Control sequencing reads in the gene 
	
		3: Experiment sequencing reads in the gene
	
		4: Mann Whitney U U value
	
		5: Mann Whitney U p-value 
	
		6: Experiment:Control read ratio 
	
		7: Length of the gene 
	
		8: Control reads per TA site 
	
		9: Experiment reads per TA site 
	
		10: Gene start 
	
		11: Gene stop  

	OUTPUT: (printed) tab-delimited information with no header and columns 0) locus tag and 1) FDR q-value

	EX: >> python fdr.py mwu_sample.txt > fdr_sample_output.txt

Galaxy-TransposonBarcodes.txt: A tab-delimited text file with '#' comment key. First line is a commented header. Next six lines contain a transposon construct identifier in the first column and the barcode for that transposon construct in the second column.

HG003_GeneList.txt: Tab-delimited text file with no header. Columns: 
	0. Genome reference number
	1. Locus tag
	2. Gene start
	3. Gene stop
	4. Gene name
	5. Gene product description

HG003_TAsites.txt: Tab-delimited text file with several lines of comments. Comment character ';'. Obtained from tool published in  van Helden et al. (2000). Yeast 16(2), 177-187. Data has the following columns: 
	0. Pattern ID 
	1. DNA strand (D for direct)	
	2. Pattern we searched the genome for (TA)
	3. Reference ID for genome we searched 
	4. Start of pattern instance
	5. End of pattern instance
	6. Matching pattern (TA)
	7. Match score

igv_to_clust.py: Python script that normalizes experimental IGV file read counts to a control file, determine the fitness values of genes, and output the 'clust_*.csv' files that are used in the KNN algorithm. 
	DEPENDENCIES: Need to have in same folder the following python scripts: 
		- normalization.py
		- mannwhitneyu.py 
		- fdr.py
		- make_clust_file.py
		- allctrls.igv - a control IGV file to normalize to. Columns are: 
			0. Reference ID for contig 
			1. TA start site 
			2. TA stop site 
			3. Total reads 
			4. Gene, if TA site is in a gene.
			5. Plus reads
			6. Minus reads

		- As many experimental igvs as desired, labeled as treatment.igv
		- saouhsc_orfs.txt - a tab-delimited text file with the columns: 
			0. Locus tag 
			1. Gene start site 
			2. Gene stop site. 
	OUTPUT: For each experimental igv file, outputs a clust_*.txt file. For more information on clust files, see the make_clust_file.py script. 
	EX: >> python igv_to_clust.py

KNN.py: python script that creates a list of informative genes to base KNN categorization on and then performs leave-one-out validation of KNN supervised machine-learning algorithm for predicting antibiotic mechanisms of action. Optionally can also predict the mechanism of action of a series of unknown antibiotics. 
	REQUIRED INPUT: A tab-delimited txt file. Column 1 is a list of file names for the data files to be used in the training set. Column 2 is the numeric category the data file belongs in. Column 3 is empty. Column 4 is a list of categories. Column 5 is the numeric code (0 to [number of categories-1]) for the category.
		The data files are comma-delimited files with the following columns and must be in the current directory: 
			0. Gene locus tag 
			1. Number of TA sites in the gene
			2. Gene length
			3. Number of control sample reads 
			4. Number of treatment sample reads 
			5. Raw treatment:control read ratio 
			6. FDR-corrected p-value
			7. Corrected ratio. For each the control and treatment files, the minimum read count for each gene is set to 1/10000 x the total read count of the file. The treatment:control read ratio is then taken. 
			8. Treatment sample reads per length of gene 
			9. Normalized rank - col 7*8
	OPTIONAL INPUT: A txt file listing the data files with unknown mechanism of action that we want to predict a category for. 
	OUTPUTS: 
		- predictions.txt: a file listing the predictions for the leave-one-out validation for the training set and their respective probabilities. 
		- filesused.txt: a list of the files that were used for the analysis. Note that some input files are not used if they do not have enough genes that are sufficiently changed compared to the control or if too many of the genes are depleted of sequencing reads compared to the control. 
		- genematrix.txt: a binary matrix showing which files in the training set contributed which genes to the list of gene fitness values fed into the KNN algorithm. 
		- OPTIONAL: PRED_<name of unknowns list>.txt: a list of the predictions for the files with unknown mechanisms of action and their respective probabilities. 
	EX: >>python KNN.py KNN_data_files.txt unknown_data_files.txt

make_clust_file.py: python script that creates the input for our KNN algorithm, known as a clust file, converting read ratios into fitness values. 

	REQUIRED INPUTS:  

		- Mann-Whitney U file: output from mannwhitneyu.py script. A tab-delimited, headerless file with the following columns: 
	
			0: Locus tag 
	
			1: Number of TA sites in the gene 
	
			2: Control sequencing reads in the gene 
	
			3: Experiment sequencing reads in the gene
	
			4: Mann Whitney U U value
	
			5: Mann Whitney U p-value 
	
			6: Experiment:Control read ratio 
	
			7: Length of the gene 
	
			8: Control reads per TA site 
	
			9: Experiment reads per TA site 

			10: Gene start 
	
			11: Gene stop  

		- FDR file: output from fdr.py script. A tab-delimited, headerless file with columns containing 0) locus tag and 1) FDR q-value 
		OUTPUT: a comma-delimited .txt file with no heading and the following columns: 
	
		0: Locus Tag
	
		1: Number of TA sites in gene 
	
		2: Length of gene 
	
		3: Control sequencing reads 
	
		4: Experiment sequencing reads 
	
		5: Experiment:Control reads ratio 
	
		6: FDR q-value 
	
		7: Modified read ratio - if either the control reads or experiment reads were less than 1*10-5 times the total read count for the file, the read count was increased to this minimum value before the ratio was calculate. 
	
		8: Index - reads per length.
	
		9: Fitness value - 	gene rank in range (0,1] with interval 1/# of genes based on the product of the index and the modified read ratio. 

	EX: >> python make_clust_file.py mwu_sample.txt fdr_sample.txt > clusteroutput.txt 

mannwhitneyu.py: python script that performs a Mann Whitney U test to compare two igv files. 
	REQUIRED INPUTS: 

		- orf list - a tab-delimited .txt file with no heading and the following columns: 
	
			0: Locus tag
	
			1: Gene start site 
	
			2: Gene stop site

		- A control .igv file with no heading and the following columns: 
	
			0: A reference ID for the contig
	
			1: TA start site 
	
			2: TA stop site 
	
			3: Number of sequencing reads 
	
			4: Locus tag, if TA site is within a gene 

		- The experment .igv file

	OUTPUT: A tab-delimited text file with no heading and the following columns: 
	
		0: Locus tag 
	
		1: Number of TA sites in the gene 
	
		2: Control sequencing reads in the gene 
	
		3: Experiment sequencing reads in the gene
	
		4: Mann Whitney U U value
	
		5: Mann Whitney U p-value 
	
		6: Experiment:Control read ratio 
	
		7: Length of the gene 
	
		8: Control reads per TA site 
	
		9: Experiment reads per TA site 
	
		10: Gene start 
	
		11: Gene stop  

	EX: >> python mannwhitneyu.py orfs.txt ctrl.igv exp.igv > output.txt

NCTC8325.fasta: Fasta for S. aureus NCTC8325 downloaded from NCBI

NCTC8325.genbank: Genbank for S. aureus NCTC8325 downloaded from NCBI

normalization.py: python script that creates a normalized control IGV file to compare an experimental IGV file to.
	REQUIRED INPUTS: two igv files. 

		- IGV file is a tab-delimited file with the following columns and no heading: 
	
			0. Reference ID for contig 
	
			1. TA start site 
	
			2. TA stop site 
	
			3. Locus tag of gene if TA site is located within a gene. 

	OUTPUT: Prints a normalized igv to replace the first input igv file. Normalized based on number of TA sites hit in the first igv compared to the second.  

	EX: >>python normalization.py file1.igv file2.igv > normalizedfile1.igv

promoter_analysis_driver.py: python script that performs promoter analysis on tabular files in a directory. 
	REQUIRED INPUTS:  
	1. A gene list tab-delimited text file with no header and the following fields: 
		0. Contig reference ID 
		1. Locus Tag
		2. Gene Start
		3. Gene Stop 
		4. Gene Name 
		5. Gene description 
	2. A TA site file. This can be downloaded from tagc.univ-mrs.fr/rsa-tools/. Tab-delimited text file, comment character ';', and the following fields: 
		0. Search ID  
		1. DNA strand searched ('D' for direct)
		2. Query sequence pattern ('TA')
		3. Contig reference ID 
		4. Sequence start site 
		5. Sequence stop site 
		6. Sequence found 
		7. Matching score 
	DIRECTORY MUST CONTAIN: 
	Tabular HopCount data files from Galaxy. Names should include treatment (for controls, treatment should be untr) and transposon construct ('blunt','cap,'dual','erm','pen','tuf') separated by '_'. Tabular files are tab-delimited text files with one heading line and the following fields: 
			0. Contig reference ID 
			1. Position 
			2. Locus tag 
			3. Gene name 
			4. Plus direction reads 
			5. Minus direction reads 
			6. Total reads 
			7. Gene product 
			8. Protein ID 
			9. Notes 
			10. Sequence 
	Outputs: 
	1. IGV files for each tabular file and a cpt tabular file that has the reads from all of the transposons with a single outward-facing promoter. IGV is a tab-delimited text file with no header and the following fields: 
		0. NCBI genome reference number
		1. Start of TA site
		2. End of TA site
		3. # Total Reads
		4. Locus Tag (if in gene)
		5. # Plus direction reads 
		6. # Minus direction reads 
	2. A promoter analysis text file (starts with PRO_) for each of the transposon files. This is a tab-delimited, headerless text file with the following fields: 
		0. Genomic window number 
		1. TA start
		2. Direction of hit 
		3. Log of the ratio of plus-direction control and experimental reads
		4. Log of the ratio of minus-direction control and experimental reads 
		5. Log of the ratio of plus direction and minus direction experimental reads
	EX: >> python promoter_driver.py genelist.txt TAsites.txt

sam_to_tabular.py: python script that converts SAM files from Bowtie alignment to tabular files that consolidate reads at the same site. 
	Directory must contain SAM files. 
	OUTPUT: Tabular file for each SAM file - the tabular file will not have annotations. Tabular is a tab-delimited with header line.  

saouhsc_orfs.txt: tab-delimited text file with no header with the columns: 
	0. locus tag
	1. gene start
	2. gene stop

tabular_to_igv.py: python script that converts tabular files downloaded from Tufts-Galaxy to IGV format. 
	DIRECTORY MUST CONTAIN: 
		- HG003_TAsites.txt
		- HG003_Genelist.txt
		- Tabular files you want to convert
	OUTPUT: IGV for each tabular file. Columns: 
		0. Reference ID
		1. TA start
		2. TA stop
		3. Total reads
		4. Locus tag (if TA site is in a gene)
		5. Plus reads
		6. Minus reads