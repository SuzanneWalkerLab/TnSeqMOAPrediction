# # PURPOSE: Creates the input for our KNN algorithm, known as a clust file, converting read ratios into fitness values. 
# # REQUIRED INPUTS:  
# # - Mann-Whitney U file: output from mannwhitneyu.py script. A tab-delimited, headerless file with the following columns: 
	# # 0: Locus tag 
	# # 1: Number of TA sites in the gene 
	# # 2: Control sequencing reads in the gene 
	# # 3: Experiment sequencing reads in the gene
	# # 4: Mann Whitney U U value
	# # 5: Mann Whitney U p-value 
	# # 6: Experiment:Control read ratio 
	# # 7: Length of the gene 
	# # 8: Control reads per TA site 
	# # 9: Experiment reads per TA site 
	# # 10: Gene start 
	# # 11: Gene stop  
# # - FDR file: output from fdr.py script. A tab-delimited, headerless file with columns containing 0) locus tag and 1) FDR q-value 
# # OUTPUT: a comma-delimited .txt file with no heading and the following columns: 
	# # 0: Locus Tag
	# # 1: Number of TA sites in gene 
	# # 2: Length of gene 
	# # 3: Control sequencing reads 
	# # 4: Experiment sequencing reads 
	# # 5: Experiment:Control reads ratio 
	# # 6: FDR q-value 
	# # 7: Modified read ratio - if either the control reads or experiment reads were less than 1*10-5 times the total read count for the file, the read count was increased to this minimum value before the ratio was calculate. 
	# # 8: Index - reads per length.
	# # 9: Fitness value - 	gene rank in range (0,1] with interval 1/# of genes based on the product of the index and the modified read ratio. 
# # EX: >> python make_clust_file.py mwu_sample.txt fdr_sample.txt > clusteroutput.txt
# # ----------

# # Packages
# # ----------
import operator, random, sys

# # Perform calculations
# # ----------
# # Create Dictionary with genes to store all of the information needed to make the clust file. 
gene_dict = {}

# # Get total read counts for control and experiment 
readcountctrl = 0
readcountexp = 0
for line in open(sys.argv[1]):
    split = line.split('\t')
    if split[0][0:3]=='SAO':
        try:
            readcountctrl = readcountctrl + int(split[2])
            readcountexp = readcountexp + int(split[3])
        except Exception:
            pass

# # Create a minimum read count for the control and experiment genes
mincountctrl = float(readcountctrl)/10000
mincountexp = float(readcountexp)/10000

# # Put important info in gene dictionary: 0-#TA sites, 1-length, 2-readsctrl,3-readsexp, 4-rawratio, 5-modratio calculated using mincountctrl and mincountexp as the lowest possible read count, 6-index: experiment reads per length, 7-rawfitscore: modratio*index
for line in open(sys.argv[1]):
    split = line.split('\t')
    if split[0][0:3]=='SAO' and int(split[1])>=0:
        try:
            readsctrl=int(split[2])
        except Exception:
            pass
		# # If the number of reads is less than the minimum set above, reset the number of reads to the minimum 
        if readsctrl < mincountctrl: readsctrl = mincountctrl
        try:
            readsexp = int(split[3])
        except Exception:
            pass
        if readsexp <mincountexp: readsexp = mincountexp
        ratio = float(readsexp)/readsctrl
        try:
            TAsites = int(split[1])
            length = int(split[7])
            rawratio = float(split[6])
            index = float(readsexp)/length
            rawfit = ratio*index
            gene_dict[split[0]] = [TAsites,length,int(split[2]),int(split[3]),rawratio,ratio,index,rawfit]
        except Exception:
            pass
        

#Add corrected pval from fdr file into dictionary as the 9th element (index 8)
for line in open(sys.argv[2]):
    split = line.split('\t')
    if split[0] in gene_dict:
        gene_dict[split[0]].append(float(split[1].rstrip()))

#Make fitvals a list instead of a dict for calculating norm-fit value
fitvals = []
for k,v in gene_dict.iteritems():
    fitvals.append([v[7],k])

# # Rank each raw fitness value in the range (0,1] with intervals of 1/# of genes to get the final fitness value. 
fitvals.sort()
count=0
for f in fitvals:
    f.append(count)
    count+=1
    normorder = float(count)/len(fitvals)
    f.append(normorder)


# # Append fitness value to gene dictionary values
for f in fitvals:
    gene = f[1]
    if gene in gene_dict:
        gene_dict[gene].append(f[3])

# # Print output clust file 
# ----------

genes = gene_dict.keys()
genes.sort()

for g in genes:
    if len(gene_dict[g]) > 6:
        print "%s,%d,%d,%d,%d,%0.10f,%0.100f,%0.10f,%0.10f,%0.6f" % (g,gene_dict[g][0],gene_dict[g][1],gene_dict[g][2],gene_dict[g][3],gene_dict[g][4],gene_dict[g][8],gene_dict[g][5],gene_dict[g][6],gene_dict[g][9])

 