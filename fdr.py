# # PURPOSE: Perform Benjamini-Hochberg FDR analysis on the output from the mannwhitneyu.py script.
# # REQUIRED INPUTS: mannwhitneyu.py output file - A tab-delimited text file with no heading and the following columns: 
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
# # OUTPUT: tab-delimited text file with no header and columns 0) locus tag and 1) FDR q-value
# # EX: >> python fdr.py mwu_sample.txt > fdr_sample_output.txt
# # ----------

# # Packages
# # ----------
import operator, random, sys

# # Functions
# # ----------
def bh_qvalues(pv):
    """
    Return Benjamini-Hochberg FDR q-values corresponding to p-values C{pv}.

    This function implements an algorithm equivalent to L{bh_rejected} but
    yields a list of 'adjusted p-values', allowing for rejection decisions
    based on any given threshold.

    @type pv: list
    @param pv: p-values from a multiple statistical test

    @rtype: list
    @return: adjusted p-values to be compared directly with the desired FDR
      level
    """
    if not pv:
        return []
    m = len(pv)
    args, pv = zip(*sorted(enumerate(pv), None, operator.itemgetter(1)))
    #print args
    qvalues = m * [0]
    mincoeff = pv[-1]
    qvalues[args[-1]] = mincoeff
    for j in xrange(m-2, -1, -1):
        coeff = m*pv[j]/float(j+1)
        if coeff < mincoeff:
            mincoeff = coeff
        qvalues[args[j]] = mincoeff
    return qvalues

# # Collect data
# # ----------

# # Create dictionary with genes and p-values from Mann-Whitney U output file  {gene:[p-value],...}
gene_dict = {}
for line in open(sys.argv[1]):
    split = line.split('\t')
    if len(split) >= 7:
        gene_dict[split[0]] = [split[5]]

# # Create ordered list of p_values with their corresponding genes

genes_list = []

for k,v in gene_dict.iteritems():
    genes_list.append([v[0],k])

genes_list.sort()

p_vals = []; genes = []

for i in genes_list:
    p_vals.append(float(i[0]))
    genes.append(i[1])

# # Perform analysis and print the results information
# # ----------

# # Assign q values for each p value
q_vals = bh_qvalues(p_vals)

# # Match q values to gene in gene dictionary. 
for i in range(len(p_vals)):
    if genes[i] in gene_dict:
        gene_dict[genes[i]].append(q_vals[i])

# # Print the results 
for line in open(sys.argv[1]):
    split = line.split('\t')
    print "%s\t" % split[0],
    if split[0] in gene_dict: print gene_dict[split[0]][1],
    print
