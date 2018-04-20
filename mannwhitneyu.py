# # PURPOSE: Performs a Mann Whitney U test to compare two igv files. 
# # REQUIRED INPUTS: 
# # - orf list - a tab-delimited .txt file with no heading and the following columns: 
	# # 0: Locus tag
	# # 1: Gene start site 
	# # 2: Gene stop site
# # - A control .igv file with no heading and the following columns: 
	# # 0: A reference ID for the contig
	# # 1: TA start site 
	# # 2: TA stop site 
	# # 3: Number of sequencing reads 
	# # 4: Locus tag, if TA site is within a gene 
# # - The experment .igv file
# # OUTPUT: A tab-delimited text file with no heading and the following columns: 
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
# # EX: >> python mannwhitneyu.py orfs.txt ctrl.igv exp.igv > output.txt
# # ----------

# # Packages
# # ----------
import sys, random, scipy, numpy
from math import *
from scipy import stats

# # Collect data 
# # ----------
# # Open orf table and put information in a dictionary. {locus tag: [start,stop]}
# # Also set up another dictionary which will be {locus tag: [[reads at each TA site for control],[reads at each TA site for experiment]]}

genes = {}; genereads = {}

for line in open(sys.argv[1]):
    split = line.split('\t')
    name = split[0]; start = int(split[1]); end = int(split[2])
    genes[name] = [start, end]
    genereads[name] = [[],[]]
  
genelist = genes.keys()
genelist.sort()

# # Obtain read counts for control file and add to dictionary, ignoring TA sites in the first or final 3% of a gene, as these are unlikely to have an effect on the gene. 
for line in open(sys.argv[2]):
    split = line.split()
    if len(split) > 4 and split[3][0].isdigit() == True and split[4] in genereads:
        site = int(split[1]); start = int(genes[split[4]][0]); end = int(genes[split[4]][1])
        length = end - start
        if (start+(0.03*length)) <= site <= (end-(0.03*length)):
            genereads[split[4]][0].append(float(split[3]))

# # Obtain read counts for control file and add to dictionary, ignoring TA sites in the first or final 3% of a gene, as these are unlikely to have an effect on the gene. 
for line in open(sys.argv[3]):
    split = line.split()
    if len(split) > 4 and split[3][0].isdigit() == True and split[4] in genereads:
        site = int(split[1]); start = int(genes[split[4]][0]); end = int(genes[split[4]][1])
        length = end - start
        if (start+(0.03*length)) <= site <= (end-(0.03*length)):
            genereads[split[4]][1].append(float(split[3]))

# # Perform Mann-Whitney U test and print to output file.  
# # ----------

for i in range(len(genelist)):
	# # Get numbers for each gene
    SAO = genelist[i]
    start = genes[SAO][0]; end = genes[SAO][1]
    Lib1Counts = genereads[SAO][0]
    Lib2Counts = genereads[SAO][1]
    TA = len(Lib1Counts)
    readslib1 = sum(Lib1Counts)
    readslib2 = sum(Lib2Counts)
    length = end - start
    if length > 0:
        index1 = readslib1/length
        index2 = readslib2/length
	# # If there are no TA sites, ignore the gene 
    if TA == 0:
        print "%s\t%d\t%s" % (SAO,TA,'Gene Has No TAs')
    else:
		# # If there are no reads, ignore the gene
        if sum(Lib1Counts) == 0 and sum(Lib2Counts) ==0:
            print "%s\t%d\t%s\t%s\t%s\t%s\t%s\t%d\t%s\t%s\t%r\t%r" % (SAO,TA,'0','0','0','0','0',length,'0','0',start,end)
		# # The gene SAOUHSC_01447 is very long and gives us a lot of trouble, so we ignore it 
        elif SAO == "SAOUHSC_01447":
            print "%s\t%d\t%d\t%d\t%d\t%s\t%0.10f\t%d\t%0.10f\t%0.10f\t%r\t%r" % (SAO,TA,readslib1,readslib2,U,'0.5',CountRatio,length,index1,index2,start,end)
		# # For all other genes, perform Mann-Whitney U analysis on the reads at each  TA site comparing the control and the experiment. 
        else:
            try:
                U, p_val = scipy.stats.mannwhitneyu(Lib1Counts,Lib2Counts)
                CountRatio = float(sum(Lib2Counts)+1)/float(sum(Lib1Counts)+1)
                print "%s\t%d\t%d\t%d\t%d\t%0.100f\t%0.10f\t%d\t%0.10f\t%0.10f\t%r\t%r" % (SAO,TA,readslib1,readslib2,U,p_val,CountRatio,length,index1,index2,start,end)
            except (ValueError):
                pass

      


        

