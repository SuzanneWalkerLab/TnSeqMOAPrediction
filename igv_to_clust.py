# #PURPOSE: Normalize experimental IGV file read counts to a control file, determine the fitness values of genes, and output the 'clust_*.csv' files that are used in the KNN algorithm. 
# # DEPENDENCIES: Need to have in same folder the following python scripts: 
# # - normalization.py
# # - mannwhitneyu.py 
# # - fdr.py
# # - make_clust_file.py
# # REQUIRED INPUTS: Folder must contain the following files: 
# # - allctrls.igv - a control IGV file to normalize to. Columns are: 
		# # 0. Reference ID for contig 
		# # 1. TA start site 
		# # 2. TA stop site 
		# # 3. Total reads 
		# # 4. Gene, if TA site is in a gene.
# # - As many experimental igvs as desired, labeled as treatment.igv
# # - saouhsc_orfs.txt - a tab-delimited text file with the columns: 
		# # 0. Locus tag 
		# # 1. Gene start site 
		# # 2. Gene stop site. 
# # OUTPUT: For each experimental igv file, outputs a clust_*.txt file. For more information on clust files, see the make_clust_file.py script. 
# # ----------

# # Packages
# # ----------
import os
import sys
from subprocess import call

# # Functions 
# # ----------
def without_extension(filename):
	# # Removes the extension of a filename. 
	return os.path.splitext(filename)[0]

def get_hit_sites(path):
	# # Determine the number of TA sites there are in an IGV file, the number of those TA sites that have a Tn insertion, and the number of total reads in the IGV file. 
	TA = 0; Reads = 0; HitSites = 0
	input_file = open(path)
	for line in input_file:
	    split = line.split()
	    if len(split) > 3:
	        if split[3][0].isdigit() == True:
	            TA += 1
	            Reads += float(split[3])
	            if float(split[3]) > 0: HitSites += 1
	input_file.close()
	return HitSites

def get_normname(filename):
	# # Get the name of the treatment from the normalized igv file names. <treatment>_<transposon>.<extension>
	name = without_extension(filename)
	return name.split('_')[0]

def get_clustname(filename):
	# # The input filenames are of the form <analysis>_<treatment>.<extension>. This returns the treatment. 
	name = without_extension(filename)
	return name.split('_')[-1]

# # Create normalized controls to compare experimental IGV files to 
# # ----------
# # Gather all igv files from folder. 
inputs = [path for path in os.listdir('.') if path.endswith('.igv')]
igv_files = []
# # Find the control file
for data in inputs:
	if without_extension(data) == "allctrls":
		ctrl = data
	else: 
		igv_files.append(data)

# # Create a normalized control for each experimental igv using normalization.py. 
for igv_file in igv_files:
	destination = without_extension(igv_file) + "_normedctrl.igv"
	os.system("python normalization.py %s %s > %s" % (ctrl, igv_file, destination))

# # Mann-Whitney U analysis (MWU)
# # ----------
# # Gather the normalized control files and add to list of igv files.  
mwuinputs = [path for path in os.listdir('.') if path.endswith('_normedctrl.igv')]
for normfiles in mwuinputs:
	igv_files.append(normfiles)

# # Group IGV files by treatment. This pairs the normalized control and experimenal files. 
group_labels = set(map(get_normname, igv_files))
grouped_items = []
for group_label in group_labels:
	group = []
	for igv_file in igv_files:
		noext = without_extension(igv_file)
		splitigv = noext.split('_')
		if group_label in splitigv[0]: group.append(igv_file)
	grouped_items.append(group)

# # Run the Mann-Whitney U analysis for the files using the mannwhitneyu.py script, comparing between the experimental igv and the normalized control igv.  
for groups in grouped_items:
	if groups[0].endswith('_normedctrl.igv'):
		destination = 'mwu_' + without_extension(groups[1]) + '.txt'
		os.system("python mannwhitneyu.py saouhsc_orfs.txt %s %s > %s" % (groups[0], groups[1], destination))
	else:
		destination = 'mwu_' + without_extension(groups[0]) + '.txt'
		os.system("python mannwhitneyu.py saouhsc_orfs.txt %s %s > %s" % (groups[1], groups[0], destination))

# # Perform false discovery rate analysis (FDR)
# # ----------
# # Find all of the MWU output files
fdrinputs = [path for path in os.listdir('.') if path.startswith('mwu_')]

# # For each file, perform Benjamini-Hochberg FDR using the fdr.py script
mwufiles = []
for fdrinput in fdrinputs:
 	destination = 'fdr_' + without_extension(fdrinput) + '.txt'
 	os.system("python fdr.py %s > %s" % (fdrinput, destination))

# # Create clust files 
# # ----------
# # Gather MWU and FDR output files
mwuclustinputs = [path for path in os.listdir('.') if path.startswith('mwu_')]
fdrclustinputs = [path for path in os.listdir('.') if path.startswith('fdr_')]

clustinputs = []
for i in mwuclustinputs:
	clustinputs.append(i)
for j in fdrclustinputs:
	clustinputs.append(j)

# # Group the MWU and FDR files by treatment
clust_labels = set(map(get_clustname, clustinputs))
grouped_items = []
for clust_label in clust_labels:
	clust = []
	for clustinput in clustinputs:
		noext = without_extension(clustinput)
		splitclust = noext.split('_')
		if clust_label in splitclust[-1]: clust.append(clustinput)
	grouped_items.append(clust)

# # Create the clust output using the make_clust_file.py script. 
for groups in grouped_items:
	destination = 'clust_' + get_clustname(groups[0]) + '.txt'
	os.system("python make_clust_file.py %s %s > %s" % (groups[0], groups[1], destination))

