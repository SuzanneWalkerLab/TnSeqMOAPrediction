# # PURPOSE: Creates a list of informative genes to base KNN categorization on and then performs leave-one-out validation of KNN supervised machine-learning algorithm for predicting antibiotic mechanisms of action. Optionally can also predict the mechanism of action of a series of unknown antibiotics. 
# # REQUIRED INPUT: A tab-delimited txt file. Column 1 is a list of file names for the data files to be used in the training set. Column 2 is the numeric category the data file belongs in. Column 3 is empty. Column 4 is a list of categories. Column 5 is the numeric code (0 to [number of categories-1]) for the category.
# # The data files are comma-delimited files with the following columns and must be in the current directory: 
# # 0. Gene locus tag 
# # 1. Number of TA sites in the gene
# # 2. Gene length
# # 3. Number of control sample reads 
# # 4. Number of treatment sample reads 
# # 5. Raw treatment:control read ratio 
# # 6. FDR-corrected p-value
# # 7. Corrected ratio. For each the control and treatment files, the minimum read count for each gene is set to 1/10000 x the total read count of the file. The treatment:control read ratio is then taken. 
# # 8. Treatment sample reads per length of gene 
# # 9. Normalized rank - 7*8
# # OPTIONAL INPUT: A txt file listing the data files with unknown mechanism of action that we want to predict a category for. 
# # OUTPUTS: 
	# # - predictions.txt: a file listing the predictions for the leave-one-out validation for the training set and their respective probabilities. 
	# # - filesused.txt: a list of the files that were used for the analysis. Note that some input files are not used if they do not have enough genes that are sufficiently changed compared to the control or if too many of the genes are depleted of sequencing reads compared to the control. 
	# # - genematrix.txt: a binary matrix showing which files in the training set contributed which genes to the list of gene fitness values fed into the KNN algorithm. 
	# # - OPTIONAL: PRED_<name of unknowns list>.txt: a list of the predictions for the files with unknown mechanisms of action and their respective probabilities. 
# # EX: >>python KNN.py KNN_data_files.txt unknown_data_files.txt 

# # Packages 
# ----------
import glob
import sys
from sklearn import neighbors

# # Open output files 
# # ----------
out = open("predictions.txt","w") # lists the predicted categories for leave-one-out validation. 
out2 = open("filesused.txt","w") # Lists the files used in the analysis
out3 = open("genematrix.txt","w") # A binary matrix showing which files contributed which genes to the gene list. 

# # Retrieve information from inputs 
# # ----------
# # Get the possible list of files and their categories and store in dicts. 
unknownlist = None
if len(sys.argv)>2:
	unknownlist = sys.argv[2]
filedict = {}
genedict = {}
classdict = {}
filecats = open(sys.argv[1],"r")
for line in filecats: 
	info = line.split("\t")
	filedict[info[0]]=info[1].rstrip()
	if len(info)>2 and info[4][0].isdigit(): 
		classdict[info[4].rstrip()]=info[3]
classes = classdict.keys()
classes.sort()

# # Determine which files are suitable and which genes will be informative 
# # ----------
# # See if the files listed in input1 are actually good to use, and if they are, collect their gene lists: 
genes = []
files = glob.glob("clust_*.txt")
for file in files: 
	if file in filedict:
		dhits = 0
		# Gather the genes that are depleted or enriched under that treatment
		dep = []
		enr = []
		for line in open(file): 
			info = line.split(",")
			# # ctrl is the number of control reads, exp is the number of experimental reads, and ratio is the exp:ctrl ratio
			tag = info[0]; ctrl = int(info[3]); exp = int(info[4]); ratio = float(info[5]); fitness = float(info[-1].rstrip())
			if ctrl >= 100 or exp >= 100:
				# # Get a list of locus tags for the depleted genes. 
				if ratio <= 0.25: 
					dep.append(tuple([fitness,tag]))
				# # Get a list of locus tags for the control genes
				elif ratio >= 4:
					enr.append(tuple([fitness,tag]))
		# Check to make sure the treatment is not too harsh (more than 1000 depleted genes) or not harsh enough (fewer than 10 depleted genes).
		dhits = len(dep)
		if dhits <1000 and dhits>=10: 
			if dhits>10: 
				# # Sort hits by fitness value
				dep.sort()
				# # Only take top 10
				dep = dep[0:10]
			# # Similarly, only take top 15 most depleted genes. 
			if len(enr) > 15: 
				enr.sort()
				enr.reverse()
				enr = enr[0:15]
			# Catalog these genes as ones to use be input into the KNN analysis. 
			curgenes = dep+enr
			curgenes = [g[1] for g in curgenes]
			genes += curgenes
			for gene in curgenes: 
				if gene in genedict: 
					genedict[gene].append(file)
				else: 
					genedict[gene] = [file]
		# If the treatment is too harsh or not harsh enough, discard it. 
		else: 
			filedict.pop(file,None)

# # Perform leave-one-out validation
# # ----------
genes = list(set(genes))
fileset = filedict.keys()
fileset.sort()
for f in fileset:
	print f
	# # Leave one file out of training set. 
	trainingset = fileset[:]
	trainingset.remove(f)
	# # From the genes whose fitness is usually used in the KNN algorithm, remove the ones that are only in the gene set because of the treatment that was removed for this step of the leave-one-out validation.  
	traininggenes = genes[:]
	for g in genedict: 
		if genedict[g]==[f]: 
			traininggenes.remove(g)
	
	  
	categories = []
	trfitvals = []
	unkfitvals = []
	# # Gather the gene fitness values for the left out file.
	for line in open(f):
		split = line.split(',')
		if split[0] in traininggenes:
			unkfitvals.append(float(split[9].rstrip()))
	# # Gather the category indicators and gene fitness values for the training set files.
	for t in trainingset:
		fitvals = []
		categories.append(int(filedict[t]))
		for line in open(t):
			split = line.split(',')
			if split[0] in traininggenes:
				fitvals.append(float(split[9].rstrip()))
		trfitvals.append(fitvals)

	# # Set up the KNN algorithm
	neigh = neighbors.KNeighborsClassifier(n_neighbors=3, weights='distance', metric='minkowski', p=2)
	neigh.fit(trfitvals,categories)

	# # Perform the prediction for the left out file
	out.write(f+'\t'+filedict[f]+"\n")
	probs = neigh.predict_proba(unkfitvals)[0]
	probnclass = []
	for i in range(0,len(probs)): 
		pc = tuple([round(float(probs[i]),2),classes[i]])
		probnclass.append(pc)
	sortpc = sorted(probnclass)
	sortpc.reverse()
	# # Write to file the predicted categories and their probability values. 
	for prob in sortpc:
		if prob[0]>0: 
			out.write(str(prob[1])+"\t"+classdict[prob[1]]+"\t"+str(prob[0])+"\n")
            
out.close()

# # Write the files used and their categories to a second file
out2.write("File\tCategory\n")
for f in filedict: 
	out2.write(f+'\t'+str(filedict[f])+'\n')

# # Create a binary gene matrix showing which genes were selected from each file to be included in the analysis 
out3.write("Locus Tag")
for f in fileset: 
	out3.write("\t" + f)
out3.write("\n")
genesources = {}

for g in genedict:
	matrix = [0 for i in range(len(fileset))]
	for f in genedict[g]: 
		matrix[fileset.index(f)] = 1
	out3.write(g)
	for i in matrix: 
		out3.write("\t"+str(i))
	out3.write("\n")
	gfiles = list(set(genedict[g]))
	gfiles.sort()
	gfiles = tuple(gfiles)
	if gfiles in genesources: 
		genesources[gfiles].append(g)
	else: 
		genesources[gfiles]=[g]

# # Predict categories of unknowns 
# # ----------

# # Gather list of unknown data files 
if unknownlist:
	unkfiles = []
	out4 = open('PRED_'+unknownlist,"w")
	for line in open(unknownlist,'r'):
		sp = line.split('\t')
		if len(sp) == 1: 
			curunk = tuple([sp[0].rstrip(),'NA'])
		else: 
			curunk = tuple([sp[0],sp[1].rstrip()])
		unkfiles.append(curunk)
	
	# # Training set and gene list need to be updated because it is currently missing a file from leave-one-out validation...
	trainingset = fileset[:]
	traininggenes = genes[:]
	
	# # Gather the category indicator values and the gene fitness values for each training set file. 
	categories = []
	trfitvals = []
	for t in trainingset:
		fitvals = []
		categories.append(int(filedict[t]))
		for line in open(t):
			split = line.split(',')
			if split[0] in traininggenes:
				fitvals.append(float(split[9].rstrip()))
		trfitvals.append(fitvals)
	
	# # Set up the KNN algorithm
	neigh = neighbors.KNeighborsClassifier(n_neighbors=3, weights='distance', metric='minkowski', p=2)
	neigh.fit(trfitvals,categories)
	
	# # Get the fitvals for the genes in the unknown file. 
	for unkfile in unkfiles: 
		unkfitvals = []
		for line in open(unkfile[0]):
			split = line.split(',')
			if split[0] in genes:
				unkfitvals.append(float(split[9].rstrip()))
		
		# # Perform prediction and write output file
		out4.write(unkfile[0]+'\t'+unkfile[1]+"\n")
		probs = neigh.predict_proba(unkfitvals)[0]
		probnclass = []
		for i in range(0,len(probs)): 
			pc = tuple([round(float(probs[i]),2),classes[i]])
			probnclass.append(pc)
		sortpc = sorted(probnclass)
		sortpc.reverse()
		for prob in sortpc: 
			if prob[0]>0: 
				out4.write(str(prob[1])+"\t"+classdict[prob[1]]+"\t"+str(prob[0])+"\n")
            
	out4.close()
