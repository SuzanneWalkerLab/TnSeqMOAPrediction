'''
This performs promoter analysis on tabular files in a directory. 

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
EX: 
>> python promoter_driver.py genelist.txt TAsites.txt
'''

# # Setup
# # ----------
import glob
import sys
import time
import decimal as d 
import numpy as np

# # Functions
# # ----------
def combine_files(files,out=False):
	
	"""
	Adds the reads of multiple IGV or tabular files. 

	REQUIRED INPUTS: 
	files: a list of IGV or tabular files that you want to combine.

	OPTIONAL INPUTS: 
	out: name of an output file. Default is combined.igv or combined.tabular depending on the type of file input. 

	OUTPUT: 
	An IGV or tabular file with the added read counts. 
	"""
	# # Determine whether we're dealing with IGVs or Tabulars
	datatype = files[0].split('.')[-1]
	genes = {}
	# # For IGVs, we want to add the reads at each TA site in the genome. (IGV is a tab delimited file with NCBI Genome Reference #, TA Start, TA End, # of Reads, and Locus Tag if applicable). Make a dictionary to hold the information, sort it, and then save it to an output file. 
	if datatype == 'igv':
		for file in files:
			if file.endswith('igv'):
				currfile = open(file,'r')
				for line in currfile.readlines(): 
					info = line.split('\t')
					genRef = info[0]; start = int(info[1]); end = int(info[2]); reads = int(info[3].strip())
					if len(info)>4:
						locusTag = info[4].strip()
					else: locusTag = ''
					if (genRef,start) not in genes:
						genes[(genRef,start)] = [end,reads,locusTag]
					else: 
						genes[(genRef,start)][1] += reads
				currfile.close()
			else: 
				print "File types do not match. Please only use all IGVs or all tabular files."
				return
		if not out:
			output = open('combined.igv','w')
		else: 
			output = open(out,'w')
		keys = genes.keys()
		keys.sort()
		for key in keys:
			output.write(key[0]+'\t'+str(key[1])+'\t'+str(genes[key][0])+'\t'+str(genes[key][1])+'\t'+genes[key][2]+'\n')
		output.close()
	
	# # For tabulars, we want to add the plus direction, minus direction, and summed reads for each file. Tabulars are tab-delimited files downloaded from Galaxy after HopCount analysis. Again, store the information in a dictionary until you've gotten the numbers from all of the files, and then save to an output file. 
	elif datatype == 'tabular':
		for file in files: 
			if file.endswith('tabular'):
				currfile = open(file,'r')
				for line in currfile.readlines()[1:]: 
					info = line.split('\t')
					# # There can be some differences in the NCBI genome reference # format. Some measures are taken to simplify the number as much as possible. 
					if '|' in info[0]:
						genRefraw = info[0].split('|')[-2]
					else: genRefraw = info[0]
					if '.' in genRefraw:
						genRef = genRefraw.split('.')[-2]
					else: genRef = genRefraw
					pos = int(info[1]); plus = int(info[4]); minus = int(info[5]); total = int(info[6]); reads=[plus,minus,total]
					if (genRef,pos) not in genes: 
						genes[(genRef,pos)]=reads
					else: 
						for i in range(0,3):
							genes[(genRef,pos)][i]+=reads[i]
				currfile.close()						
			else:
				print "File types do not match. Please only use all IGVs or all tabular files."
				return
		if not out:
			output = open('combined.tabular','w')
		else: 
			output = open(out,'w')
		output.write('Reference\tPosition\tLocus\tGene\tPlusCount\tMinusCount\tTotalCount\tProduct\tProteinID\tNote\tSequence\n')
		keys = genes.keys()
		keys.sort()
		for key in keys:
			output.write(key[0]+'\t'+str(key[1])+'\t\t\t'+str(genes[key][0])+'\t'+str(genes[key][1])+'\t'+str(genes[key][2])+'\t\t\t\t\n')
		output.close()
				
	else:
		print "Data type not accepted. Please only input igv or tabular files."
		return
	
	return
def igv(tabular,TAdict,direction ='sum'): 
	""" 
	Makes an IGV file from a tabular file.  
	
	REQUIRED INPUTS: 
	tabular: a tabular HopCounts file downloaded from Galaxy with TnSeq data. 
	TAdict: A TA dictionary (see make_TAdict() output)
	
	OPTIONAL INPUTS: 
	direction: Can make an IGV file for just plus direction, just the minus direction, the plus and minus summed, or all three. Options: 'plus','minus','sum','all'
	
	OUTPUTS:
	Saves an IGV file with the same name as the tabular file, but with the .igv extension. (e.g. 'untr_blunt.tabular' makes an igv with the name 'untr_blunt.igv')
	The IGV file has the following tab-delimited data columns: 
	0. Genome reference number
	1. TA start site
	2. TA end site
	3. Number of reads (for direction='all' this will be total reads)
	4. Locus tag
	5. If 'all', plus direction reads (put here for backwards compatibility)
	6. If 'all', minus direction reads (put here for backwards compatibility)
	"""
	# # In the case of having all the read information (plus, minus, and total), we need two extra data slots
	if direction == 'all':
		for site in TAdict: 
			TAdict[site]+=[0,0]
	
	# # Loop through tabular file
	file = open(tabular,'r')
	for line in file: 
		info = line.split("\t")
		pos = 0 
	# # Make sure you're looking at a line with real data
		if len(info)>= 6 and info[1].isdigit():

		# # Look for reads in the positive direction. Because we sequence from the MmeI cut site, actual insertion site is some distance away. 
			if direction != 'minus':
				if int(info[4]) != 0:
					pos = int(info[1])
					if str(pos+15) in TAdict:
						TAdict[str(pos+15)][0] += int(info[4])
						if direction == 'all':
							TAdict[str(pos+15)][-2] += int(info[4])
					elif str(pos+14) in TAdict:
						TAdict[str(pos+14)][0] += int(info[4])
						if direction == 'all':
							TAdict[str(pos+14)][-2] += int(info[4])
		##Do the same with the reads that are in the reverse direction. 
			if direction != 'plus':
				if int(info[5]) != 0:
					pos = int(info[1])
					if str(pos-15) in TAdict: 
						TAdict[str(pos-15)][0] += int(info[5])
						if direction == 'all':
							TAdict[str(pos-15)][-1] += int(info[5])
					elif str(pos-16) in TAdict:
						TAdict[str(pos-16)][0] += int(info[5])
						if direction == 'all':
							TAdict[str(pos-16)][-1] += int(info[5])
	file.close()
	
	# # Create file names
	if direction=='sum' or direction=='all':
		output = open(tabular[:-7]+'igv','w')
	elif direction=='plus': 
		output = open(tabular[:-8]+'_plus.igv','w')
	elif direction=='minus':
		output = open(tabular[:-8]+'_minus.igv','w')
	else:
		print "Invalid direction. Direction must be 'plus','minus', 'sum', or 'all'."
		
	keys = TAdict.keys()
	keys = [int(k) for k in keys]
	keys.sort()

	## Everything gets written to the output IGV file.  
	for k in keys:
		output.write('SAOUHSC\t%d\t%d\t%d' %(int(k),int(k)+2,TAdict[str(k)][0]))
		if len(TAdict[str(k)]) > 1: 
			output.write('\t'+TAdict[str(k)][1])
			if direction == 'all':
				output.write('\t%d\t%d\n' %(TAdict[str(k)][-2],TAdict[str(k)][-1]))
			else: output.write('\n')
		else: 
			output.write('\t\n')
			
	# # The TAdict is reset so that it can be used again	
		TAdict[str(k)][0]=0
		TAdict[str(k)] = TAdict[str(k)][0:2]
	
	output.close()
	return
	
def make_dicts(genelistfile=None,talist=None,type='both'):
	"""
	Creates TA and gene dictionaries that can be used later. 
	
	REQUIRED INPUTS: 
	List of genes and List of TA sites (see beginning of script for descriptions of these files)
	
	OPTIONAL INPUTS: 
	type: They type of dictionary you want to make. Options are 'both', 'TA', or 'gene'. Default = 'both'
	
	OUTPUT: 
	If type = 'both', returns TA dictionary, Gene dictionary. 
	Otherwise, it just returns the dictionary you asked for.
	The TA dictionary is in the format {start:[0,locus tag]}. The gene dictionary is in the format {locus tag: [start,end,source,name,desc]}	
	"""
	# # Make TA dictionary. Keys: start. Values: [0,locus tag]. The 0 is a placeholder so that data for that TA site can later be filled. Start come from the TAsites list. The locus tag comes from the GeneList file. 
	if type == 'both' or type == 'TA':
		TAsites = {}
		TAfile = open(talist,'r')
		for line in TAfile:
			info = line.split('\t')
			if info[0]=='TA':
				TAsites[str(abs(int(info[4])))] = [0]
			else: 
				pass
		TAfile.close()
		
		genelist = open(genelistfile,'r')
		for line in genelist:
			info = line.split('\t')
			locusTag = info[1]; start = int(info[2]); end = int(info[3])+1
			for i in range(start, end):
				if str(i) in TAsites: 
					TAsites[str(i)].append(locusTag)
		genelist.close()
	
	for TA in TAsites:
		if len(TAsites[TA]) == 1: 
			TAsites[TA].append('')
	
	# # Make genes dictionary. Keys: Locus tags. Values: [start,end,source,name,description]. Start and end are the start and end of the gene. Source is the NCBI Genome Reference #, name is the name of the gene(Ex: lyrA), and description is the annotation for that gene (Ex: CAAX-protease like protein)
	if type == 'both' or type == 'gene': 
		genes = {}
		genelist = open(genelistfile,'r')
		for line in genelist:
			info = line.split('\t')
			source = info[0]; locusTag = info[1]; start = int(info[2]); end = int(info[3]); name = info[4]; desc = info[5].strip('\n')
			genes[locusTag] = [start, end, source, name, desc]
		genelist.close()
			
	if type=='both':
		return TAsites, genes
	elif type=='TA': 
		return TAsites
	elif type=='gene': 
		return genes
		
def promoter_analysis(files,stringency=False,windowsize=135,maxhits=150):
	"""
	Performs promoter analysis to identify windows of TA sites where there is an enrichment of reads in one direction. 
	
	REQUIRED INPUTS: 
	FILES: Can be a list of 4 files: [ctrlplus, ctrlminus,expplus,expminus], or 2 files: [ctrl,exp], depending on whether you created your IGV using the "plus" and "minus" directions or using the "all" direction (see igv() fucntion for more information)
	
	OPTIONAL INPUTS: 
	stringency: the number of standard deviations above the mean a read ratio should be to be considered a hit. Should be 4 or above. 
	windowsize: the size of the window to lump TA sites into. 
	
	OUTPUTS: 
	PRO_*.igv file that includes the NCBI genome reference number, the TA window, the start of the TA site, the direction of the transposon, the logs of the exp:ctrl ratios for transposons in the plus and minus direction, and the log of the plus:minus ratio for the experimental file. 
	Also returns the number of hits. 
	"""
	# # Set up decimal so that all of our math gives 4 decimal places
	d.getcontext().prec = 4
	
	if stringency:
		strin = stringency
	else: strin = 3
	
	# # Get data from the IGV files
	
	readdict = {}
	totalreads = [0,0,0,0]
	for file in files: 
		curfile = open(file,'r')
		for line in curfile.readlines():
			info = line.split('\t')
			if len(files) == 4:	
				totalreads[files.index(file)] += int(info[3])
				if (info[0],int(info[1])) not in readdict: 
					readdict[(info[0],int(info[1]))]=[info[4].rstrip('\n'),float(info[3])+1]
				else: 
					readdict[(info[0],int(info[1]))].append(float(info[3])+1)
			elif len(files) == 2:
				if files.index(file) == 0:
					totalreads[0] += int(info[-2])
					totalreads[1] += int(info[-1])
				else: 
					totalreads[2] += int(info[-2])
					totalreads[3] += int(info[-1])
				if (info[0],int(info[1])) not in readdict:
					readdict[(info[0],int(info[1]))]=[info[4].rstrip('\n'),float(info[-2])+1,float(info[-1].rstrip('\n'))+1]
				else: 
					readdict[(info[0],int(info[1]))].append(float(info[-2])+1)
					readdict[(info[0],int(info[1]))].append(float(info[-1].rstrip('\n'))+1)
			else: 
				print "Incorrect number of input files for promoter analysis. Must either input 4 files if the plus and minus direction reads for the control and experiment are stored in separate IGV files or two files if plus and minus reads are stored in a single IGV."
				return
	
	# # Create a correction factor in case one sample had more reads than another for some reason
	pluscorr = d.Decimal(totalreads[0])/totalreads[2]
	minuscorr = d.Decimal(totalreads[1])/totalreads[3]
	
	# # Calculate read ratios between ctrl and exp and between exp plus and minus
	count = 0
	datapoints = [[] for i in range(4)]
	totals = 4*[0]
	for TA in readdict:
		count += 1
		plusratio = pluscorr*d.Decimal(readdict[TA][3])/d.Decimal(readdict[TA][1])
		minusratio = minuscorr*d.Decimal(readdict[TA][4])/d.Decimal(readdict[TA][2])
		expratio = d.Decimal(readdict[TA][3])/d.Decimal(readdict[TA][4])
		plusratiolog = np.log10(plusratio)
		minusratiolog = np.log10(minusratio)
		expratiolog = np.log10(d.Decimal(expratio))
		variables = [plusratiolog,minusratiolog,expratiolog]
		for variable in variables:
			readdict[TA].append(variable)
			totals[variables.index(variable)]+=variable
			datapoints[variables.index(variable)].append(variable)
	
	# # Create cutoffs for determining what counts as a hit
	stdevs = []
	means = []
	cutoffs = []
	for i in range(len(datapoints)):
		means.append(totals[i]/count)
		stdevs.append(np.std(datapoints[i],axis=0))
		cutoffs.append(means[i]+strin*stdevs[i])
	# # Go through TAs, assign windows, and determine hits. 
	TAs = readdict.keys()
	TAs.sort()
	# # If a stringency is provided, use that to determine the cutoffs.
	if stringency: 
		hitdict = {}
		window = 1
		windowmax = 1+windowsize
		for TA in TAs: 
			if TA[1]>windowmax: 
				window+=1
				windowmax+=windowsize
			# # Identify hits in the plus direction
			if readdict[TA][5]>cutoffs[0] and readdict[TA][7]>cutoffs[2]: 
				hitdict[TA]='%s\t%d\t%d\tplus\t%f\t%f\t%f\n' %(TA[0],window,TA[1],readdict[TA][5],readdict[TA][6],readdict[TA][7])
			# # Identify hits in the minus direction
			elif readdict[TA][6]>cutoffs[1] and readdict[TA][7]<-1*cutoffs[2]:  
				hitdict[TA]='%s\t%d\t%d\tminus\t%f\t%f\t%f\n' %(TA[0],window,TA[1],readdict[TA][5],readdict[TA][6],readdict[TA][7])
			else: 
				pass
	# # If no stringency is provided, do the analysis with increasing stringency, starting at a stringency of 4, until you have fewer than 150 hits. 			
	else: 
		curhits = 5000
		stringency = 3
		hitdict = {}
		while curhits>maxhits:
			stringency += 1
			cutoffs = [i*stringency/(stringency-1) for i in cutoffs]
			curhits = 0
			hitdict = {}
			window = 1
			windowmax = 1+windowsize
			for TA in TAs: 
				if TA[1]>windowmax: 
					window+=1
					windowmax+=windowsize
				# # Identify hits in the plus direction
				if readdict[TA][5]>cutoffs[0] and readdict[TA][7]>cutoffs[2]: 
					hitdict[TA]='%s\t%d\t%d\tplus\t%f\t%f\t%f\n' %(TA[0],window,TA[1],readdict[TA][5],readdict[TA][6],readdict[TA][7])
				# # Identify hits in the minus direction
				elif readdict[TA][6]>cutoffs[1] and readdict[TA][7]<-1*cutoffs[2]:  
					hitdict[TA]='%s\t%d\t%d\tminus\t%f\t%f\t%f\n' %(TA[0],window,TA[1],readdict[TA][5],readdict[TA][6],readdict[TA][7])
				else: 
					pass
			curhits = len(hitdict)
	
	# # Print hits to output file
	hits = hitdict.keys()
	hits.sort()
	numhits = len(hits)
	if len(files)==4:
		outfile=open('PRO_'+files[2][:-9]+'.txt','w')
	else: 
		outfile=open('PRO_'+files[1][:-3]+'txt','w')
	outfile.write('GenRef#\tWindow\tTA Start\tDirection\tlog(+exp/+ctrl)\tlog(-exp/-ctrl)\tlog(+exp:-exp)\n')
	for hit in hits: 
		outfile.write(hitdict[hit])
	outfile.write('Stringency = %d\nWindow Size = %d' %(stringency,windowsize))
	outfile.close()
	return numhits
		
# # Get Ready!
# # ----------

# # Start the timer
now = time.strftime('%I:%M %p',time.localtime())
print "Started at %s." %(now)
starttime = time.clock()

# # Collect File Name Info
# # ----------

# # Find all the tabular files in the folder. Determine what treatments have been used, and what transposons are present. Expecting files to be named "<Treatment>_<transposon>.tabular", with control named 'untr_<transposon>.tabular
tabulars = glob.glob('*.tabular')
if len(tabulars)==0:
	print "No tabular files found in folder."
	sys.exit()
treatments = set()
promoters = set()
procombineds = set()
filelist = []

for tabular in tabulars: 
	info = tabular.split('_')
	if len(info) <2 or len(info)>2: 
		print "Tabular files are inappropriately named. Need to be either <treatment>_<promoter>.tabular"
		sys.exit()
	treatment = info[0]; promoter = info[1][:-8]
	currfile = [treatment, promoter]
	if promoter == 'cpt':
		procombineds.add(treatment)
	# # IMPORTANT VARIABLE: Filelist keeps track of all of the information acquired above and is used to sort files for different aspects of the driver
	filelist.append(currfile)

# # Combine Files
# # ----------
chrfiles = []
plasmids = set()
conditions = set()

# # Catalog the treatments
for file in filelist: 
	conditions.add(file[0])
		
# # Combine the cap, pen, and tuf promoter files that have the same treatment and add the combined file to the file list. 
for condition in conditions: 
	prolist = []
	if condition not in procombineds: 
		for file in filelist: 
			if file[1] in ['cap','pen','tuf']:
				if condition == file[0]:
					filename = '%s_%s.tabular' %(file[0],file[1])
					prolist.append(filename)
					proout = '%s_cpt.tabular'%(file[0])
		combine_files(prolist,proout)
		filelist.append([condition,'cpt'])

# # Progress update
time1 = time.clock()

print "Appropriate files have been combined. This step took %s seconds. Now making IGVs. This may take a few minutes." %(int(time1-starttime))

# # Make IGV Files 
# # ----------

curIGVs = glob.glob('*.igv')
# # For every organism, make a TAsite dictionary and a genes dictionary. The genes dictionary is going to be added to throughout the rest of the script and is then the basis for the output file. 
TAsites, genes = make_dicts(sys.argv[1],sys.argv[2])
# # Make the IGV files from blunt and combined tabular files and sort them into either the controls list or the experiments list based on the treatment information (control treatment MUST be 'untr').

igvs = []
controls = []
experiments = []
for file in filelist: 
	filename = '%s_%s.tabular' %(file[0],file[1])
	igvs.append(filename[:-7]+'igv')
	if file[0]=='untr': 
		controls.append(file)
	else: 
		experiments.append(file)
	if filename[:-7]+'igv' not in curIGVs:
		igv(filename,TAdict=TAsites,direction='all')
		
# # Progress update
time2 = time.clock()
print "IGVS created. This step took %s minute(s) and %s seconds. Now performing the promoter analysis." %(int(time2-time1)/60,int(time2-time1)%60)

# # Perform Analyses
# # ----------		
# # Find corresponding experiments and controls (files for which the only difference is the treatment)
for control in controls: 
	for experiment in experiments: 
		time3 = time.clock()
		if experiment[1]==control[1]:
			controlfile = '%s_%s.igv' %(control[0],control[1])
			experimentfile = '%s_%s.igv' %(experiment[0],experiment[1])
			
			# # Perform promoter analysis 
			promoter_analysis([controlfile,experimentfile])
			
			# # Progress update
			print "Analyses complete for %s. This step took %s minute(s) and %s seconds." %(experimentfile,int(time.clock()-time3)/60,int(time.clock()-time3)%60)

# # Final progress update  					
end = time.strftime('%I:%M %p',time.localtime())
endtime = time.clock()				
print "Completed at %s. Run time was approximately %d minutes." %(end,(endtime-starttime)/60)