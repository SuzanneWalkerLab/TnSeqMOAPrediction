# # PURPOSE: Convert tabular files downloaded from Tufts-Galaxy to IGV format. 
# # ----------

# # PACKAGES:
# # ----------
import glob 


# # Functions
# # ----------
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
	
def make_TAdict(tafile,genefile):
	"""
	Creates TA dictionary that can be used later. 
	
	REQUIRED INPUTS: 
	- tafile: a tab-delimited file with comment character ; and columns (the pattern is 'TA': 
		0. PatternID
		1. Strand pattern found in 
		2. Pattern searched for 
		3. Sequence ID 
		4. Start of Pattern 
		5. End of Pattern 
		6. Actual Sequence 
		7. Matching score
	- genelist: a tab-delimited, headerless file with columns: 
		0. Sequence reference ID 
		1. Locus tag 
		2. Gene start 
		3. Gene stop 
		4. Gene 
		5. Protein description 
	
	OUTPUT: 
	TA dictionary in the format {TAstart:[0,locus tag]}
	"""

	TAsites = {}
	# # Get TA start sites from TA site file 
	TAfile = open(tafile,'r')
	for line in TAfile:
		if not line.startswith(';'):
			info = line.split('\t')
			TAsites[str(int(info[4]))] = [0]
	TAfile.close()
	
	# # Label the TA sites that are found within a gene. 
	genelist = open(genefile,'r')
	for line in genelist:
		info = line.split('\t')
		locusTag = info[1]; start = int(info[2]); end = int(info[3])+1
		for i in range(start, end):
			if str(i) in TAsites: 
				TAsites[str(i)].append(locusTag)
	genelist.close()
	
	# # TA sites not found in a gene are labeled with an empty string. 
	for TA in TAsites:
		if len(TAsites[TA]) == 1: 
			TAsites[TA].append('')
	
	return TAsites

# # ----------
# # RUN
# # ----------
if __name__ == "__main__":
	files = glob.glob('*.tabular')
	tas = make_TAdict('HG003_TAsites.txt','HG003_GeneList.txt')
	for f in files: 
		igv(f,tas)