# # PURPOSE: Convert SAM files to tabular files. 
# # Directory must have SAM files. Tabular files are tab-delimited with a header line. 

import glob 

sams = glob.glob('*.sam')
datadict = {}
for s in sams:
	for line in open(s,'r').readlines():
		if not line.startswith('@') and not line.startswith('/'):
			info = line.split('\t')
			flag = info[1]; ref = info[2]; site = int(info[3]); seq = info[9]
			if site == 28: 
				print line
			if flag == '0': 
				if (ref,site) in datadict: 
					datadict[(ref,site)][0]+=1
				else: datadict[(ref,site)]=[1,0]
			elif flag == '16':
				if (ref,site+15) in datadict: 
					datadict[(ref,site+15)][1]+=1
				else: 
					datadict[(ref,site+15)]=[0,1]
					
					
	output = open(s[:-3]+'tabular','w')
	output.write('Reference\tPosition\tLocus\tGene\tPlusCount\tMinusCount\tTotalCount\tProduct\tProteinID\tNote\tSequence\n')
	datakeys = datadict.keys()
	datakeys.sort(key=lambda x: x[1])
	for k in datakeys: 
		output.write(k[0]+'\t'+str(k[1])+'\t\t\t')
		output.write('%i\t%i\t%i\t\t\t\t\n' %(datadict[k][0],datadict[k][1],datadict[k][0]+datadict[k][1]))