# # PURPOSE: Create a normalized control IGV file to compare an experimental IGV file to. 
# # REQUIRED INPUTS: two igv files. 
# # - IGV file is a tab-delimited file with the following columns and no heading: 
	# # 0. Reference ID for contig 
	# # 1. TA start site 
	# # 2. TA stop site 
	# # 3. Locus tag of gene if TA site is located within a gene. 
# # OUTPUT: Prints a normalized igv to replace the first input igv file. Normalized based on number of TA sites hit in the first igv compared to the second.  
#EX: >>python normalization.py file1.igv file2.igv >normalizedfile1.igv
# # ----------

# # Packages
# # ----------
import sys, numpy
from math import *
from numpy import *

# # Normalization
# # ----------
# # Set up variables
TActrl = 0; Readsctrl = 0; HitSitesctrl = 0
inputreads = []
saouhsc= []
start = []
end = []
window = []

# # Count TA sites in genome, sequencing reads, and number of TA sites with insertions in the first file. 
for line in open(sys.argv[1]):
    split=line.split('\t')
    if len(split) > 4:
        if split[3][0].isdigit() == True:
            inputreads.append(float(split[3]))
            saouhsc.append(split[0])
            start.append(split[1])
            end.append(split[2])
            window.append(split[4].rstrip('\n'))
            TActrl += 1
            Readsctrl += float(split[3])
            if float(split[3]) > 0: HitSitesctrl += 1

# # Count TA sites in genome, sequencing reads, and number of TA sites with insertions in the second file. 
TAexp = 0; Readsexp = 0; HitSitesexp = 0
for line in open(sys.argv[2]):
    split=line.split()
    if len(split) > 3:
        if split[3][0].isdigit() == True:
            TAexp += 1
            Readsexp += float(split[3])
            if float(split[3]) > 0: HitSitesexp += 1

# # Calculate the ratio between the number of TA sites hit in the experiment and the control. 
TAproportion =  float(HitSitesexp)/HitSitesctrl   

# # Read counts in the first igv are divided by the total number of reads and then multiplied by the TAs hit proportion. This represents the probability of finding a read at the site. 
inputproportion = []
for i in inputreads:
	inputproportion.append(float(i)/Readsctrl)
inputproportiontanorm = []
for i in inputproportion:
    inputproportiontanorm.append(float(i)*TAproportion)

# # Append a number so that the total of inputproportiontanorm sums to a probability of 1
inputproportiontanorm.append(TAproportion-1)

# # Use a multinomial distribution to simulate where we would expect reads to be in the first file if the number of TA sites hit were the same for both file 1 and file 2. Peform simulation 100 times. 
multinominputsample = numpy.random.multinomial(Readsexp,inputproportiontanorm,100)

# # Calculate the number of reads in the last simulation 
multinominputsample = multinominputsample.T
for i in range(len(multinominputsample[0])):
    summulti = sum(multinominputsample[:,i])
# # Calculate a correction factor based on how many reads were put in real TA sites vs put in the bin that represents the leftover probability in inputproportiontanorm. 
multisum = numpy.repeat(summulti,100)
difference=numpy.repeat(multisum-multinominputsample[-1,0],100)
correctionfactor = round(float(Readsexp)/difference[0],4)
correctedinput = [i*correctionfactor for i in multinominputsample]
bootstrapcontrol = numpy.delete(correctedinput, (-1), axis=0)
# # Calculate the average number of simulated reads at each TA site for the 100 simulations. 
avgbootstrapcontrol = bootstrapcontrol.mean(axis=1)

# # Print a new IGV file with the normalized read counts for first file. 
for x in range(0,len(avgbootstrapcontrol)):
    print saouhsc[x], '\t', start[x], '\t', end[x], '\t', avgbootstrapcontrol[x], '\t', window[x]


