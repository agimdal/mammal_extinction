#!/usr/bin/env python 
# Created by Daniele Silvestro on 02/03/2012 => dsilvestro@senckenberg.de 
from numpy import *
import numpy as np
import sys, os, csv
from scipy.special import gdtr, gdtrix
from biopy.bayesianStats import hpd as calcHPD
np.set_printoptions(suppress=True) # prints floats, no scientific notation
np.set_printoptions(precision=3) # rounds all array elements to 3rd digit



#################### SIMULATION ##########################################
def simulate(L,M,root,s_species,scale=1.,maxSP=np.inf):
	ts=list()
	te=list()
	L,M,root=L/scale,M/scale,int(root*scale)
	
	for i in range(s_species): 
		ts.append(root)
		te.append(0)
	
	### ADD PROBABILITIES OF STATUS CHANGE
	for t in range(root,0): # time		
		TE=len(te)
		if TE>maxSP: break
		for j in range(TE): # extant lineages
			if te[j]==0:
				ran=np.random.random()
				#print ran, L,M[j]
				if ran<L: 
					te.append(0) # add species
					ts.append(t) # sp time
				elif ran < (L+M[j]): # extinction
					te[j]=t
	
	te=array(te)
	return -array(ts)/scale, -(te)/scale


############### SIMULATION SETTINGS ########################################
def write_to_file(f, o):
	sumfile = open(f , "wb") 
	sumfile.writelines(o)
	sumfile.close()

logfile = open("sim.txt" , "wb") 
wlog=csv.writer(logfile, delimiter='\t')



##########################################################################
###########                 SIMULATION SETTINGS                 ##########
##########################################################################




# size data set
minSP=5
maxSP=7000
n_reps = 1000 # number of simulations
scale=100.
time_range=50 # time range of future projection (years)
print_LTT=True

############ (REPLACE WITH REAL VALUES) ####################################
L0 = 0.00001   # speciation rate at the present
Mcat = np.array([0.00001, 0.2, 0.4]) # probability of extinction in 1 year (LC, VU, EN, ...) 
P_better = 0.01 # per-lineage probability of change in status (per year): Mcat -= 1
P_worse  = 0.05 # per-lineage probability of change in status (per year): Mcat += 1

# IUCN DATA
s_species=50 # no. of starting species
species_cat = np.random.randint(0,3,s_species) # random assignment to IUCN categories (REPLACE WITH REAL DATA)

M_species = Mcat[species_cat] # extinction probabilities (per year)




forward_sim=np.zeros((time_range,n_reps))	

for sim in range(n_reps):
	sys.stdout.write(".")
	sys.stdout.flush()
	i=0
	FA,LO=simulate(L0,M_species,-time_range,s_species,scale)
	ltt= "\nlineages through time (spindle diagram):"
	N=np.zeros(time_range)
	for i in range(int(max(FA))):
		n=len(FA[FA>i])-len(LO[LO>i])
		#nlog=int((n))
		ltt += "\n%s\t%s\t%s" % (i, n, "*"*n)
		N[i]=n
	
	forward_sim[:,sim]=N
	if print_LTT: print ltt
	for i in range(len(FA)):
		o+= "%s\t%s\t%s\n" % (i+1,FA[i],LO[i])



wlog.writerow(["year","mean","95m","95M","50m","50M","25m","25M"])

wlog.writerow([2015,s_species,s_species,s_species,s_species,s_species,s_species,s_species])


for i in range(root_r[0])[::-1]:
	hpd95 = CI(forward_sim[i,:],.95)
	hpd50 = CI(forward_sim[i,:],.50)
	hpd25 = CI(forward_sim[i,:],.25)
	meanN = mean(forward_sim[i,:])

	wlog.writerow([2015+(root_r[0]-i), meanN, hpd95[0],hpd95[1], hpd50[0],hpd50[1], hpd25[0],hpd25[1]])






quit()
