# Created by Daniele Silvestro on 02/03/2012 => dsilvestro@senckenberg.de
from numpy import *
import numpy as np
import sys, os, csv
from scipy.special import gdtr, gdtrix
#from biopy.bayesianStats import hpd as calcHPD
np.set_printoptions(suppress=True) # prints floats, no scientific notation
np.set_printoptions(precision=3) # rounds all array elements to 3rd digit



#************************TOBI'S PART**********************
#written by Tobias Hofmann (tobias.hofmann@bioenv.gu.se)

import argparse


#XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
#%%% Input %%%


# Complete path function
class CompletePath(argparse.Action):
    """give the full path of an input file/folder"""
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, os.path.abspath(os.path.expanduser(values)))


# Get arguments
def get_args():
	parser = argparse.ArgumentParser(
		description="Simulate the extinction of IUCN categorized taxa in the future",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument(
		'--csv',
		required=True,
		action=CompletePath,
		default=None,
		help='The csv file containing the IUCN information'
	)
	parser.add_argument(
		'--speciation_rate',
        type=int,
        default=0.00001,
		help='Speciation rate of the input clade, as estimated by e.g. PyRate'
	)
	parser.add_argument(
		'--years',
        type=int,
        default=50,
		help='Number of years that you want to simulate'
	)
	return parser.parse_args()

# Preparation for calling input variables and files
args = get_args()


csv_file = args.csv
time_range = args.years


with open(csv_file) as f:
    content = [x.strip('\n') for x in f.readlines()]

header = content[0]
head_length = len(header.split('\t'))
body = content[1:]

species_list = []
species_cat = []
total_up_count = 0
total_down_count = 0

for line in body:
    elements = line.split('\t')
    dif = head_length-len(elements)
    # add potentially missing columns to end of each line
    for i in range(dif):
        elements.append('')
    statuses = elements[1:]
    iucn_code = {'LC':0, 'CD':1, 'NT':2, 'VU':3, 'EN':4, 'CR':5, 'EW':6, 'EX':7}
    print "\n", elements[0]
    species_list.append(elements[0])
    #print statuses
    list_statuses = []
    for i in range(len(statuses)):
        status = ''
        if statuses[i] != '':
            status = statuses[i]
            if i+1 >= len(statuses):
                pass
            elif statuses[i+1] == '':
                statuses[i+1] = status
        if status in iucn_code:
            print status, '-->', iucn_code[status]
            new_status = iucn_code[status]
            list_statuses.append(new_status)
    print "Current status:", status, "(%s)" %new_status
    species_cat.append(new_status)
    up_count = 0
    down_count = 0
    for n in range(len(list_statuses)):
        if n+1 >= len(list_statuses):
            pass
        elif list_statuses[n] > list_statuses[n+1]:
            down_difference = list_statuses[n]-list_statuses[n+1]
            down_count = down_count + down_difference
        elif list_statuses[n] < list_statuses[n+1]:
            up_difference = list_statuses[n+1] - list_statuses[n]
            up_count = up_count + up_difference
    print "Count of rate shifts up (more severe protection status):", up_count
    total_up_count = total_up_count + up_count
    print "Count of rate shifts down (relaxation of protection status):", down_count
    total_down_count = total_down_count + down_count
# get the number of years that are recorded in the data matrix
rec_years = int(header.split('\t')[len(header.split('\t'))-1])-int(header.split('\t')[1])
lineages_count = len(species_list)

# Calculate rates of shifts per lineage and year
rate_shifts_up = float(total_up_count)/(float(rec_years)* float(lineages_count))
rate_shifts_down = float(total_down_count)/(float(rec_years)* float(lineages_count))



print "\n\n", "=" * 50, "\nTotal statistics:\n"
print "Recorded years:", rec_years
print "Number of lineages:", lineages_count
print "Category-shifts up:", total_up_count
print "Category-shifts down", total_down_count
print "The rate of category-shifts up (per year and lineage)", rate_shifts_up
print "The rate of category-shifts down (per year and lineage)", rate_shifts_down
print "=" * 50


#******************************************************





#################### SIMULATION ##########################################
def simulate(L,Mcat,species_cat_arg,root,s_species,P_better,P_worse,scale=1.,maxSP=np.inf):
	ts=list()
	te=list()
	species_cat = species_cat_arg+0
	M = Mcat[species_cat]
	L,M,root=L/scale,M/scale,int(root*scale)
	P_better,P_worse = P_better/scale,P_worse/scale
	#print P_better, P_worse

	for i in range(s_species):
		ts.append(root)
		te.append(0)

	### ADD PROBABILITIES OF STATUS CHANGE
	for t in range(root,0): # time
		TE=len(te)
		if TE>maxSP: break
		for j in range(TE):
			ran=np.random.random()
			if ran < P_better and species_cat[j]>0: # improve status
				species_cat[j]-=1
				M[j] = Mcat[species_cat[j]-1]

			elif ran < (P_better+P_worse) and species_cat[j]<max(species_cat): # worsen status
				species_cat[j]+=1
				M[j] = Mcat[species_cat[j]]


			if te[j]==0:
				ran=np.random.random()
				#print ran, L,M[j]
				if ran<L:
					te.append(0) # add species
					ts.append(t) # sp time
				elif ran < (L+M[j]): # extinction
					te[j]=t

	te=array(te)
	return -array(ts)/scale, -(te)/scale, species_cat


############### SIMULATION SETTINGS ########################################
def write_to_file(f, o):
	sumfile = open(f , "wb")
	sumfile.writelines(o)
	sumfile.close()

logfile = open("projected_diversity.txt" , "wb")
wlog=csv.writer(logfile, delimiter='\t')

tste_file = open("tste.txt" , "wb")
tste_log=csv.writer(tste_file, delimiter='\t')


##########################################################################
###########                 SIMULATION SETTINGS                 ##########
##########################################################################




# size data set
minSP=5
maxSP=7000
n_reps = int(time_range) # number of simulations
scale=100.
print_LTT=False

############ (REPLACE WITH REAL VALUES) ####################################
L0 = args.speciation_rate   # speciation rate at the present (from pyrate output or estimated from phylogeny etc.)
Mcat = np.array([0.0, 0.00001, 0.0001, 0.00105305, 0.011095167, 0.066967008, 0.1, 1]) #!!!IMPORTANT!!!!:values 0,1,2,6 and 7 need review
# probability of extinction in 1 year (LC, VU, EN, ...)
P_better = rate_shifts_down # per-lineage probability of change in status (per year): Mcat -= 1
P_worse  = rate_shifts_up # per-lineage probability of change in status (per year): Mcat += 1

# IUCN DATA
s_species = lineages_count
print s_species
# no. of starting species (no. of species present today)
#M_species = Mcat[species_cat] # extinction probabilities (per year)


forward_sim=np.zeros((time_range,n_reps))
print len(forward_sim[0])
species_status=np.zeros((n_reps,max(species_cat)+1))
ts_te_array = np.zeros((time_range,n_reps*2+2))
print len(ts_te_array[:,1])
ts_te_array[:,1]=np.arange(s_species)

for sim in range(n_reps):
	sys.stdout.write(".")
	sys.stdout.flush()
	i=0
	FA,LO,species_cat_mod=simulate(L0,Mcat,species_cat,-time_range,s_species,P_better,P_worse,scale)
	ltt= "\nlineages through time (spindle diagram):"
	N=np.zeros(time_range)
	for i in range(int(max(FA))):
		n=len(FA[FA>i])-len(LO[LO>i])
		#nlog=int((n))
		ltt += "\n%s\t%s\t%s" % (i, n, "*"*n)
		N[i]=n

	forward_sim[:,sim]=N
	sp_cat_hist = np.array([len(species_cat_mod[species_cat_mod==j]) for j in range(0,max(species_cat)+1)])
	species_status[sim,:]=sp_cat_hist
	if print_LTT:
		print ltt
		print species_cat_mod
	#for i in range(len(FA)):
		#ts_te+= "0\t%s\t%s\t%s\n" % (i+1,FA[i],LO[i])
	ts_te_array[:,2+2*sim]=FA
	ts_te_array[:,3+2*sim]=LO



wlog.writerow(["year","mean","95m","95M","50m","50M","25m","25M"])

wlog.writerow([2015,s_species,s_species,s_species,s_species,s_species,s_species,s_species])

def CI(x,p):
	x= sort(x)
	dp = int(len(x)*((1- p)/2))
	m = x[dp]
	M = x[len(x)-1-dp]
	return [m, M]

for i in range(time_range)[::-1]:
	hpd95 = CI(forward_sim[i,:],.95)
	hpd50 = CI(forward_sim[i,:],.50)
	hpd25 = CI(forward_sim[i,:],.25)
	meanN = mean(forward_sim[i,:])
	wlog.writerow([2015+(time_range-i), meanN, hpd95[0],hpd95[1], hpd50[0],hpd50[1], hpd25[0],hpd25[1]])



tste_head="clade\tspecies"
for i in range(s_species): tste_head+= "\tts\tte"
tste_log.writerow(tste_head.split("\t"))

for i in range(s_species):
	S = list(ts_te_array[i,:])
	tste_log.writerow(S)



quit()
