#!/usr/bin/python2.7
import os
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
		'--output',
		required=True,
		action=CompletePath,
		default=None,
		help='The name of the output folder'
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


csv = args.csv
output = args.output
time_range = args.years


with open(csv) as f:
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


# get the beginning and end year difference and the current category
