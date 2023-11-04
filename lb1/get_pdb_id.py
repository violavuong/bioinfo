#!/usr/bin/env python3

import sys 

def get_pdb_id(f):
	with open(f, 'r') as fh:
		lines = fh.readlines()[1:] #reading all lines except the first (header)
		for line in lines: #for each line of all of them which are read
			line = line.rstrip() #strip them 
			entry = line.split(';') #and split them by ;
			new_entry = ''
			if len(entry[1]) > 1: #those PDB IDs which has more than one chain
				entry[1] = entry[1].split() #splitting the chain
				entries = ''
				for chain in entry[1]:
					new_entry += entry[0] + ' ' + chain + '\n' #assign to each ID a single chain 
				print(new_entry)
			else: #len(entry[1]) == 1 meaning that there is only a chain
				new_entry = entry[0] + ' ' + entry[1] + '\n'
				print(new_entry)
	


if __name__ == '__main__':
	f = sys.argv[1]
	print(get_pdb_id(f))	
