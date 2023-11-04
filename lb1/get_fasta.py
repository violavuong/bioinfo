#!/usr/bin/env python3

import sys
	
def parse_fasta(fasta):
	'''This function takes as input a file and creates a dict which keys are the PDB IDs and whose values are the relative fasta seqs'''
	d = {}
	with open(fasta, 'r') as f:
		for line in f:
			if line[0] == '>':
				p_id = line.split('|')[1] #p_id, the protein id, is the key of the dictionary
			else:
				d[p_id] = d.get(p_id, '') + line.rstrip() #getting the fasta of each p_id
	return d


	

if __name__ == '__main__':
	ids = sys.argv[1]
	fasta = sys.argv[2]
	d =  parse_fasta(fasta)
	l_id = open(ids).read().rstrip().split('\n') #list of ids: reading the file, stripping the \n of each line, so that you obtain a single long line, and split it by the \n char
	for p_id in l_id:
		if p_id in d:
			print('>' + p_id)
			print(d[p_id])


