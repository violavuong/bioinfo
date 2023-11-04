#!/usr/bin/env python3

import sys
import numpy as np

def get_cm(f, t = 1.0):
	'''function that build a confusion matrix'''
	cm = np.zeros((2, 2)) #defining an 2x2 empty matrix
	#to compute the matrix, we need TP, TN, FP, FN
	#listed as 0s: negatives (TN + FN); listed as 1s: positives (TP + FP)
	#tuple (0,0): TN
	#tuple (1,0): FP
	#tuple (0,1): FN
	#tuple (1,1): TP
	with open (f, 'r') as fh:
		for line in fh:
			v = line.split()
			pos_col = int(v[-1]) #the class of the ID is in the last column
			e_val = float(v[1]) #the e-value is in the second colum (first for list)
			pos_row = 0
			if e_val < t: pos_row = 1 
			cm[pos_row, pos_col] += 1 
	return cm		
			
def get_accuracy(cm):
	'''function which returns the accuracy of the confusion matrix'''
	return (cm[0,0] + cm[1, 1]) / np.sum(cm) #cm[0, 0] are the TN; cm[1, 1] are the TP

def get_mcc(cm):
	n = cm[0, 0] * cm[1, 1] - cm[0, 1] * cm[1, 0] #numerator of the MCC eq (same tuples above)
	d = np.sqrt((cm[0, 0] + cm[1, 1]) * (cm[0, 0] + cm[1, 0]) * (cm[1, 1] + cm[0, 1]) * (cm[1, 1] + cm[1, 0])) #denominator of MCC
	return n/d

if __name__ == '__main__':
	f = sys.argv[1] #file containing all info 
	t = float(sys.argv[2]) #threshold
	cm = get_cm(f, t)
	acc = get_accuracy(cm)
	print('TH:', t)
	print('ACC:', acc)
	print('MCC:', get_mcc(cm))
	print('TN:', cm[0, 0], 'FN:', cm[0, 1])
	print('TP:', cm[1, 1], 'FP:', cm[1, 0])
	print()
