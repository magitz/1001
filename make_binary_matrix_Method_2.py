#!/usr/bin/env python

import argparse
import re

########################################################################################################################
#
#  make_binary_matrix.py
#
#  This script takes a DNA matrix as the input and outputs 4 files:
#	1) The binary matrix, optionally as phylip file: binary.outfile
#	2) The binary matrix as CSV file: binary.outfile.csv
#	3)		: outfile
#	4) The UIB file: uib.outfile 
#	Usage: python make_binary_matrix.py -i infile -o outfile -g outgroup -p y
#		-i: Input file, may be in phylip (non-interleaved) or simple text format.
#		-o: output base see above for full names.
#		-p: Phylip formated input: y or n: Default is y. Output will also be phylip formated.
#		-g: Outgroup taxon name, must be identical to name in file. If not supplied, last taxon is used.
########################################################################################################################
#
#	Written by: Matt Gitzendanner
#				University of Florida
#				Department of Biology
#				magitz@ufl.edu
#
#	License: This script is provided as-is with no guarantees. You are free to copy, modify or distribute, as long as
#			proper acknowledgement is provided and citated in publications.
#
# 	Versions:
#		1.0: Sept. 24, 2013
#
#
########################################################################################################################


#Parse commandline options.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input matrix file")
parser.add_argument("-o", help="Output matrix file")
parser.add_argument("-g", help="Outgroup taxon name, if not used, last taxon in dataset is used.")
parser.add_argument("-p",default='y', help="Input in phylip format (Note: Does not support interleaved data? [y] or n")

args = parser.parse_args()

in_file = args.i
out_file = args.o
outgroup = args.g
phylip= args.p

#Set output file names
out_bin_file="binary." + out_file
out_uib_file="uib." + out_file
out_csv_file="binary." + out_file + ".csv"

#Open the inputfile
try:
	IN = open (in_file, 'r')
except OSError as exception:
	pass


#Setup empty dictionaries to store sequences.
seq_dict={}
trans_seq_dict={}
trans_bin_dict={}
trans_bin_csv_dict={}
trans_uib_dict={}

seq_len=0

if phylip == "y" or phylip == None :
	IN.readline()   # skip the first line

for Line in IN :		#Read in_file into a dictionary with keys as taxa and values as sequences (all sequences are same length)
	Line = Line.strip('\n')
	Line_bits=re.split("\s+",Line)
	
	if Line_bits[0] == '' :		#Skip blank lines
		continue
	if seq_len == 0 :				#Set Length from first sequence
		seq_len = len(Line_bits[1])
		print "All input sequences should be %s bases long." %(seq_len)
	else :
		if seq_len != len(Line_bits[1]) :
			print "Error, %s has different sequence length (%d is not %d)" %(Line_bits[0], seq_len, len(Line_bits[1]) )	#if sequences are not the same exit.
			break
		
	seq_dict[Line_bits[0]]=Line_bits[1]
	trans_seq_dict[Line_bits[0]]=''			#Add the taxon names to all the dictionaries for output.
	trans_bin_dict[Line_bits[0]]=''
	trans_bin_csv_dict[Line_bits[0]]=''
	trans_uib_dict[Line_bits[0]]=''
	Last_taxon=Line_bits[0] #keep resetting until the true last taxon.

IN.close 

if outgroup == None:  #If outgroup was not specified on command line, set outgroup to Last_taxon.
	outgroup=Last_taxon
	
print "The outgroup taxon is %s." %(outgroup)

#Go through the sequence dictionary character by character, and make the output dictionaries from this.
base_count=1
while base_count <= seq_len :
	#Get sequence of outgroup
	out_base=seq_dict[outgroup][(base_count-1):base_count].upper()  #Get the base and convert to uppercase
	#print out_base
	#out_base=out_base[:base_count]
	
	#Go through rest of dictionary and get the bases present in other taxa.
	bases=[]
	for key in seq_dict :
		tax_base=seq_dict[key][(base_count-1):base_count].upper()
		if (tax_base != out_base) and (tax_base in ['G', 'A', 'T', 'C']) and (tax_base not in bases) :
			bases.append(tax_base)
		
	print "Character %d of outgroup is base %s, other bases are: %s" %(base_count,out_base,bases)
	
	#For each base that is present for this character comput the output as needed.
	for nuc in bases :
		#print "Looking at %s in char %d" %(nuc,base_count)
		#Go through rest of dictionary and get the transformed bases.
		for key in seq_dict :		# For each taxon...
			tax_base=seq_dict[key][(base_count-1):base_count].upper()	#Get the base
			if tax_base == out_base :
				trans_seq_dict[key]=str(trans_seq_dict[key]) + out_base		#Add base to growing transformed sequence
				trans_uib_dict[key]=str(trans_uib_dict[key]) + out_base
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "0"
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",0"
		
			elif tax_base == nuc :
				trans_seq_dict[key]=str(trans_seq_dict[key]) + nuc
				trans_uib_dict[key]=str(trans_uib_dict[key]) + nuc
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "1"
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",1"
				
			else :
				trans_seq_dict[key]=str(trans_seq_dict[key]) + "?"
				trans_uib_dict[key]=str(trans_uib_dict[key]) + "?"
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "?"
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",?"
	
	for key in seq_dict :			#Add extra column with Outgroup base or ? to the UIB dataset
		tax_base=seq_dict[key][(base_count-1):base_count].upper()
		if tax_base == out_base :
			trans_uib_dict[key]=str(trans_uib_dict[key]) + out_base
		else :
			trans_uib_dict[key]=str(trans_uib_dict[key]) + "?"
			
	base_count+=1
			

#to print everything in the same order, go back through the original file to get taxa to pull out of dictionary.
try:
	IN = open (in_file, 'r')
except OSError as exception:
	pass

#Open all the output files.
try:
	OUT_SEQ = open (out_file, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN = open (out_bin_file, 'w')
except OSError as exception:
	pass

try:
	OUT_CSV = open (out_csv_file, 'w')
except OSError as exception:
	pass
	
	
try:
	OUT_UIB = open (out_uib_file, 'w')
except OSError as exception:
	pass

if phylip == "y" or phylip == '':
	num_tax= len(trans_seq_dict)

	OUT_SEQ.write ("%d\t%d\n" %(num_tax,len(trans_seq_dict[Last_taxon])))
	OUT_BIN.write ("%d\t%d\n" %(num_tax,len(trans_bin_dict[Last_taxon])))
	OUT_UIB.write ("%d\t%d\n" %(num_tax,len(trans_uib_dict[Last_taxon])))
	#CSV file shouldn't have phylip header.
	
	IN.readline()   # skip the first line

for Line in IN :
	Line = Line.strip('\n')
	Line_bits=re.split("\s+",Line)
	if Line_bits[0] == '' :
		continue
	else : 
		OUT_SEQ.write ("%s\t%s\n" %(Line_bits[0],trans_seq_dict[Line_bits[0]]))
		OUT_BIN.write ("%s\t%s\n" %(Line_bits[0],trans_bin_dict[Line_bits[0]]))
		OUT_UIB.write ("%s\t%s\n" %(Line_bits[0],trans_uib_dict[Line_bits[0]]))
		OUT_CSV.write ("%s%s\n" %(Line_bits[0],trans_bin_csv_dict[Line_bits[0]]))

#OUT_SEQ.write ("Out\t%s\n" %(trans_seq_dict[Last_taxon]))
#OUT_BIN.write ("Out\t%s\n" %(trans_bin_dict[Last_taxon]))
#OUT_UIB.write ("Out\t%s\n" %(trans_uib_dict[Last_taxon]))
