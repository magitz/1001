#!/usr/bin/env python

import argparse
import re

########################################################################################################################
#
#  1001.py
#
#  This script takes a character matrix as the input and outputs 8 files:
#	1,2) The binary matrix with all characters as phylip and csv
#	3,4) A reduced version of the matrix with invariant characters removed as phylip and csv
#	5,6) A polarized version polarized with respect to outgroup as phylip and csv.
#			For polarizing: Outgroup=0; if character in taxon is same as outgroup then 0, else 1. 
#	7,8) Same as 5&6, but removing invariant characters.
# 
#	Usage: python 1001.py -i infile -o outfile -g outgroup -c DNA 
#		-i: Input file, may be in phylip (non-interleaved) or simple text format.
#		-c: Matrix character type (DNA, AA, or MULTI). Default is DNA
#		-a: Characters used for ambiguous or missing data
#			* If -c DNA, the assumption is -a RYKMSWNBDHV? if others are used, use -a with list of characters
#			* If -c AA, the assumption is -a X? if others are used, use -a with list of characters
#			* If -c MULTI, the assumption is only -a ? again, if others are used, use -a with list of characters
#		-o: output base see above for full names.
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
#		2.0: Dec. 18, 2013: Modified to accept AA and multi-state character datasets and changed outputs a bit.
#		2.1: Dec. 19, 2013: Bug fixes in handling of characters.
#		2.2: Jan. 1, 2014: Bug fix in polarizing characters.
#		2.3: Feb 3, 2015: First public release, renamed to 1001.py
#
########################################################################################################################


#Parse commandline options.
parser = argparse.ArgumentParser()
parser.add_argument("-i", help="Input matrix file")
parser.add_argument("-o", help="Output matrix file")
parser.add_argument("-c",default="DNA", help="Matrix character type (DNA, AA, MULTI)")
parser.add_argument("-a", help="Characters used for ambiguous or missing data")
parser.add_argument("-g", help="Outgroup taxon name, if not used, last taxon in dataset is used.")

args = parser.parse_args()

in_file = args.i
out_file = args.o
outgroup = args.g
data_type=args.c
ambig_chars=args.a


#Set output file names
out_bin_all = out_file + ".binary.all.phy"
out_bin_all_csv = out_file + ".binary.all.csv"
out_bin_reduced = out_file + ".binary.reduced.phy"
out_bin_reduced_csv = out_file + ".binary.reduced.csv"
out_bin_polarized = out_file + ".binary.polarized.phy"
out_bin_polarized_csv = out_file + ".binary.polarized.csv"
out_bin_polarized_reduced = out_file + ".binary.polarized.reduced.phy"
out_bin_polarized_reduced_csv = out_file + ".binary.polarized.reduced.csv"

#Open the inputfile
try:
	IN = open (in_file, 'r')
except OSError as exception:
	pass


#Setup empty dictionaries to store sequences.
char_dict={}
trans_char_dict={}
trans_bin_dict={}
trans_bin_csv_dict={}
trans_bin_reduced_dict={}
trans_bin_reduced_csv_dict={}
trans_polarized_dict={}
trans_polarized_csv_dict={}
trans_polarized_reduced_dict={}
trans_polarized_reduced_csv_dict={}


#Report Matrix type and (assumed) ambiguity codes

if data_type == "DNA" :
	print "Matrix data type is: DNA"
	if ambig_chars == None :
		ambig_chars = list("RYKMSWNBDHV?")
		print "No set of ambiguous characters was passes with -a option, so assuming %s." %(ambig_chars)
	else :
		ambig_chars=list(ambig_chars)
		print "The following characters are treated as ambiguous or missing: %s " %(ambig_chars)

elif data_type == "AA" :
	print "Matrix data type is: AA"
	if ambig_chars == None :
		print "No set of ambiguous characters was passes with -a option, so assuming 'X?'."
		ambig_chars = ["X","?"]
	else :
		ambig_chars = list(ambig_chars)
		print "The following characters are treated as ambiguous or missing: %s " %(ambig_chars)

elif data_type == "MULTI" :
	print "Matrix data type is: MULTI"
	if ambig_chars == None :
		print "No set of ambiguous characters was passes with -a option, so assuming '?'."
		ambig_chars = ["?"]
	else :
		ambig_chars=list(ambig_chars)
		print "The following characters are treated as ambiguous or missing: %s " %(ambig_chars)

else :
	print "Matrix character type was not set correctly %s is not one of the three options (DNA, AA, MULTI)." %(data_type)
	quit()

char_len=0	

phylip=re.compile('\d*\s\d*')	#String to matching header in phylip file.

for Line in IN :		#Read in_file into a dictionary with keys as taxa and values as sequences (all sequences are same length)
	Line = Line.strip('\n')
	if phylip.match(Line) : # If first line matches phylip header string defined above.
		print "Input file detected as phylip format."
	else :
		Line_bits=re.split("\s+",Line)
	
		if Line_bits[0] == '' :		#Skip blank lines
			continue
		if char_len == 0 :				#Set Length from first sequence
			char_len = len(Line_bits[1])
			print "All input character sets should be %s characters long." %(char_len)
		else :
			if char_len != len(Line_bits[1]) :
				print "Error, %s has different number of characters (%d is not %d)" %(Line_bits[0], char_len, len(Line_bits[1]) )	#if character lengths are not the same exit.
				break

		char_dict[Line_bits[0]]=Line_bits[1]
		trans_char_dict[Line_bits[0]]=''			#Add the taxon names to all the dictionaries for output.
		trans_bin_dict[Line_bits[0]]=''
		trans_bin_csv_dict[Line_bits[0]]=''
		trans_bin_reduced_dict[Line_bits[0]]=''
		trans_bin_reduced_csv_dict[Line_bits[0]]=''
		trans_polarized_dict[Line_bits[0]]=''
		trans_polarized_csv_dict[Line_bits[0]]=''
		trans_polarized_reduced_dict[Line_bits[0]]=''
		trans_polarized_reduced_csv_dict[Line_bits[0]]=''
	
		Last_taxon=Line_bits[0] #keep resetting until the true last taxon.

IN.close 

if outgroup == None:  #If outgroup was not specified on command line, set outgroup to Last_taxon.
	outgroup=Last_taxon
	
print "The outgroup taxon is %s." %(outgroup)

#Go through the character dictionary character by character, and make the output dictionaries from this.
char_count=1
while char_count <= char_len :
	#Get character of outgroup
	out_char=char_dict[outgroup][(char_count-1):char_count].upper()  #Get the character and convert to uppercase
	
	#Go through rest of dictionary and get the characters present in other taxa at this position.
	chars=[]
	for key in char_dict :
		tax_char=char_dict[key][(char_count-1):char_count].upper()
		if (tax_char not in ambig_chars) and (tax_char not in chars) :
			chars.append(tax_char)
			
	#print "Character %d of outgroup is character %s, other characters are: %s" %(char_count,out_char,chars)
	
	#For each character state that is present for this position comput the output as needed.
	for character in chars :
		for key in char_dict :		# For each taxon...
			tax_char=char_dict[key][(char_count-1):char_count].upper()	#Get the character
			if tax_char == character : # if the character in this taxon is the same as the current character state
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "1" #add a 1
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",1"
				
				if len(chars) > 1 :	#if position is not invariant, add to reduced dataset.
					trans_bin_reduced_dict[key]=str(trans_bin_reduced_dict[key]) + "1" #add a 1
					trans_bin_reduced_csv_dict[key]=str(trans_bin_reduced_csv_dict[key]) + ",1"
					
			elif tax_char in ambig_chars : #if character in this taxon is ambiguous
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "?" #add a ?
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",?"
				if len(chars) > 1 :	#if position is not invariant, add to reduced dataset.
					trans_bin_reduced_dict[key]=str(trans_bin_reduced_dict[key]) + "?" #add a ?
					trans_bin_reduced_csv_dict[key]=str(trans_bin_reduced_csv_dict[key]) + ",?"
			else  : #otherwise add a 0
				trans_bin_dict[key]=str(trans_bin_dict[key]) + "0" #add a 0
				trans_bin_csv_dict[key]=str(trans_bin_csv_dict[key]) + ",0"
				if len(chars) > 1 :	#if position is not invariant, add to reduced dataset.
					trans_bin_reduced_dict[key]=str(trans_bin_reduced_dict[key]) + "0" #add a 1
					trans_bin_reduced_csv_dict[key]=str(trans_bin_reduced_csv_dict[key]) + ",0"
			
			if (out_char in ambig_chars) or (tax_char in ambig_chars) : #if outgroup or current taxon character is in the ambiguous list
				trans_polarized_dict[key]=str(trans_polarized_dict[key]) + "?"
				trans_polarized_csv_dict[key]=str(trans_polarized_csv_dict[key]) + ",?"
				if (len(chars) > 1) and (out_char not in ambig_chars) :	#if position is not invariant, and outgroup isn't ambig add to reduced dataset.
					trans_polarized_reduced_dict[key]=str(trans_polarized_reduced_dict[key]) + "?" #add a ?
					trans_polarized_reduced_csv_dict[key]=str(trans_polarized_reduced_csv_dict[key]) + ",?"	
			elif tax_char == out_char :		#polarize for outgroup taxon: if = outgroup, then 0, else 1
				trans_polarized_dict[key]=str(trans_polarized_dict[key]) + "0"
				trans_polarized_csv_dict[key]=str(trans_polarized_csv_dict[key]) + ",0"
				if (len(chars) > 1)  :	#if position is not invariant, add to reduced dataset.
					trans_polarized_reduced_dict[key]=str(trans_polarized_reduced_dict[key]) + "0" #add a 0
					trans_polarized_reduced_csv_dict[key]=str(trans_polarized_reduced_csv_dict[key]) + ",0"	
			else :
				trans_polarized_dict[key]=str(trans_polarized_dict[key]) + "1"
				trans_polarized_csv_dict[key]=str(trans_polarized_csv_dict[key]) + ",1"
				
				if (len(chars) > 1):	#if position is not invariant, add to reduced dataset.
					trans_polarized_reduced_dict[key]=str(trans_polarized_reduced_dict[key]) + "1" #add a 1
					trans_polarized_reduced_csv_dict[key]=str(trans_polarized_reduced_csv_dict[key]) + ",1"	
					

# 	for key in char_dict :			#Add extra column with Outgroup characters or ? to the UIB dataset
# 		tax_char=char_dict[key][(char_count-1):char_count].upper()
# 		if tax_char == out_char :
# 			trans_uib_dict[key]=str(trans_uib_dict[key]) + out_char
# 		else :
# 			trans_uib_dict[key]=str(trans_uib_dict[key]) + "?"
			
	char_count+=1
			

#to print everything in the same order, go back through the original file to get taxa to pull out of dictionary.
try:
	IN = open (in_file, 'r')
except OSError as exception:
	pass

#Open all the output files.
try:
	OUT_BIN = open (out_bin_all, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_CSV = open (out_bin_all_csv, 'w')
except OSError as exception:
	pass
	
try:
	OUT_BIN_REDUCED = open (out_bin_reduced, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_REDUCED_CSV = open (out_bin_reduced_csv, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_POLARIZED = open (out_bin_polarized, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_POLARIZED_CSV = open (out_bin_polarized_csv, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_POLARIZED_REDUCED = open (out_bin_polarized_reduced, 'w')
except OSError as exception:
	pass

try:
	OUT_BIN_POLARIZED_REDUCED_CSV = open (out_bin_polarized_reduced_csv, 'w')
except OSError as exception:
	pass


#Add phylip headers to phylip output files, csv don't need header
num_tax= len(trans_char_dict)

OUT_BIN.write ("%d\t%d\n" %(num_tax,len(trans_bin_dict[Last_taxon])))
OUT_BIN_REDUCED.write ("%d\t%d\n" %(num_tax,len(trans_bin_reduced_dict[Last_taxon])))
OUT_BIN_POLARIZED.write ("%d\t%d\n" %(num_tax,len(trans_polarized_dict[Last_taxon])))
OUT_BIN_POLARIZED_REDUCED.write ("%d\t%d\n" %(num_tax,len(trans_polarized_reduced_dict[Last_taxon])))


for Line in IN :
	if phylip.match(Line) : # If first line matches phylip header string defined above.
		continue
	else:
		Line = Line.strip('\n')
		Line_bits=re.split("\s+",Line)
		if Line_bits[0] == '' :
			continue
		else : 
			OUT_BIN.write ("%s\t%s\n" %(Line_bits[0],trans_bin_dict[Line_bits[0]]))
			OUT_BIN_REDUCED.write ("%s\t%s\n" %(Line_bits[0],trans_bin_reduced_dict[Line_bits[0]]))
			OUT_BIN_CSV.write ("%s\t%s\n" %(Line_bits[0],trans_bin_csv_dict[Line_bits[0]]))
			OUT_BIN_REDUCED_CSV.write ("%s\t%s\n" %(Line_bits[0],trans_bin_reduced_csv_dict[Line_bits[0]]))
			OUT_BIN_POLARIZED.write ("%s\t%s\n" %(Line_bits[0],trans_polarized_dict[Line_bits[0]]))
			OUT_BIN_POLARIZED_CSV.write ("%s\t%s\n" %(Line_bits[0],trans_polarized_csv_dict[Line_bits[0]]))
			OUT_BIN_POLARIZED_REDUCED.write ("%s\t%s\n" %(Line_bits[0],trans_polarized_reduced_dict[Line_bits[0]]))
			OUT_BIN_POLARIZED_REDUCED_CSV.write ("%s\t%s\n" %(Line_bits[0],trans_polarized_reduced_csv_dict[Line_bits[0]]))


