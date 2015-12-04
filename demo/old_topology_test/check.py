#! /usr/bin/env python2.7
import os

def vals(filename):

	# Open file object
	fobj = open(filename,'r')

	# Initialise flag that stays true only if all read values
	# are "correct"
	all_correct = True

	size = os.path.getsize(filename)
	if (size == 0):
		all_correct = False
#		string = 'File is empty '+filename 
#		print(string) 

	# Read cells and values
	for line in iter(fobj.readline,""):

		cells = line.split("(")[1].split(")")[0].split(",")
		value = float(line.split("=")[1].strip())	
		ixyz  = float(cells[0])
		icell = float(cells[1])
		jcell = float(cells[2])
		kcell = float(cells[3])

		# Calculate expected "value" and error
		expec_val = 0.1*ixyz + 1.0*icell + 1000*jcell + 1000000*kcell  
		error = abs(value - expec_val)

		# If error is nonzero print message and set all_correct
		# flag to false
		if (error != 0.0):
			all_correct = False
			string = 'Error in recvarray('+str(int(ixyz))+','
			string += str(int(icell))+','+str(int(jcell))+','
			string += str(int(kcell))+'). Expected: '
			string += str(expec_val)+'      Read: '
			string += str(value)+'.'
			print(string)

	fobj.close()

	# Success message
	if (all_correct == True):
		#string = 'All printed results are as expected in '+filename 
		#print(string)
		return 0
	else:
		#string = 'Printed results not all as expected in '+filename 
		#print(string)
		return 1
	
