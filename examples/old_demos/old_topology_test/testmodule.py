#! /usr/bin/env python2.7
"""
	TOPOLOGY TEST MODULE

	Contents:
		
		prevjobID (integer)   -   Parameter incremented every time an object of
		                          type RunClass is created, then saved to
		                          that object's self.jobID

		( various input parameters )
		
		RunClass  (class)     -   Everything required for a job object.

"""
infile = 'TOPOL.in'
prevjobID = 0

npx_md = 8
npy_md = 4
npz_md = 1
npx_cfd = 2
npy_cfd = 2
npz_cfd = 1
ncx = 128
ncy = 128
ncz = 8
xL_cfd = 1280.0
yL_cfd = 1280.0 
zL_cfd = 80.0
icmin_olap = 1
icmax_olap = 128
jcmin_olap = 1
jcmax_olap = 23
kcmin_olap = 1
kcmax_olap = 8

class RunClass:

	# Constructor
	def __init__(self, npx_md,  npy_md,  npz_md, 
	                   npx_cfd, npy_cfd, npz_cfd,
	                   ncx,     ncy,     ncz,
	                   xL_cfd,  yL_cfd,  zL_cfd,
	                   icmin_olap, icmax_olap, 
	                   jcmin_olap, jcmax_olap,
	                   kcmin_olap, kcmax_olap):
		from re import sub 
		from os import system

		global prevjobID

		# Remove any previous "info" files
		system('rm info_*')

		# Create new ID for every job initialised
		self.jobID = prevjobID + 1
		prevjobID  = self.jobID

		# Read and store input file in "lines"
		with open(infile,'r') as source:
			lines = source.readlines()

		# Clean input file and rewrite with new procs
		with open(infile,'w') as output:
			# Store inputs in an array for loop below
			inputs = ([ npx_md,  npy_md,  npz_md,  
						npx_cfd, npy_cfd, npz_cfd,
						ncx,     ncy,     ncz,    
						xL_cfd,  yL_cfd,  zL_cfd, 
						icmin_olap, icmax_olap,   
						jcmin_olap, jcmax_olap,   
						kcmin_olap, kcmax_olap    ])

			for line in range(len(lines)):
				subcount = 1
				output.write( 
				              sub ( '\S{1,}'   , str(inputs[line]), 
				                    lines[line], subcount           )
				            )

		# Store nprocs
		self.nproc = ( npx_md  * npy_md  * npz_md 
		             + npx_cfd * npy_cfd * npz_cfd )

	# Print statement
	def __str__(self):
		
		string =  '-----------------------------------------\n'
		string += 'jobID: ' + str(self.jobID) + '\n'
		string += '-----------------------------------------'

		# Add input file contents
		#with open(infile, 'r' ) as source:
		#	lines = source.readlines()
		#for line in lines:
		#	string += line

		#string += '-----------------------------------------'

		return( string )

	
	def execute(self):
		import os	


		if (os.path.exists('./logs')):

			logfile = './logs/' + str(self.jobID) + '_log'
			errfile = './logs/' + str(self.jobID) + '_errlog'

			cmd = 'cat input > ' + logfile
			os.system(cmd)

			cmd = ('mpiexec -n ' + str(self.nproc) + ' ./a.out >> ' + logfile 
				   + ' 2>> ' + errfile)

			print(cmd)	
			os.system(cmd)

		else:

			raise


	def concatenate(self):
		from os import system 

		cmd = 'cat fort.1* > info_realms     2> /dev/null && rm fort.1*'
		system(cmd)
		cmd = 'cat fort.2* > info_MD_recv    2> /dev/null && rm fort.2*'
		system(cmd)
		cmd = 'cat fort.3* > info_graph      2> /dev/null && rm fort.3*'
		system(cmd)
		cmd = 'cat fort.4* > info_MD_send    2> /dev/null && rm fort.4*'
		system(cmd)
		cmd = 'cat fort.5* > info_CFD_recv   2> /dev/null && rm fort.5*'
		system(cmd)
		cmd = 'cat fort.6* > info_maps       2> /dev/null && rm fort.6*'
		system(cmd)
		cmd = 'cat fort.7* > info_scatter_md 2> /dev/null && rm fort.7*'
		system(cmd)
		cmd = 'cat fort.8* > info_gather_cfd 2> /dev/null && rm fort.8*'
		system(cmd)
		cmd = 'cat fort.9* > info_CFD_send   2> /dev/null && rm fort.9*'
		system(cmd)

	def checkvalues(self):
		import check

		all_success = True

		success = check.vals('info_scatter_md')
		if success is not 0:
			print('Error in gather values')
			all_success = False
		success = check.vals('info_gather_cfd')
		if success is not 0:
			print('Error in scatter values')
			all_success = False
		success = check.vals('info_CFD_send')
		if success is not 0:
			print('Error in CFD send values')
			all_success = False
		success = check.vals('info_MD_send')
		if success is not 0:
			print('Error in MD send values')
			all_success = False
		success = check.vals('info_MD_recv')
		if success is not 0:
			print('Error in MD recv values')
			all_success = False
		success = check.vals('info_CFD_recv')
		if success is not 0:
			print('Error in CFD recv values')
			all_success = False

	
		if all_success:
			print('jobID ' + str(self.jobID) + ' value checks successful!')



def set_defaults():

	global npx_md
	global npy_md
	global npz_md
	global npx_cfd
	global npy_cfd
	global npz_cfd

	global ncx
	global ncy
	global ncz 

	global xL_cfd
	global yL_cfd
	global zL_cfd
	
	global icmin_olap
	global icmax_olap
	global jcmin_olap
	global jcmax_olap
	global kcmin_olap
	global kcmax_olap

	npx_md = 8
	npy_md = 4
	npz_md = 1

	npx_cfd = 2
	npy_cfd = 2
	npz_cfd = 1

	ncx = 128
	ncy = 128
	ncz = 8

	xL_cfd = 1280.0
	yL_cfd = 1280.0 
	zL_cfd = 80.0
	
	icmin_olap = 1
	icmax_olap = 128
	jcmin_olap = 1
	jcmax_olap = 23
	kcmin_olap = 1
	kcmax_olap = 8
