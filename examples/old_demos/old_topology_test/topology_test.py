#! /usr/bin/env python2.7
from testmodule import *
import os

# Clean and/or setup logs directory
if (os.path.exists('./logs')):
	os.system('rm -r ./logs')
os.mkdir('./logs')

####	P A R A M E T E R    S T U D I E S

# Number of cfd procs in y
print('========================================')
print(' -- Study: numbers of cfd procs in y...')
print('========================================')
set_defaults()
params = [1,2,3,4,6,8]
for PARAM in params:
	job = RunClass(npx_md,      npy_md,      npz_md,  
				   npx_cfd,     PARAM,       npz_cfd,
				   ncx,         ncy,         ncz,    
				   xL_cfd,      yL_cfd,      zL_cfd, 
				   icmin_olap,  icmax_olap,    
				   jcmin_olap,  jcmax_olap,    
				   kcmin_olap,  kcmax_olap            )
	print(job)
	job.execute()
	job.concatenate()
	job.checkvalues()

# Odd numbers of cfd procs
print('========================================')
print(' -- Study: odd numbers of cfd procs')
print('========================================')
set_defaults()
params = [1,3,5,7]
for PARAM in params:
	job = RunClass(PARAM*4,     npy_md,      npz_md,  
				   PARAM,       npy_cfd,     npz_cfd,
				   PARAM*64,    ncy,         ncz,    
				   xL_cfd,      yL_cfd,      zL_cfd, 
				   1, PARAM*64,    
				   1, jcmax_olap,    
				   1, kcmax_olap            )
	print(job)
	job.execute()
	job.concatenate()
	job.checkvalues()

# jcmax_olap 
print('========================================')
print(' -- Study: jcmax_olap...')
print('========================================')
set_defaults()
params = [1,2,3,5,10,20,100,150]
for PARAM in params:
	job = RunClass(npx_md,  npy_md,  npz_md,  
				   npx_cfd, npy_cfd, npz_cfd,
				   ncx,     ncy,     ncz,    
				   xL_cfd,  yL_cfd,  zL_cfd, 
				   icmin_olap, icmax_olap,   
				   jcmin_olap, PARAM,   
				   kcmin_olap, kcmax_olap )
	print(job)
	job.execute()
	job.concatenate()
	job.checkvalues()

# overlap region smaller than domain in x?
print('========================================')
print(' -- Study: small olap region')
print('========================================')
set_defaults()
job = RunClass(npx_md,  npy_md,  npz_md,  
			   npx_cfd, npy_cfd, npz_cfd,
			   ncx,     ncy,     ncz,    
			   xL_cfd,  yL_cfd,  zL_cfd, 
			   (2*16+1), (6*16) ,   
			   jcmin_olap, jcmax_olap,   
			   kcmin_olap, kcmax_olap )
print(job)
job.execute()
job.concatenate()
job.checkvalues()
