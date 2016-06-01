#! /usr/bin/env python
from os import system 
import check

fobj    = open('TOPOL.in','r')

npx_md  = int(fobj.readline()[0:4])
npy_md  = int(fobj.readline()[0:4])
npz_md  = int(fobj.readline()[0:4])
npx_cfd = int(fobj.readline()[0:4])
npy_cfd = int(fobj.readline()[0:4])
npz_cfd = int(fobj.readline()[0:4])

fobj.close()	

nproc = npx_md*npy_md*npz_md + npx_cfd*npy_cfd*npz_cfd
cmd = 'mpiexec -n ' + str(nproc) + ' ./a.out' 
print(cmd)
system(cmd)

cmd = 'cat fort.1* > info_realms     2> /dev/null && rm fort.1*'
system(cmd)
cmd = 'cat fort.2* > info_MD_recv    2> /dev/null && rm fort.2*'
#cmd = 'cat fort.2* > info_olap_MD    2> /dev/null && rm fort.2*'
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

check.vals('info_scatter_md')
check.vals('info_gather_cfd')
check.vals('info_CFD_send')
check.vals('info_MD_send')
check.vals('info_MD_recv')
check.vals('info_CFD_recv')
