
cd openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

cplexec -c 2 "CPLIcoFoam -case ./openfoam -parallel" -m 2 "./lmp_cpl -i lammps.in" 

