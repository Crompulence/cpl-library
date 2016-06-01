cd ./openfoam
    python clean.py
    blockMesh
    decomposePar
cd -  
rm cpl/coupler_header cpl/map_*
rm lammps/cplchunk
