cd ./test_forces_case
    python clean.py
    blockMesh
    decomposePar
cd -  
rm cpl/coupler_header cpl/map_*
