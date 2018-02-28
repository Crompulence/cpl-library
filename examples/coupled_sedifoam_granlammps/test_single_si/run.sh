#!/bin/bash

#Set OpenFOAM and LAMMPS folder and binaries. Default values provided below, unless overwriten by argument.

OPEN_FOAM_CASE=openfoam/
LAMMPS_CASE=lammps/
OPEN_FOAM_BIN=../bin/CPLSediFOAM
LAMMPS_BIN=../bin/lmp_cpl

iFLag=0
while getopts ":o:O:l:L:i:" opt; 
do
	case $opt in
	    o)
	      	echo "-o was triggered!"
	      	OPEN_FOAM_CASE=$OPTARG
	      	;;
	  	O)
		  	echo "-O was triggered!"
		  	OPEN_FOAM_BIN=$OPTARG
		  	;;
	  	l)
	      	echo "-l was triggered!"
	      	LAMMPS_CASE=$OPTARG
	      	;;
	  	L)
		  	echo "-L was triggered!"
		  	LAMMPS_BIN=$OPTARG
		  	;;
	  	i)
			echo "-i was triggered!"
			LAMMPS_INPUT="${LAMMPS_CASE}$OPTARG"
			iFlag=1
		  	;;
	    \?)
	      	echo "Invalid option: -$OPTARG"
	      	exit 1
	      	;;
	  	:)
		  	echo "Option $OPTARG requires an argument."
		  	exit 1
		  	;;
	esac
done
shift $((OPTIND - 1))

if [ $iFlag -eq 0 ];
then
	echo "LAMMPS input file must be provided."
	exit 1
fi

#Clean CFD
cd ${OPEN_FOAM_CASE}
rm -f ../log.openfoam
python clean.py -f
blockMesh
decomposePar
cd ../

#Clean MD
cd ${LAMMPS_CASE}
rm -f ../log.lammps vmd_out.dcd thermo_output.txt particle_dump/*
cd ../

#Run simulation
cplexec -c 1 "${OPEN_FOAM_BIN} -case ${OPEN_FOAM_CASE} -parallel > log.openfoam" -m 1 "${LAMMPS_BIN} < ${LAMMPS_INPUT}"
