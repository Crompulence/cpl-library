
#######################
#
#  Prerequists
#
#######################
ENVNAME=cplrun
INSTALL_DIR=$HOME/codes
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

#Apt get some things
sudo apt-get update
sudo apt-get install -y gfortran
sudo apt-get install -y git-core

#Get git so we can get all the codes
git clone https://edwardsmith999@bitbucket.org/granlammps/gitlammps.git $INSTALL_DIR/granlammps
git clone https://github.com/Crompulence/cpl-library.git $INSTALL_DIR/cpl-library
git clone https://github.com/Crompulence/cpl_conda_builds.git $INSTALL_DIR/cpl_conda_builds
git clone https://github.com/Crompulence/CPL_APP_OPENFOAM-3.0.1.git $INSTALL_DIR/CPL_APP_OPENFOAM-3.0.1
git clone https://github.com/Crompulence/CPL_APP_LAMMPS-DEV.git $INSTALL_DIR/CPL_APP_LAMMPS-DEV

#Also get utilities and post processing
git clone https://github.com/edwardsmith999/pyDataView.git $INSTALL_DIR/pyDataView
git clone https://github.com/edwardsmith999/SimWrapPy.git $INSTALL_DIR/SimWrapPy
git clone https://github.com/edwardsmith999/drag-utils.git $INSTALL_DIR/dragutils

#Get miniconda
unset PYTHONPATH
unset LD_LIBRARY_PATH
wget https://repo.continuum.io/miniconda/Miniconda2-latest-Linux-x86_64.sh -O miniconda.sh
bash miniconda.sh -b -p $INSTALL_DIR/miniconda
export PATH="$INSTALL_DIR/miniconda/bin:$PATH"
rm -f miniconda.sh

#Configure conda, get numpy and scipy
conda config --set always_yes yes --set changeps1 no
conda update -q conda
conda info -a
conda create -q -n $ENVNAME python=2.7 numpy scipy
source activate $ENVNAME
conda config --add channels edu159

#We need to explicitly get gcc/gfortran as miniconda is not available for 5+ yet
conda install -y gxx_linux-64
conda install -y gfortran_linux-64

#Here we install MPI version mpich and mpi4py
conda remove mpich           #Note sure why but mpich 1.4 appeared here
conda install -c edu159 -y mpich
conda install -c edu159 -y mpi4py

#Need a few other utils for Python
conda install -y pytest
conda install -y matplotlib
conda install -y wxpython 


#Go to virtual enviroment directory
ENVPREFIX=$INSTALL_DIR/miniconda/envs/$ENVNAME
cd $ENVPREFIX

#######################
#
#  CPL library
#
#######################

#Build CPL library using conda version of code
cd $INSTALL_DIR/cpl-library
make PLATFORM=gcc clean
make PLATFORM=gcc
ldd lib/*
cp lib/* $ENVPREFIX/lib
source SOURCEME.sh
cd ../

#######################
#
#      OpenFOAM
#
#######################

#Get conda build so we can copy required stuff
RECIPE_DIR=$ENVPREFIX/cpl_conda_builds/cplapp-openfoam3.0.1/
cp -Rf ${RECIPE_DIR}/linux64Gcc $ENVPREFIX/opt/$FOAM_DIR_NAME/wmake/rules

#Install OpenFOAM from conda
conda install -c edu159 -y openfoam

#Setup aliases
FOAM_VERSION=3.0.1
FOAM_DIR_NAME=OpenFOAM-$FOAM_VERSION

#Get all of CPL APP
cd $INSTALL_DIR/CPL_APP_OPENFOAM-3.0.1
echo "$ENVPREFIX/opt" > CODE_INST_DIR
source SOURCEME.sh

#We need to make Pstream and patch the OpenFOAM version
make clean
make pstream
mv $FOAM_LIBBIN/$FOAM_MPI/libPstream.so $FOAM_LIBBIN/$FOAM_MPI/libPstream.so.orig
cp lib/libPstream.so $FOAM_LIBBIN/$FOAM_MPI/
make
cp lib/* $ENVPREFIX/lib
cp bin/* $ENVPREFIX/bin

make test-hydrostatic
make test-fcc_dummy

#######################
#
#  GranLAMMPS
#
#######################

#Build version of LAMMPS
cd $INSTALL_DIR/granlammps
git checkout common
pwd > ../CPL_APP_LAMMPS-DEV/CODE_INST_DIR
cd ../CPL_APP_LAMMPS-DEV
echo granular >> config/lammps_packages.in
cd config
sh ./enable-packages.sh make
cd ../
make patch-lammps
make -j 4
source SOURCEME.sh

make test-single
make test-simwrap





