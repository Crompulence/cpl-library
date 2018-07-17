
#######################
#
#  Prerequists
#
#######################
INSTALL_DIR=$HOME/codes
mkdir -p $INSTALL_DIR
cd $INSTALL_DIR

#Apt get some things
sudo add-apt-repository universe
sudo apt-get update
sudo apt-get install -y gfortran
sudo apt-get install -y git-core
sudo apt-get purge -y --auto-remove openmpi-bin
sudo apt-get install -y mpich

#Get all the codes on github
git clone https://edwardsmith999@bitbucket.org/granlammps/gitlammps.git $INSTALL_DIR/granlammps

#OpenFOAM github seems to be missing third party stuff
#git clone https://github.com/OpenFOAM/ThirdParty-3.0.x.git ThirdParty-3.0.x
#git clone https://github.com/OpenFOAM/OpenFOAM-3.0.x.git OpenFOAM-3.0.x
wget http://downloads.sourceforge.net/foam/OpenFOAM-3.0.1.tgz
tar -xvf OpenFOAM-3.0.1.tgz
wget http://downloads.sourceforge.net/foam/ThirdParty-3.0.1.tgz
tar -xvf ThirdParty-3.0.1.tgz

#Crompulence cpl library and two apps
git clone https://github.com/Crompulence/cpl-library.git $INSTALL_DIR/cpl-library
git clone https://github.com/Crompulence/CPL_APP_OPENFOAM-3.0.1.git $INSTALL_DIR/CPL_APP_OPENFOAM-3.0.1
git clone https://github.com/Crompulence/CPL_APP_LAMMPS-DEV.git $INSTALL_DIR/CPL_APP_LAMMPS-DEV

#Also get utilities and post processing
git clone https://github.com/edwardsmith999/pyDataView.git $INSTALL_DIR/pyDataView
git clone https://github.com/edwardsmith999/SimWrapPy.git $INSTALL_DIR/SimWrapPy
git clone https://github.com/edwardsmith999/drag-utils.git $INSTALL_DIR/dragutils

#Need a few other utils for Python
sudo apt-get install -y ipython
sudo apt-get install -y python-pip
sudo pip install numpy scipy matplotlib pytest
sudo apt-get install -y python-tk
sudo apt-get install -y python-wxgtk3.0
wget https://bitbucket.org/mpi4py/mpi4py/downloads/mpi4py-2.0.0.tar.gz
tar -xvf mpi4py-2.0.0.tar.gz
cd mpi4py-2.0.0
sudo python setup.py install
cd ../

#######################
#
#  CPL library
#
#######################

#Build CPL library using conda version of code
cd $INSTALL_DIR/cpl-library
make PLATFORM=gcc clean
make PLATFORM=gcc
source SOURCEME.sh

make test-examples
make test-valgrind
make test-gtests
make test_Dragmodels

cd ../

#######################
#
#  GranLAMMPS
#
#######################

#Build version of LAMMPS
cd $INSTALL_DIR/granlammps
git checkout common
pwd > ../CPL_APP_LAMMPS-DEV/CODE_INST_DIR

#Build CPL APP
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

cd ../

#######################
#
#      OpenFOAM
#
#######################

#Prerequists
sudo apt-get install -y flex-old
sudo apt-get install -y libz-dev
sudo apt-get install -y libboost-system-dev libboost-thread-dev bison libreadline-dev libncurses-dev libxt-dev

#Set some aliasa
FOAM_VERSION=3.0.1
FOAM_SRC_DIR=$INSTALL_DIR/OpenFOAM-$FOAM_VERSION
APP_DIR=$INSTALL_DIR/CPL_APP_OPENFOAM-$FOAM_VERSION

#We copy this pref file to build OpenFOAM with system MPICH instead of OpenMPI
cp $APP_DIR/config/prefs_system_mpich.sh $FOAM_SRC_DIR/etc/pref.sh

#Build from CPL APP file
cd $INSTALL_DIR/CPL_APP_OPENFOAM-3.0.1
echo $FOAM_SRC_DIR > $APP_DIR/CODE_INST_DIR
source SOURCEME.sh  # Also calls source $FOAM_SRC_DIR/etc/bashrc

# Build on multiple processes
export WM_NCOMPPROCS=4

#Build Third Party code
cd $INSTALL_DIR/ThirdParty-$FOAM_VERSION
./Allwmake

# -- COMPILE -- 
cd $FOAM_SRC_DIR
./Allwmake -j

#We need to make Pstream and patch the OpenFOAM version
cd $INSTALL_DIR/CPL_APP_OPENFOAM-3.0.1
make clean
make pstream
mv $FOAM_LIBBIN/$FOAM_MPI/libPstream.so $FOAM_LIBBIN/$FOAM_MPI/libPstream.so.orig
cp lib/libPstream.so $FOAM_LIBBIN/$FOAM_MPI/
make
cp lib/* $ENVPREFIX/lib
cp bin/* $ENVPREFIX/bin

make test-hydrostatic
make test-fcc_dummy

cd ../


