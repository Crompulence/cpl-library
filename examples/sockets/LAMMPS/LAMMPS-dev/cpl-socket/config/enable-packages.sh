cd $LAMMPS_PATH/src
while read p; do
  make "no-$p"
  make "yes-$p"
done <$LAMMPS_PATH/../cpl-socket/lammps_packages.in
make mode=lib cpl
