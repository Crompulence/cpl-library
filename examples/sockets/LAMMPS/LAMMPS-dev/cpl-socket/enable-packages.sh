if [ ! -d "$LAMMPS_PATH" ]; then
    echo "LAMMPS PATH not found, have you sourced SOURCEME file?"
    exit
fi

cd $LAMMPS_PATH/src
while read p; do
  make "no-$p"
  make "yes-$p"
done <$LAMMPS_PATH/cpl-socket/lammps_packages.in
make mode=lib cpl
