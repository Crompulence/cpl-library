# Modify scripts from base directory in ~/Document and then rsync to
# respective folders in CPL-library, LAMMPS-APP and OpenFOAM-APP
rsync -arv --exclude='*pyc' --exclude='__pycache__/' /home/asufian/Documents/Documentation/cfd-dem-theory-validation-examples-documentation/examples/python_scripts/ python_scripts