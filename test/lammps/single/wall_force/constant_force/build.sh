MYPWD=${PWD}
rm -f ./lmp_cpl
cd ../../../../GranLAMMPS/src/
make package-update
make cpl
cp ./lmp_cpl ${MYPWD}/
cd ${MYPWD}
