mypwd=$(pwd)
cd ../../
if [ $1 = "1" ]; then
	make clean
fi
make BUILD=debug 
cd $mypwd
if [ $1 = "1" ]; then
	make clean
fi
make
