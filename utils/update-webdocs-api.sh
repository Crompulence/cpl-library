#!/bin/bash
WEBSITE_DIR=$CPL_PATH/website
WEBSITE_DOC_DIR=$WEBSITE_DIR/api-docs/
cd $WEBSITE_DOC_DIR

#Sphinx is a Python documentation tool so we use
# 1) a Fortran convertor called "sphinx-fortran" for the Fortran code
# 2) Doxygen for the C++ part and breathe to convert to sphinx format
pip install -U sphinx sphinx-fortran breathe
REQUIRED_PKG="doxygen"
PKG_OK=$(dpkg-query -W --showformat='${Status}\n' $REQUIRED_PKG|grep "install ok installed")
echo Checking for $REQUIRED_PKG: $PKG_OK
if [ "" = "$PKG_OK" ]; then
  echo "No $REQUIRED_PKG. Setting up $REQUIRED_PKG."
  sudo apt-get --yes install $REQUIRED_PKG
fi

#Generate documents using Sphinx, etc
make clean
make html

#Copy all html folders to server side format
cp -R build/html/* $WEBSITE_DOC_DIR
mv python_api.html python_api.shtml
mv fortran_api.html fortran_api.shtml
mv cpp_api.html cpp_api.shtml
mv index.html index.shtml
mv search.html search.shtml
mv genindex.html genindex.shtml
mv f-modindex.html f-modindex.shtml
for file in *.shtml; 
	do echo $file; 
	sed -i 's/\.html/\.shtml/' $file;
done

#Update module documents to server side
cd _modules
mv cplpy.html cplpy.shtml
mv index.html index.shtml
for file in *.shtml; 
	do echo $file; 
	sed -i 's/\.html/\.shtml/' $file;
done

