#!/bin/bash
WEBSITE_DIR=$CPL_PATH/website
WEBSITE_DOC_DIR=$WEBSITE_DIR/api-docs/
DOCS_DIR=$CPL_PATH/doc

cd $DOCS_DIR
make html
cp -R build/html/* $WEBSITE_DOC_DIR
make clean
cd $WEBSITE_DOC_DIR
mv python_api.html python_api.shtml
mv fortran_api.html fortran_api.shtml
mv cpp_api.html cpp_api.shtml
