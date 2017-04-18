# Clean include and lib dirs
rm -rf $CPL_THIRD_PARTY/include $CPL_THIRD_PARTY/lib > /dev/null

# Clean json-fortran
json_builddir=$CPL_THIRD_PARTY/json-fortran/build
if [ -d  $json_builddir ]; then
    echo "Cleaning json-fortran project..."
	rm -rf $json_builddir
fi
