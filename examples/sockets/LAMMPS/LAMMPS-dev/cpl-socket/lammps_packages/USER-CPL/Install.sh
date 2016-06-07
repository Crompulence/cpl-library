# Install/unInstall package files in LAMMPS
# mode = 0/1/2 for uninstall/install/update

# this is default Install.sh for all packages
# if package has an auxiliary library or a file with a dependency,
# then package dir has its own customized Install.sh

mode=$1

# arg1 = file, arg2 = file it depends on

action () {
  if (test $mode = 0) then
    rm -f ../$1
  elif (! cmp -s $1 ../$1) then
    if (test -z "$2" || test -e ../$2) then
      cp $1 ..
      if (test $mode = 2) then
        echo "  updating src/$1"
      fi
    fi
  elif (test -n "$2") then
    if (test ! -e ../$2) then
      rm -f ../$1
    fi
  fi
}

# all package files with no dependencies
for file in *.cpp *.h; do
  action $file
done

# Add CPL includes and lib
if (test $1 = 1) then

  if (test -e ../Makefile.package) then
    sed -i -e 's/[^ \t]*atc[^ \t]* //' ../Makefile.package
    sed -i -e 's|^PKG_INC =[ \t]*|&-I$(CPL_PATH)/include |' ../Makefile.package
    sed -i -e 's|^PKG_PATH =[ \t]*|&-L$(CPL_LIBRARY_PATH) -Wl,-rpath=$(CPL_LIBRARY_PATH) |' ../Makefile.package
    sed -i -e 's|^PKG_LIB =[ \t]*|&-lcpl -lmpifort -lgfortran |' ../Makefile.package
#    sed -i -e 's|^PKG_SYSINC =[ \t]*|&$(user-cpl_SYSINC) |' ../Makefile.package
#    sed -i -e 's|^PKG_SYSLIB =[ \t]*|&$(user-cpl_SYSLIB) |' ../Makefile.package
#    sed -i -e 's|^PKG_SYSPATH =[ \t]*|&$(user-cpl_SYSPATH) |' ../Makefile.package
  fi

# if (test -e ../Makefile.package.settings) then
#   sed -i -e '/^include.*cpl.*$/d' ../Makefile.package.settings
#   # multiline form needed for BSD sed on Macs
#    sed -i -e '4 i \
#include ..\/..\/lib\/cpl\/Makefile.lammps
#' ../Makefile.package.settings
#  fi

elif (test $1 = 0) then

  if (test -e ../Makefile.package) then
    sed -i -e 's|[^ \t]*(CPL_PATH)/include*[^ \t] ||g' ../Makefile.package
    sed -i -e 's|[^ \t]*CPL_LIBRARY_PATH*[^ \t] ||g' ../Makefile.package
    sed -i -e 's|[^ \t]*cpl*[^ \t] ||' ../Makefile.package
    sed -i -e 's|[^ \t]*mpifort*[^ \t] ||' ../Makefile.package
    sed -i -e 's|[^ \t]*gfortran*[^ \t] ||' ../Makefile.package
  fi

#  if (test -e ../Makefile.package.settings) then
#    sed -i -e '/^include.*cpl.*$/d' ../Makefile.package.settings
#  fi

fi
