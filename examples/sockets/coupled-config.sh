#!/bin/bash


# Path to examples
SOCKETS_PATH=$CPL_PATH/examples/sockets
EXAMPLES_PATH=$SOCKETS_PATH/examples

# CFD software 
SUPPORTED_CFD=( "OpenFOAM-3.0.1" )
declare -A URLS_CFD
URLS_CFD=( ["OpenFOAM-3.0.1"]="https://sourceforge.net/projects/foam/files/foam/3.0.1/OpenFOAM-3.0.1.tgz" )
declare -A MD5_CFD
MD5_CFD=( ["OpenFOAM-3.0.1"]="304e6a14b9e69c20989527f5fb1ed724" )
declare -A URLS_FOAM_3rdP
URLS_FOAM_3rdP=( ["OpenFOAM-3.0.1"]="https://sourceforge.net/projects/foam/files/foam/3.0.1/ThirdParty-3.0.1.tgz" )
declare -A MD5_FOAM_3rdP
MD5_FOAM_3rdP=( ["OpenFOAM-3.0.1"]="4665072d7d29ab9af5ced402f667185a" )
declare -A DOWNLOAD_CFD
DOWNLOAD_CFD=( ["OpenFOAM-3.0.1"]=wget )
declare -A DECOMPRESS_CFD
DECOMPRESS_CFD=( ["OpenFOAM-3.0.1"]=tgz )
declare -A CPL_SOCKET_CFD
CPL_SOCKET_CFD=( ["OpenFOAM-3.0.1"]="$CPL_PATH/examples/sockets/OpenFOAM/OpenFOAM-3.0.1" )

# MD software
SUPPORTED_MD=( "LAMMPS-dev" )
declare -A URLS_MD
URLS_MD=( ["LAMMPS-dev"]="https://github.com/lammps/lammps" )
declare -A MD5_MD
MD5_MD=( ["LAMMPS-dev"]="" )
declare -A DOWNLOAD_MD
DOWNLOAD_MD=( ["LAMMPS-dev"]=git )
declare -A DECOMPRESS_MD
DECOMPRESS_MD=( ["LAMMPS-dev"]="")
declare -A CPL_SOCKET_MD
CPL_SOCKET_MD=( ["LAMMPS-dev"]="$SOCKETS_PATH/LAMMPS/LAMMPS-dev" )
declare -A URLS_LAMMPS_PACKAGE
URLS_LAMMPS_PACKAGE=( ["USER-CPL"]="$SOCKETS_PATH/LAMMPS/LAMMPS-dev/cpl-socket/lammps_packages/USER-CPL")






# Pretty printing for "info", "error" and "warning" mensages.
pprint() {
    if [ "$2" == info ]; then
        printf "\e[1m\e[32m$1\e[0m"
    elif [ "$2" == info2 ]; then
        printf "\e[1m\e[35m$1\e[0m"
    elif [ "$2" == error ]; then
        printf "\e[1m\e[31m$1\e[0m"
    elif [ "$2" == warning ]; then
        printf "\e[1m\e[33m$1\e[0m"
    fi
}

exit_error() {
    echo
    pprint "Error:\n" error
    pprint "   $1\n" error
    exit 1
}

# This function download, md5sum and extract files from
download() {
    name=$1
    url=$2
    fname=$2
    fname=${fname##*/}
    download_m=$3
	# Optional parameters
    md5_hash=${4:-""}
    decompress_m=${5:-""}

    # Download version according to the specified method
    if [ $download_m == wget ]; then
        pprint "[Downloading ${name}...]\n" info
        wget $url --no-check-certificate
        [ $? -ne 0 ] && exit_error "Download of $fname failed!\n" || pprint "[Download... Done]\n" info
        if [ $md5_hash != "" ]; then
            md5sum $fname | grep $md5_hash --quiet
            [ $? -ne 0 ] && exit_error "MD5 check failed! Download again $fname.\n" || pprint "*MD5 check... [OK]\n" info
        fi
        # Decompress if necessary
        if [ $decompress_m == tgz ]; then
            pprint "*Decompressing files... " info
            tar -xzf $fname
            [ $? -ne 0 ] && exit_error "[FAIL]\n" || pprint "[OK]\n" info
        fi
    elif [ $download_m == git ]; then
        pprint "[Downloading ${name}...]\n" info
        git clone $url $name
        [ $? -ne 0 ] && exit_error "Download of $fname failed!\n" || pprint "[Download... Done]\n" info
    elif [ $download_m == cp ]; then
        pprint "[Downloading ${name}...]\n" info
        cp -R $url $name
        [ $? -ne 0 ] && exit_error "Download of $name failed!\n" || pprint "[Download... Done]\n" info

    fi
}
    

    

case $1 in
create)
    echo
    pprint "****************************************************\n" info2
    pprint "[--[CONFIGURATION SCRIPT FOR COUPLED SIMULATIONS]--]\n" info2
    pprint "****************************************************\n" info2
    echo
    pprint "[*] Select supported CFD:\n" info
    index=0
    for item in ${SUPPORTED_CFD[*]}
    do
        printf "   [%d] %s\n" $index $item
        index+=1
    done
    printf "Option: "
    read $opt_cfd
    cfd_name=${SUPPORTED_CFD[$opt_cfd]}
    echo
    pprint "[*] Select supported MD:\n" info
    index=0
    for item in ${SUPPORTED_MD[*]}
    do
        printf "   [%d] %s\n" $index $item
        index+=1
    done
    printf "Option: "
    read $opt_md
    md_name=${SUPPORTED_MD[$opt_md]}
    ROOT_DIR="${cfd_name}_${md_name}"
    mkdir $ROOT_DIR
    [ $? -ne 0 ] && exit_error "Directory $ROOT_DIR already exists."
    cd $ROOT_DIR
    MD_DIR="${md_name}_coupled"
    CFD_DIR="${cfd_name}_coupled"

    
    echo
    pprint "[---[Downloading Sockets]---]\n" info2
    echo
    pprint "[1]" info
    download $MD_DIR ${CPL_SOCKET_MD[$md_name]} cp
    echo
    pprint "[2]" info
    download $CFD_DIR ${CPL_SOCKET_CFD[$cfd_name]} cp

    cd $CFD_DIR
    echo
    pprint "[---[Downloading CFD and MD software]---]\n" info2
    echo
    url_cfd=${URLS_CFD[$cfd_name]}
    download_m_cfd=${DOWNLOAD_CFD[$cfd_name]}
    md5_hash_cfd=${MD5_CFD[$cfd_name]}
    decompress_m_cfd=${DECOMPRESS_CFD[$cfd_name]}
    pprint "[1]" info
    download $cfd_name $url_cfd $download_m_cfd $md5_hash_cfd $decompress_m_cfd
    # If the CFD is an OpenFOAM version, the Third party libraries have to be downloaded.
    if [[ "$cfd_name" =~ "FOAM" ]]; then
        url_3rdP=${URLS_FOAM_3rdP[$cfd_name]}
        md5_hash_3rdP=${MD5_FOAM_3rdP[$cfd_name]}
        name_3rdP="3rdPartyLibraries(OpenFOAM)"
        echo
        pprint "[1.2]" info
        download ${name_3rdP} $url_3rdP  wget $md5_hash_3rdP tgz
    fi
    cd ..

    cd $MD_DIR
    url_md=${URLS_MD[$md_name]}
    download_m_md=${DOWNLOAD_MD[$md_name]}
    md5_hash_md=${MD5_MD[$md_name]}
    decompress_m_md=${DECOMPRESS_MD[$md_name]}
    echo
    pprint "[2]" info
    download $md_name $url_md $download_m_md $md5_hash_md $decompress_m_md
    if [[ "$md_name" =~ "LAMMPS" ]]; then
        cd $md_name/src
        echo 
        pprint "[2.2]" info
        download "USER-CPL" ${URLS_LAMMPS_PACKAGE["USER-CPL"]} cp
        cd ../..
    fi
    cd ..

    echo
    pprint "[---[Downloading related examples]---]\n" info2
    # This assumes example names have to be formated as MDname-CFDname 
    example_md=$(echo $md_name | cut -f 1 -d '-')
    example_cfd=$(echo $cfd_name | cut -f 1 -d '-')
    example_name="${example_md}-${example_cfd}"
    pprint "[3]" info
    download "examples" "${EXAMPLES_PATH}/${example_name}/" cp


    pprint "[Success!]\n" info
    ;;
esac
