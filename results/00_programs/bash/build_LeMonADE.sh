#!/bin/bash
set -e

## Matthew A. Dorsey
## Institute for Theory of Polymer
## 03.06.2024
## linux script for automating compiling of LeMonADE package with executables

## PARAMETERS
# default path to lemonade library, if not specified bu user
DEFAULT_LEOMONADE_PATH='/home/users/dorsey/Desktop/IPF/LeMonADE_install/'

## FUNCTIONS
# none


## OPTIONS
while getopts "l:" option; do
	case $option in
		l) # adjust default LeMonADE library
			DEFAULT_LEOMONADE_PATH="${OPTARG}" ;;
        \?) # sonstiges
            exit $NONZERO_EXITCODE
    esac
done


## ARGUMENTS
# none


## SCRIPT
# build lemonade library
if test -d build; then
    rm -r build
fi

mkdir build
cd build
cmake -DLEMONADE_DIR="${DEFAULT_LEOMONADE_PATH}" ..
make
cd ..
