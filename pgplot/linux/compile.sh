#!/bin/bash

if [ -z "$1" ]; then
    echo "provide PGPLOT source distribution directory"
    echo "usage: ./compile.sh [source directory]"
	exit
fi

if [ ! -d "$1" ]; then
    echo "The directory $1 does not exist"
	exit
fi

# quick check of the directory content
if [ ! -d "$1/sys_linux" ]; then
    echo "The directory $1 is not the PGPLOT source directory."
	exit
fi

make -f makefile_gfortran SRC=$1
make -f makefile_gfortran SRC=$1 clean

