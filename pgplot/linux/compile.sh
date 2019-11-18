#!/bin/bash

if [ -z "$1" ]; then
    echo "Usage:"
	echo
	echo "./compile.sh [all,prog] SRC"
	echo "       SRC = path to PGPLOT source distribution"
	echo
	echo "./compile.sh [clean,erase]"
	echo
	echo "./compile.sh install TGT PGD "
	echo "       TGT = folder for binary files"
	echo "       PGD = folder for grfont.dat etc. (what should be PGPLOT_DIR variable)"
	echo
	exit
fi

tgt="$1"
if [[ $tgt ==  "all" || $tgt == "prog" ]]; then
	if [ ! -d "$2" ]; then
		echo "The directory $1 does not exist"
		exit
	fi
	# quick check of the directory content
	if [ ! -d "$2/sys_linux" ]; then
		echo "The directory $2 is not the PGPLOT source directory."
		exit
	fi
	make -f makefile_gfortran $tgt SRC="$2"
elif [[ $tgt == "clean" || $tgt == "erase" ]]; then
    echo "target=$tgt"
	make -f makefile_gfortran $tgt
elif [ $tgt == "install" ]; then
    echo "target=$tgt [$2] [$3]"
	if [ ! -d "$2" ]; then
		echo "The directory $2 does not exist"
		exit
	fi
	if [ ! -d "$3" ]; then
		echo "The directory $3 does not exist"
		exit
	fi
	make -f makefile_gfortran $tgt BIN="$2" PGD="$3"
else
	echo "Unknown target: $1"
fi



