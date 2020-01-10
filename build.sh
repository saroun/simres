#!/bin/bash
#   Setup PATH and ANT_HOME to make required toolchain available
#   See also submodules/simresUI/README.md for details about Java3D setting
#

# verify tools:
cmds=( "perl" "make" "gcc" "gfortran" "javac" "ant" "java")
for c in "${cmds[@]}";
do
	command -v $c >/dev/null 2>&1 || { echo >&2 "$c required but not found. Aborting."; exit 1;}
done

if [ -z "$ANT_HOME" ];then
  echo "Warning: variable ANT_HOM is not defined. Ant is required to build simres GUI"
fi

#   Configure makefile (release options)
perl configure.pl gfortran
#   Use following for compiling with debugging options:
# perl configure.pl -dbg gfortran

#   Compile SIMRES core
make -f makefile 

#   Compile GUI for SIMRES
make -f makefile simresUI

#   Cleanup 
make -f makefile clean

#   install in ./distr
make -f makefile distr
