#!/bin/bash
# set up environment for your McStas distribution
MCSTAS=/usr/share/mcstas/2.5
MCSTAS_CC=gcc
INS=$1
if [ -z "$1" ]
then
	echo Provide file name without extension as a parameter.
	echo Trying BEER_MCPL_PWD.
	INS=BEER_MCPL_PWD
fi
else
   INS=$1
fi
echo Compiling $INS
mcstas -o $INS.c $INS.instr
$MCSTAS_CC -o $INS.exe $INS.c -lm -I $MCSTAS/libs/mcpl $MCSTAS/libs/mcpl/libmcpl.a
rm $INS.c

