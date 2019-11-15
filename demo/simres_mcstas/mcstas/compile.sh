#!/bin/bash
MCSTAS=/usr/share/mcstas/2.5
MCSTAS_CC=gcc

if [ -z "$1" ]
then
	echo Provide file name without extension as a parameter.
else
	echo Compiling $1
	mcstas -o $1.c $1.instr
	$MCSTAS_CC -o $1.exe $1.c -lm -I $MCSTAS/libs/mcpl $MCSTAS/libs/mcpl/libmcpl.a
fi


