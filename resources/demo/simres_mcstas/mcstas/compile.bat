@echo off
rem set up environment for your McStas distribution
set PATH=C:/mcstas-2.6/bin;%PATH%
set MCSTAS=C:/mcstas-2.6/lib
set MCSTAS_CC=gcc
set MCSTAS_TOOLS=C:/mcstas-2.6/lib/tools/Perl

set name=%1
if [%name%]==[] (
	echo Provide file name without extension as a parameter.
	echo Trying BEER_MCPL_PWD
	set name=BEER_MCPL_PWD
)
echo Compiling %name%

del %name%.c
mcstas -o %name%.c %name%.instr
%MCSTAS_CC% -o %name%.exe %name%.c -I %MCSTAS%/libs/mcpl %MCSTAS%/libs/mcpl/libmcpl.a
del /Q %name%.c 

:end

