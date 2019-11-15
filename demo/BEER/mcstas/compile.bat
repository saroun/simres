@echo off
rem set PATH=C:/mcstas-2.4.1/bin;%PATH%
set MCSTAS=C:/mcstas-2.4.1/lib
set MCSTAS_CC=gcc
set MCSTAS_TOOLS=C:/mcstas-2.4.1/lib/tools/Perl
rem set mcmd=perl C:/mcstas-2.4.1/lib/perlbin/mcrun.pl

set name=%1
if (%name%.==.) (
	echo "Provide file name without extension as a prameter "
	exit
)
echo Compiling %name%
mcstas -o %name%.c %name%.instr
%MCSTAS_CC% -o %name%.exe %name%.c -I %MCSTAS%/libs/mcpl %MCSTAS%/libs/mcpl/libmcpl.a


