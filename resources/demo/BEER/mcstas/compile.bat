@echo off
set PATH=C:/mcstas-2.6/bin;%PATH%
set MCSTAS=C:/mcstas-2.6/lib
set MCSTAS_CC=gcc
set MCSTAS_TOOLS=C:/mcstas-2.6/lib/tools/Perl
rem set mcmd=perl C:/mcstas-2.6/lib/perlbin/mcrun.pl

set name=BEER_MCPL
if (%name%.==.) (
	echo "Edit the script and provide McStas file name (without extension)"
	exit
)
echo Compiling %name%
mcstas -o %name%.c %name%.instr
%MCSTAS_CC% -o %name%.exe %name%.c -I %MCSTAS%/libs/mcpl %MCSTAS%/libs/mcpl/libmcpl.a


pause
