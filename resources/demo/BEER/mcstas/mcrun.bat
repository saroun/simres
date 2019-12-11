@echo off
rem set PATH=C:/mcstas-2.4.1/bin;%PATH%
set MCSTAS=C:/mcstas-2.4.1/lib
set MCSTAS_CC=gcc
set MCSTAS_TOOLS=C:/mcstas-2.4.1/lib/tools/Perl
set mcmd=perl C:/mcstas-2.4.1/lib/perlbin/mcrun.pl


rem # run simulation for MR setop
set outpath=result_MR
rmdir /S/Q %outpath%
%mcmd% -n 10000 -p MR.par -d %outpath% BEER_MCPL.instr
pause

rem # run simulation for MCB setop
set outpath=result_MCB
rmdir /S/Q %outpath%
%mcmd% -n 10000 -p MCB.par -d %outpath% BEER_MCPL.instr
pause



