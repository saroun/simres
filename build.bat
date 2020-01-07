echo off
rem Setup PATH and ANT_HOME to make required toolchain available:
set ANT_HOME=%programfiles%\Java\apache-ant
set PATH=%programfiles%\Java\jdk-12.0.2\bin;%ANT_HOME%\bin;%PATH%
set PATH=C:\lazarus;%PATH%
set PATH=C:\mingw-w64\mingw64\bin;%PATH%
set MK=mingw32-make

rem Configure makefile
perl configure.pl mingw_w64

rem Compile SIMRES core
%MK% -f makefile 

rem Compile windows server for PGPLOT
%MK% -f makefile jsdriv

rem Compile GUI for SIMRES
%MK% -f makefile simresUI

rem Cleanup 
%MK% -f makefile clean

rem install in ./distr
%MK% -f makefile distr

pause
