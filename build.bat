echo off
rem   Setup PATH and ANT_HOME to make required toolchain available
rem   See also submodules/simresUI/README.md for details about Java3D setting
rem
rem   Caution: if this batch is launched from a 32-bit file manager,
rem   then %programfiles%  translates to 'Program Files (x86)' instead of 'Program Files'
rem   even if Windows is 64bit. 


rem   Example setting: uncomment and edit as needed:
rem set ANT_HOME=\Java\apache-ant
rem set PATH=%programfiles%\Java\jdk-12.0.2\bin;%ANT_HOME%\bin;%PATH%
rem set PATH=C:\lazarus;%PATH%
rem set PATH=C:\mingw-w64\mingw64\bin;%PATH%

rem   Configure makefile (release options)
perl configure.pl mingw_w64
rem   Use following for compiling with debugging options:
rem perl configure.pl -dbg mingw_w64

rem   Compile SIMRES core
mingw32-make -f makefile 

rem   Compile windows server for PGPLOT
mingw32-make -f makefile jsdriv

rem   Compile GUI for SIMRES
mingw32-make -f makefile simresUI

rem   Cleanup 
mingw32-make -f makefile clean

rem   install in ./distr
mingw32-make -f makefile distr

pause
