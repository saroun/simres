echo off
set ANT_HOME=%programfiles%\Java\apache-ant
set PATH=%programfiles%\Java\jdk-12.0.2\bin;%ANT_HOME%\bin;%PATH%
set PATH=C:\lazarus;%PATH%
set PATH=C:\mingw-w64\mingw64\bin;%PATH%
rem "C:\WINDOWS\system32\cmd.exe"

perl configure.pl mingw_w64
mingw32-make -f makefile 
mingw32-make -f makefile jsdriv
mingw32-make -f makefile simresUI
mingw32-make -f makefile clean
perl Install.pl -dist

pause
