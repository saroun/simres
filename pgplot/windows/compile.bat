@echo off



echo "Current directory is %CD%"

if %1.==. (
	echo "Provide *full* path to PGPLOT source distribution as the 1st parameter"
	goto end:
)

if %2.==. (
	echo "Provide folder with makefile_gfortran and system dependent PGPLOT sources for mingw-w64 system"
	goto end:
)

set src=%1
set tdir=%2

if EXIST %src% (
    echo SRC=%src% 
) ELSE (
	echo "Source directory [%src%] does not exist"
	goto end:
)

set sysdir=%tdir%/../sys_mingw

if EXIST %sysdir% (
    echo SYSDIR=%sysdir% 
) ELSE (
	echo "System source directory [%sysdir%] does not exist"
	goto end:
)

make -C %tdir% -f makefile_gfortran SRC=%src%
make -C %tdir% -f makefile_gfortran SRC=%src% clean 
:end

