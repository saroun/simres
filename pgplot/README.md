# PGPLOT binding for RESTRAX

This folder contains code needed to bind PGPLOT graphics library with RESTRAX.

- For PGPLOT, see http://www.astro.caltech.edu/~tjp/pgplot/
- Provided makefiles are configured for gcc/gfortran.
- On Windows, use the mingw-w64 compiler suite
- A native graphics server is required to use PGPLOT on Windows, see ./jsdriv_server/README.md

===============================================================

## Contains:

`copyright.notice`, PGPLOT copyright conditions

`./src`	
   - PGPLOT 5.2 SOURCE DISTRIBUTION SHOULD BE COPIED HERE

`./linux`
   - `make_gfortran`,  makefile to compile libpgplot.so 
   - `grexec.f`,  GREXEC function 
   - `compile.sh`,  compilation script for Linux

`./windows`
   - `compile.bat`, compilation script for Windows

`./jsdriv_server`
   - pgplot binding code + JSDRIV graphics device for windows 
   - see `README.md`
	
## Compilation quick guide:

1. Make sre you have PGPLOT ver. 5.2 source files copied in ./src
2. Change current directory according to your system and run the compilation script: 
   - on linux: build `libpgplot.so`
     - `cd ./linux`
     - `./compile.sh ../src`
   - on windows: build `libpgplot.dll`
	 - `cd ./windows`
	 - `compile.bat %CD%/../src  ../jsdriv_server/pgplot_binding/tgt `
	 - read `./jsdriv_server/README.md` for more information 
	 - build `jsdriv_server.exe` and `jsdrivlib.dll`, see `./jsdriv_server/README.md`
	   NOTE: These files are needed to run PGPLOT graphics on Windows, but they are not required for building libpgplot.dll
	 - NOTE 2: the 1st argument = path to `src` must be either absolute or relative to the 2nd argument. 
	   The 2nd argument = path to `tgt` can be absolute or relative to the current directory.

