# SIMRES - neutron ray-tracing simulation

SIMRES is an application for simulation of neutron beams using Monte Carlo ray-tracing method. It has been developed for instrument scientists as a tool for realistic simulation of neutron beam characteristics, optimization of instrument parameters, planning of experiments as well as simulation of synthetic data for validation of data analysis methods. 

Homepage: http://neutron.ujf.cas.cz/restrax 

User guide: http://neutron.ujf.cas.cz/restrax/download/simres-guide.pdf


-----------------------------------------------------------------
# TODO

## Terms of use
see COPYING.md

## Binary distribution

### Linux

### Windows

## Source distribution

### Get source files from git repository
- Clone SIMRES from GitHub repository:  
`git clone --recurse-submodules https://github.com/saroun/simres`

- Add PGPLOT source files:  
get source distribution at http://www.astro.caltech.edu/~tjp/pgplot, Version 5.2.2.
Unpack the archive to `./3dparty/pgplot`.

### Get source files from archive
Alternatively, it is possible to download source distribution as a zip archive at `http://neutron.ujf.cas.cz/restrax/download` and unpack it to a new folder. 


### Compile on Windows
Edit the `build.bat` script from the source distribution to make sure that the required toolchain is available. Specifically, the following tools are required:

To build SIMRES core:
- `perl` for makefile configuration
- `mingw-w64` compiler suite including `gfortran`

To build JSDRIV Windows server for PGPLOT:
- `lazbuild` command line compiler from Lazarus IDE. See submodules/jsdriv_server/README.md for details.

To build GUI:
- Java SDK (version 1.8 or later)
- `ant` builder
- `Java3D` package, see submodules/simresUI/README.md for instructions. 

Run `build.bat`. This will perform all necessary steps: configure, compile and install SIMRES to ./distr. At this stage, it should be possible to start SIMRES by executing `./distr/startGUI_win32.bat`. 

Use the command `perl Install.pl [target directory]` to install SIMRES. 

### Compile on Linux








