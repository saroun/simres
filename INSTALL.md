# SIMRES - neutron ray-tracing simulation

SIMRES is an application for simulation of neutron beams using Monte Carlo ray-tracing method. It has been developed for instrument scientists as a tool for realistic simulation of neutron beam characteristics, optimization of instrument parameters, planning of experiments as well as simulation of synthetic data for validation of data analysis methods. 

Homepage: http://neutron.ujf.cas.cz/restrax  
SIMRES repository: https://github.com/saroun/simres  
User guide: http://neutron.ujf.cas.cz/restrax/download/simres-guide.pdf

-----------------------------------------------------------------
## CONTENTS

- [Terms of use](#terms-of-use)
- [Binary distributions](binary-distributions)
- [Building from the source](building-from-the-source)
- [Building of binary distribution packages](building-of-binary-distribution-packages)

## Terms of use
see COPYING.md

## Binary distributions

### Requirements

#### Java Runtime Environment

SIMRES GUI is build in Java and requires a 64-bit version of Java JRE or SDK to execute (minimum version is 1.8). On Linux, it can be installed e.g. as (Ubuntu, Debian):  
`sudo apt install openjdk-8-jre`  

On Windows, there are two options:  

- Get the instaler from https://www.java.com/en/download/manual.jsp and choose the Windows 64 bit offline version. Note that by default, this site offers only the 32 bit version so you need to choose the 64 bit installer manually. Note also that Oracle has recently changed their license conditions for end users so that you may not be eligible for using this JRE.  

- Choose the open source GPL-licensed JDK from https://openjdk.java.net/.Binary distributions should be available at https://jdk.java.net/. 

When Java JRE or SDK is installed at `Program Files\Java\[some name]`, the SIMRES installer should find it. Otherwise you need to provide a path to it by editting the JRE variable in `./simresGUI_win32.bat` manually.

#### Other Runtime Libraries
Read additional information provided with the download links. In general, the *Windows* binary distributions are build using the Mingw-w64 package and required runtime libraries are included in the distribution.  
On *Linux*, the GNU gcc package is used and required runtime libraries (in particular, `libgfortran`) should be installed by the system admin. This also depends on the version of gcc the SIMRES binaries were compiled against. At present, the binary distribution for Linux is built using the *Windows Linux Subsystem*. By Feb 2020, the version is gcc `Ubuntu 5.4.0-6ubuntu1~16.04.12`. 
If you plan to run SIMRES on a newer system, you may need to install `libgfortran.so.3`:   
`sudo apt-get install libgfortran3`
or build the SIMRES core from the source (see [below](#building-from-the-source)).

### Linux

- Get the binary package at `http://neutron.ujf.cas.cz/restrax/download` and unpack it.
- Run command  `[sudo] perl Install.pl [target directory]` from the distribution directory to install SIMRES. If the target directory is "`.`", SIMRES is installed to `/opt/simres` and an executable link is made to `/usr/bin/simres` (requires sudo privileges).
- Run the program by executing `[target directory]/simresGUI`. See the user guide at `[target directory]/doc/simres-guide.pdf` for more information.

### Windows

- Get the installer at `http://neutron.ujf.cas.cz/restrax/download` (the installer is packed in a zip archive to prevent Windows from blocking its execution).
- Execute the installer and follow instructions as usual.  
- Run the program (a launch icon should be on the desktop and in the Start menu. See the user guide for more information (a link is provided in the Start menu in Simres folder).

## Building from the source

### Get source files from git repository
- Clone SIMRES from GitHub repository:  
`git clone --recurse-submodules https://github.com/saroun/simres`  
Alternatively, it is possible to download source distribution as a zip archive at `http://neutron.ujf.cas.cz/restrax/download` and unpack it to a new folder. 

### Add 3rd party software

- Get PGPLOT source files:  
get source distribution at http://www.astro.caltech.edu/~tjp/pgplot, Version 5.2.2.
Unpack the archive to `./3dparty/pgplot`.

- Get Java3D library for your OS and architecture and make sure the submodule simresUI can see it. Refer to the Requirements section of ./submodules/simresUI/README.md for details.  

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

Run `build.bat`. This will perform all necessary steps: configure, compile and install SIMRES to ./distr. See the content of `build.bat` for make commands which can also be used to build the core, PGPLOT driver and GUI separately. After execution of the script, it should be possible to run SIMRES by executing `./distr/startGUI_win32.bat`. For building an installer, see the [last section](#"building-of-binary-distribution-packages") below. 



### Compile on Linux

TODO

## Building of binary distribution packages

### Windows
TODO

### Linux

- If you want to build a self-contained binary distribution, copy all required runtime libraries to the ./rtlib/windows directory. For example, when building with mingw-w64 on Windows, the libraries `libgcc_s_seh-1.dll` and 




