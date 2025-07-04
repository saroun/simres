# How to compile and install SIMRES

## CONTENTS

- [Terms of use](#terms-of-use)
- [System requirements](#system-requirements)
- [Binary distributions](#binary-distributions)
- [Building from the source](#building-from-the-source)
- [Building of binary distribution packages](#building-of-binary-distribution-packages)

## Terms of use
see [COPYING.md](COPYING.md)

## System requirements

### At runtime

*Linux*: 
- tested with Ubuntu 22.04 LTE
- `java JRE`, tested with verisons <= 11
- `perl` (only for running installer script)

*Windows*:
- Windows 10+
- `java JRE`, tested with verisons <= 11

### Build toolchain

*Linux*: 
- tested with Ubuntu 22.04 LTE
- `java SDK`, tested with java-11-openjdk
- `perl` for installation and build scripts
- `gcc`, `gfortran`, tested with versions <= 9
- `ant` for compilation of CLI and GUI
- `3d party software` (see [Building of binary distribution packages](#building-of-binary-distribution-packages))
    - PGPLOT
    - Java3D

*Windows*:
- Windows 10+
- `java SDK`, tested with Oracle JDK 11
- `perl` (for the configure script)
- `mingw-w64` with gcc and gfortran (tested with version 8.1.0)
- `Lazarus` FreePascal, if you want to compile the PGPLOT graphics driver for Windows 
(jsdriv_server) 
- `Inno Setup` for making the Windows installer


## Binary distributions

### Requirements

#### Java Runtime Environment

SIMRES CLI and GUI are build in Java and require a 64-bit version of JRE or JDK to execute (minimum version is 1.8). On Linux, it can be installed e.g. as (Ubuntu, Debian):  
`sudo apt install openjdk-11-jre`  

On Windows, there are two options:  

- Get the JDK instaler from Oracle https://www.oracle.com/java/technologies/downloads/#java11  

- Choose the open source GPL-licensed JDK from https://openjdk.java.net/. Binary distributions should be available at https://jdk.java.net/. 

**NOTE** 
There are incompatibility issues related to the Java3D package when using Java version > 13. 

The SIMRES installers will try to find the latest version installed on your system. If it fails, it is possible to provide the path to the `java` executable by editting the JRE variable in the start scripts `simres.bat` (Windows) or `simres` (Linux).

#### Other Runtime Libraries
The *Windows* binary distributions are build using the Mingw-w64 package and required runtime libraries should be included in the binary distribution.  
On *Linux*, the GNU gcc and gfortran are used and required runtime libraries (in particular, `libgfortran`) should be installed by the system admin. This also depends on the version of gcc the SIMRES binaries were compiled against. At present, the binary distribution for Linux is built on Ubuntu 22.04 with gfortran version 9.1.0, therefore `libgfortran5` is required.

### Linux

- Get the binary package (see [DOWNLOAD.md](DOWNLOAD.md)) and unpack it.
- Run command  `[sudo] perl Install.pl [target directory]` from the distribution directory to install SIMRES. If the target directory is "`.`", SIMRES is installed to `/opt/simres` and an executable link is made to `/usr/bin/simres` (requires sudo privileges).
- Test the CLI by executing `[target directory]/simres -test 0` (or just `simres -test 0` for the standard installation). 
- Run the GUI by running just `simres` (or `[target directory]/simres`)

See the user guide at `[target directory]/doc/simres-guide.pdf` for more information.

### Windows

- Get the installer (see [DOWNLOAD.md](DOWNLOAD.md)).
- Execute the installer and follow instructions.  
- Run the program (a launch icon should be on the desktop and in the Start menu. See the user guide for more information (a link is provided in the Start menu in the Simres folder).
- For testing the program in command mode, open the Simres command window (a link should be available in the program group). Then execute  e.g. `simres -test 0`.
- Some antivirus programs like AVG may hinder the first launch of the program. If this happens, press RESET on the control panel to restart the kernel.

## Building from the source

### Get source files from git repository
- Clone SIMRES from GitHub repository:  
`git clone --recurse-submodules https://github.com/saroun/simres`  
Alternatively, it is possible to download source distribution from the GitHub site as a zipped archive. 

### Add 3rd party software

- Get `PGPLOT` source files:  
get the source distribution at http://www.astro.caltech.edu/~tjp/pgplot, Version 5.2.2.
Unpack the archive to `./3dparty/pgplot`.

- Get the `Java3D` library for your OS and architecture and make sure the submodule simresUI can see it. Refer to the Requirements section of [./submodules/simresUI/README.md](./submodules/simresUI/README.md) for details.  

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

Run `build.bat`. This will perform all necessary steps: configure, compile and install SIMRES to ./distr. See the content of `build.bat` for the make commands which can also be used to build the core, PGPLOT driver and GUI separately. After execution of the script, it should be possible to run SIMRES by executing `.\distr\simres.bat`. For running a test in command mode, execute e.g. `.\distr\simres.bat -test 0`. For building an installer, see the [last section](#"building-of-binary-distribution-packages") below. 

### Compile on Linux

Make sure that the required toolchain is available. Specifically, the following tools are required:
- gcc compilers suite with gfortran  
`sudo apt-get install gcc gfortran make`
- libX11 headers  
`sudo apt-get install libx11-dev`
- Java development kit  
`sudo apt-get install openjdk-11-jdk-headless`
- Ant builder  
`sudo apt-get install ant`

Run the script `build.sh`. This will perform all necessary steps: configure, compile and install SIMRES to ./distr. See the content of `build.sh` for the make commands which can also be used to build the core, PGPLOT driver and GUI separately. After execution of the script, it should be possible to run SIMRES by executing `./distr/simres`.  For running a test in command mode, execute e.g. `./distr/simres -test 0`. For building a binary package, see the [last section](#"building-of-binary-distribution-packages") below. To install, run 

`perl Install.pl [target directory]`

## Building of binary distribution packages

### Windows

- If you want to build a self-contained binary distribution, copy all required runtime libraries to the `./rtlib/windows` subdirectory. For example, when building with mingw-w64 on Windows, the libraries `libgcc_s_seh-1.dll` and `libquadmath-0.dll` are required. 
- After building and testing the distribution (see [Compile on Windows](#compile-on-windows)), run the script ZipBin.pl. Add the "-inno" option to build Windows installer using the INNO Setup utility `iscc`:  
`perl ZipBin.pl -inno`

### Linux

- After building and testing the distribution (see [Compile on Linux](#compile-on-linux)), run the script ZipBin.pl:  
`perl ZipBin.pl`

- To create source distribution, run:  
`perl ZipSrc.pl`



