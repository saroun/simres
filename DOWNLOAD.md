# Download SIMRES

The latest release of SIMRES is also available as binary package for Windows and Linux at http://neutron.ujf.cas.cz/restrax/download.
For system requirements and installation guide, refer to [INSTALL.md](INSTALL.md).


## Linux
This is 64 bit binary distribution compiled with `gcc 7.5.0` on `Ubuntu 18.04` (using Windows Linux Subsystem). Download the archive  
http://neutron.ujf.cas.cz/restrax/download/simres-6.5.3-linux-x86_64.tar.gz  
unpack it and run `perl Install.pl [target directtory]`.

## Windows
The windows 64 bit binary distribution is compiled with `mingw-w64` (`gcc` and `gfortran`, ver. 8.1.0)). It also contains PGPLOT graphics device for Windows (https://github.com/saroun/jsdriv_server) compiled with Lazarus FreePascal. The user interface is compiled with Java SDk and requires JRE ver. 1.8 or higher. Download the installer packed in  
http://neutron.ujf.cas.cz/restrax/download/simres-6.5.3-win32-x86_64.zip  
Unpack and run the installer as usually.  
You may need to set exceptions in your antivirus program if it blocks the installer or SIMRES (there have been issues with CyberCapture in Avast AVG).  
