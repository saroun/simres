# Download SIMRES

Binary distributions for Linux and Windows platforms, as well as the guide for users are currently available at the CESNET owncloud [repository](https://owncloud.cesnet.cz/index.php/s/lcyKz7ceSkiqZYD).

For system requirements and detailed installation guide, refer to [INSTALL.md](INSTALL.md).
For running the program, you should have Java Runtime Environment (JRE) installed. If not found during the installation process, you can still edit the starting scripts and define the path to JRE later. 

## Linux
Using a terminal window, unpack the distribution and make it your current directory:

`tar -xzf simres-6.5.5-linux-x86_64.tar.gz`  
`cd simres-6.5.5-linux-x86_64`

Run the installer:

`perl Install.pl [target directtory]` 

## Windows
Unpack and run the installer as usually.  
You may need to set exceptions in your antivirus program if it blocks the installer or SIMRES (there have been issues with CyberCapture in Avast AVG).  
