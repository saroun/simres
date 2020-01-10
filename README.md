# SIMRES - neutron ray-tracing simulation

SIMRES is an application for simulation of neutron beams using Monte Carlo ray-tracing method. It has been developed for instrument scientists as a tool for realistic simulation of neutron beam characteristics, optimization of instrument parameters, planning of experiments as well as simulation of synthetic data for validation of data analysis methods. 

Homepage: http://neutron.ujf.cas.cz/restrax   
SIMRES repository: https://github.com/saroun/simres  
User guide: http://neutron.ujf.cas.cz/restrax/download/simres-guide.pdf


-----------------------------------------------------------------

## Use conditions

SIMRES is provided under the terms of the *GNU General Public License v. 2.0* (see LICENSE), with the exception of the `3rd party software` listed below. 

### 3rd party software

`PGPLOT` - The Fortran graphics library written by Tim Pearson, California Institute of Technology, provides graphical representation of results.  
http://www.astro.caltech.edu/~tjp/pgplot  
*License*: see pgplot/copyright.notice

`Mersenne-Twister` - The random number generator developed by Makoto Matsumoto and Takuji Nishimura [ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, Januray pp.3-30, 1998] and coded to Fortran by Hiroshi Takano.  
http://www.math.sci.hiroshima-u.ac.jp/%7Em-mat/MT/emt.html  
*License*: see licenses/mersenne_twister.txt

`MCPL` - Monte Carlo Particle Lists library   
https://github.com/mctools/mcpl  
*License*: Creative Commons, licenses/mcpl.txt

`Java3D` - 3D Graphics Package  
http://www.java3d.org/
*License*: Sun Microsystems, see licenses/java3d.txt

`Jama` - Java Matrix Package.   
http://math.nist.gov/javanumerics/jama/  
*License*: public domain

## System requirements

Linux: 
- tested with Ubuntu 18.04 LTE
- `perl` (for running installation script only)
- `Java`, JRE or SDK (ver. 7 or higher), 64 bit
- `gcc`, `gfortran` (version >= 7), if you want to compile from sources

Windows:
- tested with Windows 10 ver. 1809 
- `perl` (for running installation script only)
- `Java`, JRE or SDK (ver. 7 or higher), 64 bit
- `mingw-w64` with gcc and gfortran (version >= 7), if you want to compile from sources
- `Lazarus` FreePascal, if you want to compile the PGPLOT driver (jsdriv_server) 






