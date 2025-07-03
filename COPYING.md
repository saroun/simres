# SIMRES - neutron ray-tracing simulation
Copyright (c) 1997-2025   
Nuclear Physics Institute, CAS, Rez, Czech Republic  
Institut Laue Langevin, Grenoble, France  

Authors:   
Jan Saroun, saroun@ujf.cas.cz  
Jiri Kulda, kulda@ill.eu

--------------------------------------------------------------

## Terms of use

The source code of SIMRES is distributed under the terms of the 
GNU General Public License ver. 2. The text of this license 
is provided in the file [LICENSE](LICENSE). 

However, note that different license conditions may apply 
for the 3rd party code, which is embedded in the SIMRES source 
distribution. Additionaly, SIMRES interfaces to other 3rd party modules 
which provide additional functionality. For details, see the   
section [3rd party software](#3rd-party-software) below.

Note also that substantial effort has been spent to develop SIMRES.
When you use this program in your work and publish results obtained 
with the help of it, please include a reference:

1. J. Saroun, J. Kulda, Physica B, 234-236, 1997, 1102-1104.
2. J. Saroun, J. Kulda, "Raytrace of Neutron Optical Systems with RESTRAX", 
in Modern Developments in X-Ray and Neutron Optics, 
eds. A. Erko, M. Idir, T. Krist, A.G. Michette, Springer Berlin 2008, p. 57-68.

------------------------------------------------------------------
## 3rd party software
See the license conditions provided with each source file/package.

### Code embeded in the SIMRES source distribution

`Mersenne-Twister` - The random number generator developed by Makoto Matsumoto and Takuji Nishimura [ACM Trans. on Modeling and Computer Simulation Vol. 8, No. 1, Januray pp.3-30, 1998] and coded to Fortran by Hiroshi Takano.  
http://www.math.sci.hiroshima-u.ac.jp/%7Em-mat/MT/emt.html  
*License*: see [licenses/mersenne_twister.txt](licenses/mersenne_twister.txt)

`MCPL` - Monte Carlo Particle Lists library   
https://github.com/mctools/mcpl   
*License*: Creative Commons, [licenses/mcpl.txt](licenses/mcpl.txt)

`Jama` - Java Matrix Package.   
http://math.nist.gov/javanumerics/jama/  
*License*: public domain

`napack` - Numerical linear algebra and optimization  
http://www.netlib.org/napack/  
*License*: public domain

### Software linked to SIMRES, distributed in binary packages

`PGPLOT` - The Fortran graphics library written by Tim Pearson, California Institute of Technology, provides graphical representation of results.  
http://www.astro.caltech.edu/~tjp/pgplot   
*License*: see [licenses/pgplot.txt](licenses/pgplot.txt)

`Java3D` - 3D Graphics Package  
https://jogamp.org/deployment/java3d/  
*License*: see [licenses/java3d.txt](licenses/java3d.txt)

