# S I M R E S - Neutron ray-tracing simulations
Copyright (c) 1997-2020   
Nuclear Physics Institute, CAS, Rez, Czech Republic  
Institut Laue Langevin, Grenoble, France  
http://neutron.ujf.cas.cz/restrax/  

Authors:   
Jan Saroun, saroun@ujf.as.cz  
Jiri Kulda, kulda@ill.eu  

--------------------------------------------------------------

## Terms of use

The source code of SIMRES is distributed under the terms of the 
GNU General Public License ver. 2. The text of this license 
is provided in the file LICENSE. 

However, note that different license conditions may apply 
for the 3rd party code, which is embedded in the SIMRES source 
distribution. Additionaly, SIMRES interfaces to 3rd party modules 
which provide additional functionality. For details, see the   
section *3rd party software* below.

Note also that substantial effort has been spent to develop SIMRES.
When you use this program in your work and publish results obtained 
with the help of it, please include a reference:

J. Saroun, J. Kulda, Physica B, 234-236, 1997, 1102-1104.


J. Saroun, J. Kulda, "Raytrace of Neutron Optical Systems with RESTRAX", 
in Modern Developments in X-Ray and Neutron Optics, 
eds. A. Erko, M. Idir, T. Krist, A.G. Michette, Springer Berlin 2008, p. 57-68.

------------------------------------------------------------------
## 3rd party software
See the license conditions provided with each source file/package.

### Code embeded in the SIMRES source distribution

`Jama` - Java Matrix Package.   
embedded in submodules/simresUI/simresCON/src/Jama  and GUI/simresCON.jar
http://math.nist.gov/javanumerics/jama/  
*License*: public domain

`napack` - Numerical linear algebra and optimization
embedded in src/napack.f
http://www.netlib.org/napack/  
*License*: public domain

`Mersenne-Twister` - The random number generator 
Written by Makoto Matsumoto and Takuji Nishimura and coded to Fortran by Hiroshi Takano.  
embedded in src/mtmod.f
http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/VERSIONS/FORTRAN/mtfort90.f  
*License*: LGPL 

`MCPL` - Monte Carlo Particle Lists library 
embedded in submodules/mcpl  
https://github.com/mctools/mcpl  
*License*: Creative Commons

### Software linked to SIMRES, distributed in binary packages

`PGPLOT` - The Fortran graphics library written by Tim Pearson, California Institute of Technology, provides graphical representation of results.  
http://www.astro.caltech.edu/~tjp/pgplot   
*License*: see licenses/pgplot.txt

`Java3D` - 3D Graphics Package  
http://www.java3d.org  
*License*: see licenses/java3d.txt

