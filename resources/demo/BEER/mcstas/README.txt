Update info 2020-02-17:
This example has been developed with McStas v2.4.1. Therefore it has to use some customised *.comp components. In later versions, these components are distributed with McStas. In any case, this example still seems to work with McStas 2.6. 
-------------------------------------------------
This is a demonstrator of the instrument file BEER_MCPL.instr, which describes secondary part of the diffractometer BEER@ESS (including the sample).
It employs MCPL_input.comp using a *.mcpl file as a source of neutrons in front of the sample. On the output, it yields:
- input beam map, spectrum and time/wavelength map
- sampling volume projections
- simulated diffractogram 
- time/2theta map (global and a detail focused at 2 steel reflections)
In the modulation mode, (MCP.par input parameters), it also performs data reduction of the modulated data, which removes the modulation from the diffractogram. A time/2theta map showing assumed single line, overlap and empty regions is also provided (based on the provided estimates of the dhkl values in the dhkl.dat file). 

Quick test: 
------------
On windows, you can run the run.bat script.
On Linux, run the ./run.sh script.

These batch files run 2 simulations using the parameters written in MR.par and MCB.par files and require MCPL files previously created by SIMRES (see tha associated SIMRES project).

*.mcpl files
--------------
The mcpl files are needed as input. Neutrons emerging form a 1x3 mm slit at 50 mm before the sample can be simulated by SIMRES using the associated SIMRES demo project (instrument file BEER_MCPL.xml and input files BEER_MCPL_MCB.inp and BEER_MCPL_MR.inp for the modulation and single pulse modes, respectively). Neutrons are collected at 10 mm behind the slit, i.e. 40 mm before the sample. They are normalized so that the sum of weighted neutrons equals to the neutronintensity in one ESS pulse at 5 MW. Required files defined in the *.par input files are: 

BEER_MR.mcpl: simulation in medium resolution mode with pulse shaping choppers, wavelength range centred at 2 AA. 

BEER_MCB.mcpl: simulation in the modulation mode, using the modulation chopper MCB with 8 slits (4 deg wide) rotating at 280 Hz. Wavelength range centred at 2 AA. 

NOTE: the files are supposed to contain a small number of events (~ 10000) to keep the file size small, therefore a large repetition factor is used in the example. For a real use, one should provide more neutrons in the input file for better statistics. The capability of MCPL export and an example of simulation of the primary part of the BEER diffractometer shall be available in SIMRES as of version 6.4 (assumed relase Jan 2019).  

*.par files
------------
McStas input parameters for BEER_MCPL.instr, employing the BEER_MR and BEER_MCB mcpl files, respectively. 

dhkl.tab
---------
Estimates of dhkl values (duplex steel), which are needed for data reduction in modulation mode (used by NPI_tof_dhkl_detector.comp).

*.laz files
------------
Input files for the Powdewr_N sample from the McStas distribution.

NPI_tof_dhkl_detector.comp
---------------------------
A detector component, which performs data reduction in event mode yielding the 1D diffractogram in both the pulse shaping and modulation modes. It also produces the overlap region map (the modulation mode only).

NPI_tof_theta_monitor.comp
---------------------------
A component derived from TOF_cylPSD_monitor, allowing to set angular limits. 
 
PowderN_upd.comp
---------------------------
Updated PowderN.comp component, permits to define focus target and focusing range for incoherent scattering. Improves statistics especially if the sample is rotated like in this example (a cylyndric sample oriented for axial strain measurement). 

