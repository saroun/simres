Example for running a combined SIMRES -> McStas -> Simres simulation
Jan Saroun, saroun@ujf.cas.cz
---------------------------------------------------------------------
Instrument configuration: ToF diffractometer BEER@ESS

1) The simulation starts by tracing the primary part of the  in medium resolution mode
(pulse shapinng choppers PSC1 and PSC2), with input slit 1x3 mm at 50 mm before the sample. This simulation runs in reverse tracing mode (sample -> source). At the end, neutrons are exported to a MCPL file (simres.mcpl) at the last component on the instrument list, placed at 10 mm after the input slit.

2) Neutrons are then taken over by McStas and propagated to the sample (PowderN) (by default using the atatched duplex.laz, duplex steel). After scattering in the sample, they are saved in mcstas.mcpl. McStas finishes  tracing of the secondary beam (radial collimator and detector).

3) Simres then loads back mcstas.mcpl at the same monitor position and traces neutrons through its version of the secondary beam (the same one, but made from SIMRES components).

This process is controled by a script, which can be found (or loaded) in the Script window:

#----------------
# Files
#---------------
SIMRES instrument file: ./BEER_MCPL_PWD.xml
Compiled McStas instrument file (Windows): ./mcstas/BEER_MCPL_PWD.exe
Source McStas instrument file: ./mcstas/BEER_MCPL_PWD.instr
Other SIMRES configuration data: placed in ./

#--------------
# How to start
#--------------
After installation, this folder is placed in the installation directory, which may have no write access for a user. It is then better first to go to File/Projects menu and make a clone of this project. All files
will be copied to the directory selected by user.
Then just execute the script below (should be also found in teh Script window and in ./BEER_MCPL_PWD.inp.
Let me know if it works :-)

#----------------------------------------------------------------------
# Commented script fot running combined SIMRES-McStas-SIMRES simulations.
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# NOTE 1:
# McStas instrument file for this example is provided in the ./mcstas subdirectory. 
# You have to compile this file and set the mcstas.file parameter below 
# to whatever is the instrument executable file name. 
# You can use full path to override the default mcstas directory. Compilation scripts 
# for Windows and Linux (McStas 2.5) are provided (edit them to match actual environment of your system). 
#
# NOTE 2:
# Start each command block with mcstas commands by "block ui" keyword.
# This will invoke the appropriate script interpreter.
# Switch to back  standard SIMRES interpreter by "block end' when running other SIMRES commands.
#
# NOTE 3:
# The McStas output directory is erased befor each start of McStas simulation (except of the mcpl files)
# You have to copy them manually somewhere else in order to preserve them.
# Alternatively, you can define a different output directory by using the script command
# mcstas outputdir
#----------------------------------------------------------------------


#----------------------------------------------------------------------
# STEP 0
# To switch on the McStas mode in SIMRES:
# a) In the SIMRES tracing option, set McStas = yes (click Apply).
# b) Fill in appropriate input and output distances. These are the distances between:
#  for input: last SIMRES monitor and McStas MCPL_input
#  for outout: McStas MCPL_output and the SIMRES monitor (where neutrons are loaded back)
#----------------------------------------------------------------------

#----------------------------------------------------------------------
# STEP 1
# Set up mcstas parameters (file, output directory, initialize, change instrument parameters as needed
# Use mcstas.list command to get actual instrument parameters printed in the console  window.
# Make sure the mcstas instrument parameters agree with those used in SIMRES.
# In this example, lam0 (wavelength) and pulse chopper distance (6.65 m = centre between PSC1 and PSC2) 
# and other affect the calculation of dhkl by the mcstas detector component.
#----------------------------------------------------------------------

# Setting wavelength in SIMRES
set INST LAMBDA 2.1
XML UPDATE

# Setting  up mcstas 
block ui

# executable file name
mcstas.file BEER_MCPL_PWD.exe

# directory for output, can be relative to SIMRES output directory
# this will be McStas current directory for execution of the instrument file 
mcstas.outdir mcresult

# initialize (run with -h option and get parameters)
mcstas.initialize

# Optional: change McStas instrument parameters
block ui
mcstas.lam0=2.1
mcstas.Lc=6.65
mcstas.Linp=0.04

# select a sample file, can be relative to SIMRES output directory:
mcstas.sample=./data/duplex.laz
# mcstas.sample=./data/Al.laz
# mcstas.sample=./data/SiO2_quartza
# mcstas.sample=./data/Ni.laz
# mcstas.sample=./data/Cu.laz
# mcstas.sample=./data/Ti.laz


# Optional: list the actual McStas parameters in the Console window
mcstas.list

#----------------------------------------------------------------------
# STEP 2
# Run combined SIMRES+McStas+SIMRES simulation with 10 MCPL repetitions.
# Number of neutrons in one cycle is given by the value of "tracing options / counts".
#----------------------------------------------------------------------
block ui
mcstas.runall 10

# McStas simulation is now finished, the results are stored in ./mcstas subdirectory 
# of the current project output folder (as set in  File / Projects).

#----------------------------------------------------------------------
# Plot resulting diffractogram (result of secondary beam simulation by SIMRES, not McStas !)
# Start with "block end" to make sure the standard script interpreter is invoked. 
#----------------------------------------------------------------------

block end
DO DET1D

#----------------------------------------------------------------------
# It is possible to repeat just the McStas part of the instrument simulation.
# use command mcstas.run with requested number of mcpl repetitions.
#----------------------------------------------------------------------

# block ui
# mcstas.run 20

