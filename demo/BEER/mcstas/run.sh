#!/bin/bash
echo Running McStas
MCSTAS=/usr/share/mcstas/2.4.1
# MCSTAS_CC=gcc
mcmd="perl $MCSTAS/bin/mcrun.pl"

# run simulation MR
cfg="MR"
outpath=result_$cfg
rm -r $outpath
$mcmd  -n 10000 -p $cfg.par -d $outpath BEER_MCPL.instr

# run simulation MCB
cfg="MCB"
outpath=result_$cfg
rm -r $outpath
$mcmd  -n 10000 -p $cfg.par -d $outpath BEER_MCPL.instr
