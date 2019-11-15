!//////////////////////////////////////////////////////////////////////
!////  $Id: instruments.f,v 1.1 2008/07/30 21:21:31 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2008/07/30 21:21:31 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Instrument definitions
!////
!//////////////////////////////////////////////////////////////////////
      MODULE INSTRUMENTS
      USE TRACING
      use IO
      implicit none

! ID constants for instruments
      integer,parameter :: instr_bench=1
      integer,parameter :: instr_double=2
      integer,parameter :: instr_triple=3
      integer,parameter :: instr_tas=4

! ID variable for instruments
      integer :: INSTR_ID=instr_bench

      contains


      end MODULE INSTRUMENTS
