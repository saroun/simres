!//////////////////////////////////////////////////////////////////////
!////  $Id: constants.f,v 1.8 2016/02/09 18:30:40 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.8 $
!////     $Date: 2016/02/09 18:30:40 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Physical constants
!////
!//////////////////////////////////////////////////////////////////////
      MODULE CONSTANTS
      implicit none

! *** Physical constants ***
! validated by NIST tables (2011)
      REAL(KIND(1.D0)), PARAMETER :: PI       = 3.141592653589793239D0
      REAL(KIND(1.D0)), PARAMETER :: TWOPI    = 6.283185307179586477D0
      REAL(KIND(1.D0)), PARAMETER :: SQRT2PI  = 2.5066282746310005024D0
      REAL(KIND(1.D0)), PARAMETER :: deg      = PI/180.D0
      REAL(KIND(1.D0)), PARAMETER :: minute   = PI/180.D0/60.D0
      REAL(KIND(1.D0)), PARAMETER :: sec      = PI/180.D0/3600.D0
      REAL(KIND(1.D0)), PARAMETER :: rad      = 180.D0/PI
      REAL(KIND(1.D0)), PARAMETER :: R8LN2    = 2.354820045D0
      REAL(KIND(1.D0)), PARAMETER :: HBAR     = 6.58211899D-4  !  [meV*ns]  Planck / 2PI
      REAL(KIND(1.D0)), PARAMETER :: HOVM     = 6.296223627D-1 !  [m/ms*A]  HBAR/m_n
      REAL(KIND(1.D0)), PARAMETER :: HSQOV2M  = 2.072124655D0  !  [meV*A^2] HBAR^2/(2.m_n)
      REAL(KIND(1.D0)), PARAMETER :: gammaL   = 1.83247185D8   !  [1/T/s]     gyromagnetic ratio
      REAL(KIND(1.D0)), PARAMETER :: GammaNi  = 1.731D-3       !  [rad/A]     critical angle for Ni nat. = qTNi/2PI
      REAL(KIND(1.D0)), PARAMETER :: qTNi     = 1.0876D-2      !  [1/A]       critical qT for Ni nat.
      REAL(KIND(1.D0)), PARAMETER :: k_B      = 8.6173324D-2    !  [meV*K^-1]  Boltzman constant
      REAL(KIND(1.D0)), PARAMETER :: G_EARTH  = 9.81D-9/HOVM**2   !  [mm/(mm*A)^2]  earth gravity constant / HOVM^2

! *** ray-tracing ***
      REAL(KIND(1.D0)), PARAMETER :: FLXNORM  = 1.D14/100      ! flux is normalized to 10^14/cm^2/sec

      end MODULE CONSTANTS
