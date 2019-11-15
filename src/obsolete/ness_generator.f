!//////////////////////////////////////////////////////////////////////
!////  $Id: ness_generator.f,v 1.1 2007/10/08 08:35:37 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2007/10/08 08:35:37 $
!//////////////////////////////////////////////////////////////////////
!////
C////  Subroutines describing objects - GENERATOR
C////  special component used to generate initial event object (NEUTRON)
!////
!//////////////////////////////////////////////////////////////////////


C-------------------------------------------------
      Subroutine GENERATOR_DEFAULT(OBJ)
C sets default generator properties (all dimensions free)
C-------------------------------------------------
      USE COMPONENTS
      implicit none
      TYPE(GENERATOR) :: OBJ
        OBJ%SHAPE=0
        OBJ%HPROF=0
        OBJ%VPROF=0
        OBJ%KPROF=0
      end Subroutine GENERATOR_DEFAULT

C-------------------------------------------------
      Subroutine GENERATOR_INIT(OBJ)
C-------------------------------------------------
      USE COMPONENTS
      implicit none
      TYPE(GENERATOR) :: OBJ
        OBJ%COUNT=0
        OBJ%SCOUNT=0.D0
      end Subroutine GENERATOR_INIT

C--------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION GENERATOR_VOL(OBJ)
C return sampling volume - only for dimensions with applied limits
C--------------------------------------------------------------
      USE COMPONENTS
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'randvars.inc'
      TYPE(GENERATOR) :: OBJ
      real(KIND(1.D0)) :: PP,COEF
  ! impose limits
      PP=1.D0
      COEF=1.D0
  ! width x height
      if (OBJ%SHAPE.NE.0) then
        PP=PP*OBJ%WIDTH*OBJ%HEIGHT
        select case(OBJ%SHAPE)
        case(2)
          COEF=COEF*PI/4.D0
        end select
      endif
  ! horizontal divergence
      if (OBJ%HPROF.NE.0) then
        PP=PP*OBJ%DIVH
        select case(OBJ%HPROF)
        case(2)
          COEF=COEF*0.5
        case(3)
          COEF=COEF*R8LN2
        end select
      endif
  ! vertical divergence
      if (OBJ%VPROF.NE.0) then
        PP=PP*OBJ%DIVV
        select case(OBJ%VPROF)
        case(2)
          COEF=COEF*0.5
        case(3)
          COEF=COEF*R8LN2
        end select
      endif
  ! wavelength
      if (OBJ%KPROF.NE.0) then
        PP=PP*OBJ%DK
        select case(OBJ%KPROF)
        case(2)
          COEF=COEF*0.5
        case(3)
          COEF=COEF*R8LN2
        end select
      endif
      GENERATOR_VOL=PP*COEF
      end FUNCTION GENERATOR_VOL



C-------------------------------------------------
      LOGICAL*4 FUNCTION GENERATOR_GO(OBJ,NEU)
C-------------------------------------------------
      USE COMPONENTS
      USE TRACING
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'randvars.inc'

      TYPE(GENERATOR) :: OBJ
      TYPE(NEUTRON) :: NEU
      real(KIND(1.D0)) :: PP,X,Y
      REAL :: RAN1

! fill NEU
      NEU.R(1)=XRND(1)
      NEU.R(2)=XRND(2)
      NEU.R(3)=0.0
      NEU.K(1)=XRND(3)*OBJ%K0
      NEU.K(2)=XRND(4)*OBJ%K0
      NEU.K(3)=XRND(5)*OBJ%K0+OBJ%K0
      NEU.K0=SQRT(NEU.K(1)**2+NEU.K(2)**2+NEU.K(3)**2)
      NEU.S(1:3)=0.D0
      NEU.S(1)=2*NINT(RAN1())-1
      NEU.T=0.D0
1     format(a,7(G10.4,1x))
  ! impose limits
      PP=1.D0
  ! width x height
      if (OBJ%SHAPE.NE.0) then
        X=ABS(NEU.R(1))/OBJ%WIDTH
        Y=ABS(NEU.R(2))/OBJ%HEIGHT
        select case(OBJ%SHAPE)
        case(1)
          if (X.gt.0.5) PP=0.D0
          if (Y.gt.0.5) PP=0.D0
        case(2)
          if ((X**2+Y**2).gt.0.25) PP=0.D0
        end select
      endif
c      write(*,1) 'SHAPE: ',PP
  ! horizontal divergence
      if (OBJ%HPROF.NE.0) then
        X=ABS(NEU.K(1)/NEU.K0/OBJ%DIVH)
        select case(OBJ%HPROF)
        case(1)
          if (X.GT.0.5) PP=0.D0
        case(2)
          PP=PP*MAX(0.D0,1.D0-X)
        case(3)
          if (X.GT.3) then
            PP=0.D0
          else
            PP=PP*exp(-0.5*X**2)
          endif
        end select
      endif
c      write(*,1) 'HDIV: ',PP
  ! vertical divergence
      if (OBJ%VPROF.NE.0) then
        X=ABS(NEU.K(2)/NEU.K0/OBJ%DIVV)
        select case(OBJ%VPROF)
        case(1)
          if (X.GT.0.5) PP=0.D0
        case(2)
          PP=PP*MAX(0.D0,1.D0-X)
        case(3)
          if (X.GT.3) then
            PP=0.D0
          else
            PP=PP*exp(-0.5*X**2)
          endif
        end select
      endif
c      write(*,1) 'VDIV: ',PP
  ! wavelength
      if (OBJ%KPROF.NE.0) then
        X=abs(NEU.K0-OBJ%K0)/OBJ%K0/OBJ%DK
        select case(OBJ%KPROF)
        case(1)
          if (X.GT.0.5) PP=0.D0
        case(2)
          PP=PP*MAX(0.D0,1.D0-X)
        case(3)
          if (X.GT.3) then
            PP=0.D0
          else
            PP=PP*exp(-0.5*X**2)
          endif
        end select
      endif

  ! weight
      NEU.P=PP
c      write(*,1) 'GEN: ',NEU.R(1:2),NEU.K,PP
c      read(*,*)
      if (PP.GT.1.D-10) then
        OBJ%SCOUNT=OBJ%SCOUNT+1.D0
        OBJ%COUNT=OBJ%COUNT+1
        GENERATOR_GO=.TRUE.
      else
        GENERATOR_GO=.FALSE.
      endif

      END
