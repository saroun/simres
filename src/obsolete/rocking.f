!//////////////////////////////////////////////////////////////////////
!////  $Id: rocking.f,v 1.1 2008/07/30 21:21:31 saroun Exp $
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
!////  Module for simulation of a crystal rocking curve
!////
!//////////////////////////////////////////////////////////////////////
      MODULE ROCKING
      USE TRACINGDATA
      USE COMPONENTS
      use RESMAT
      USE TRACING
      USE REPORTS
      implicit none

      CONTAINS

C--------------------------------------------------------
      SUBROUTINE ROCKING_INI(CRYST,DIVH,DIVV)
C Define beamline for rocking curve measurement
C only a generator and crystal
C K0    ... mean wave vector magnitude
C DIVH  ... horizontal divergence [rad]
C DIVV  ... vertical divergence [rad]
C--------------------------------------------------------
      include 'const.inc'
      type(CRYSTAL),intent(in) :: CRYST
      REAL(KIND(1.D0)),intent(in) :: DIVH,DIVV
      REAL(KIND(1.D0)) :: K0
      integer :: i

! set defaults first
!      call TRACING_DEFAULTS
! uncomment the following to override default settings:
  ! down-stream tracing (set 1 for up-stream)
c      TRACING_MODE=1
  ! switch event logs
c      EventLogs=.TRUE.
  ! switch variance reduction
c      TRACING_VARI=.TRUE.
  ! swicth max. value optimization
c      TRACING_MAXV=.FALSE.
  ! reporting of tracing statistics
      REPOPT%STATS=.FALSE.
  ! reporting of sampling volume
      REPOPT%VOL=.FALSE.
  ! reporting of variance reduction progress
      REPOPT%VARI=.FALSE.
  ! reporting of tracing progress
      REPOPT%PROG=.FALSE.
  ! reporting basic results (intensity, E-spread)
      REPOPT%RES=.FALSE.
  ! report max. value optimization
      REPOPT%MAXV=.FALSE.

      TROPT%ELIMIT=1.D2

  ! clean beamline
      call DeleteComponents
  ! add crystal
      i=AddComponent(PCOMOBJ(CRYST),1)
      ACRYSTALS(1)%FRAME%DIST=100.D0 ! default source distance 100 mm

 ! prepare generator
  ! set narrow limits in all variables => parallel monochromatic beam
  ! when required, allow extended divergence
      call GENERATOR_SETDIV(max(1.D-6,DIVH),max(1.D-6,DIVV))
      call GENERATOR_SETSHAPE(1,1.D-3,1.D-3)
  ! addCRYSTAL calls CRYST_INI on the new object, so that we can use it's LAMBDA value now
      if (ACRYSTALS(1)%LAMBDA.gt.0) then
        K0=2.D0*PI/ACRYSTALS(1)%LAMBDA
      else
        K0=2.D0*PI/1.8
      endif
      call GENERATOR_SETK(K0,1.D-6)

  !    CALL WRITE_BEAMLINE(20)

      END SUBROUTINE ROCKING_INI

C------------------------------------------------------------------------
      SUBROUTINE ROCKING_CURVE(CR,NC,NTH,DTH,DIVH,DIVV,RTH,DRTH)
C measures rocking curve of a crystal
C INPUT:
C   CR  .. crystal object
C   NC  .. number of counts/point
C   NTH .. number of steps
C   DTH .. step in [rad]
C   DIVH,DIVV .. optional divergence of incident beam in [rad]
C OUTPUT:
C   RTH(NTH)  .. array with resulting reflectivity curve
C------------------------------------------------------------------------
      implicit none
      INCLUDE 'const.inc'
      TYPE(CRYSTAL),intent(in) :: CR
      integer, intent(in) :: NC,NTH
      REAL(KIND(1.D0)),intent(in) :: DTH,DIVH,DIVV
      REAL(KIND(1.D0)),intent(out) :: RTH(NTH),DRTH(NTH)

      INTEGER :: ISTEP
      REAL(KIND(1.D0)) :: gon3,theta,Z,DINC
1     format(a,': ',4(G10.4,2x))

C count target
      TROPT%CNT=NC
  ! start scan with rocking angle GON(3)
      gon3=CR%FRAME%GON(3)
! simulate rocking curve
      call ROCKING_INI(CR,DIVH,DIVV)
      call XML_PROGRESS(SMES,0,NTH,0,' ')
      DO ISTEP=1,NTH
        theta=(-(NTH-1)/2.+ISTEP-1)*DTH
        ACRYSTALS(1)%FRAME%GON(3)=theta+gon3
        call TRACING_PROC
        RTH(ISTEP)=TRACING_INT
        DRTH(ISTEP)=TRACING_DINT
        call XML_PROGRESS(SMES,1,NTH,ISTEP,' ')
      ENDDO
      call XML_PROGRESS(SMES,2,NTH,NTH,' ')
! delete crystal, simulate incident intensity and normalize reflectivity
      call DeleteComponents
      TROPT%CNT=NC
      call TRACING_PROC
c      write(*,2) 0.0,TRACING_INT,TRACING_DINT
      if (TRACING_INT.gt.0) then
        DINC=TRACING_DINT/TRACING_INT
        DO ISTEP=1,NTH
          Z=DINC**2
          if (RTH(ISTEP).GT.0) Z=Z+(DRTH(ISTEP)/RTH(ISTEP))**2
          RTH(ISTEP)=RTH(ISTEP)/TRACING_INT
          DRTH(ISTEP)=MAX(DINC,RTH(ISTEP)*SQRT(Z))
        ENDDO
      endif
      end SUBROUTINE ROCKING_CURVE

      end module ROCKING

