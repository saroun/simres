!//////////////////////////////////////////////////////////////////////
!////  $Id: beamline.f,v 1.1 2007/10/08 08:35:37 saroun Exp $
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
!////////////////////////////////////////////////////////////////////////////
!////
!////  Envelop subroutines for ray-tracing through the beamline.
!////  Assumes that actual beamline has already been defined in COMPONENTS module.
!////////////////////////////////////////////////////////////////////////////
      MODULE BEAMLINE
      USE COMPONENTS
      implicit none

C// ID of all defined beamlines
      integer,parameter :: BML_ROCKING=1
      integer,parameter :: BML_TAS1=2

      contains

!---------------------------------------------------------------------
! increment counters of started events
!---------------------------------------------------------------------
C--------------------------------------------------------
      SUBROUTINE BEAMLINE_INI
C General initialization procedure
C--------------------------------------------------------
      implicit NONE
      INCLUDE 'const.inc'
      INCLUDE 'inout.inc'
      INCLUDE 'randvars.inc'

      type(GENERATOR),intent(in) :: GENER
      type(CRYSTAL),intent(in) :: CRYST

      REAL(KIND(1.D0)) :: EIVALS(CRND),EIVECT(CRND,CRND)
      integer :: N,I,J


      GEN=GENER
      CR=CRYST
      CR%FRAME.DIST=100.D0 ! standard source distance = 10 cm

! default rand. var. pool settings
      DO I=1,CRND
        RNDLIST.MEAN(I)=0.D0
        RNDLIST.POOL(I)=1.1D0
        RNDLIST.ACTIVE(I)=1
        RNDLIST.LIMITS(I)=1.D0
      enddo

C// initialize transformation matrix to a diagonal
      TMAT=0.D0
      DO i=1,CRND
        TMEAN(i)=0.D0
        TMAT(i,i)=1.D0
      ENDDO

C// Initialize generator
      CALL GENERATOR_INIT(TRACING_GEN)

C// Use RESMAT to optimize sampling
      call RESMAT_INIT
      call RESMAT_ADDGENERATOR(TRACINGGEN) ! does nothing unless limits are applied
    ! add all components on the beamline
      do i=1,BEAMLINE_NC
        call RESMAT_ADDCOMPONENT(BEAMLINE(i,1),BEAMLINE(i,2))
      enddo
      call RESMAT_EVAL
      call RESMAT_EIGEN(0,EIVALS,EIVECT,CRND,N)

C// set number of random variables
      RNDLIST.DIM=N

! get new events transformation matrix and limits
  ! calculate new TMAT and limits
      DO I=1,RNDLIST.DIM
      IF (RNDLIST.ACTIVE(I).GT.0) THEN
        RNDLIST.LIMITS(I)=2*SQRT(6*ABS(EIVALS(I)))
        DO J=1,RNDLIST.DIM
          IF (RNDLIST.ACTIVE(J).GT.0) THEN
             TMAT(I,J)=EIVECT(J,I)
          ENDIF
        ENDDO
      ENDIF
      ENDDO

    ! report transformation matrix
      if( idbg.gt.0) then
        call XML_RSXDUMP(sout,' ',1)
        call XML_VALUE(sout,'THETA','float',' ',0,CR%FRAME%GON(1)/deg,' ')
        call XML_MATRIX(sout,N,N,'TMAT','x|y|kx|ky|kz|c1|c2','x|y|kx|ky|kz|c1|c2',TMAT(1,1),CRND)
        call XML_ARRAY(sout,N,'LIMIT','x|y|kx|ky|kz|c1|c2',RNDLIST.LIMITS(1))
        call XML_RSXDUMP(sout,' ',0)
      endif

      END SUBROUTINE ROCKING_INI


      end MODULE BEAMLINES
