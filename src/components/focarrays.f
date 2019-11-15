!//////////////////////////////////////////////////////////////////////
!////  $Id: focarrays.f,v 1.9 2019/08/15 15:02:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes FOCARRAY - focussing array for crystals
!////  Used as subclass for crystal components
!////  Implements data structure and basic procedures needed to set focusing parameters
!////////////////////////////////////////////////////////////////////////
      MODULE FOCARRAYS
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      implicit none

      save
      private

! FOCARRAY
      TYPE TFOCARRAY
        integer :: NSEG
        logical :: AUTO
        logical :: BENT
        logical :: STACK
        real(kind(1.D0)) :: GAP
        real(kind(1.D0)) :: RHO
        real(kind(1.D0)) :: FOC1
        real(kind(1.D0)) :: FOC2
      END TYPE TFOCARRAY

      integer, parameter :: FOCARRAYS_DIM=2
! pointer type for FOCARRAY
      type PFOCARRAY
        integer :: IDX  ! instance index
        TYPE(TFOCARRAY),pointer :: X
      end type PFOCARRAY
! instances of FOCARRAY. AFOCARRAYS(0) is always unallocated
      integer :: FOCARRAYS_NC
      type(PFOCARRAY) :: AFOCARRAYS(1:FOCARRAYS_DIM)

      public TFOCARRAY
      public FOCARRAY_INP,FOCARRAY_OUT,FOCARRAY_DEFAULT

      contains

!---------------------------------------------------------
      SUBROUTINE FOCARRAY_DEFAULT(OBJ)
      TYPE(TFOCARRAY),intent(out) :: OBJ
        OBJ%NSEG=1
        OBJ%AUTO=.false.
        OBJ%BENT=.false.
        OBJ%STACK=.false.
        OBJ%RHO=0.D0
        OBJ%GAP=0.D0
        OBJ%FOC1=2.D0
        OBJ%FOC2=2.D0
      end SUBROUTINE FOCARRAY_DEFAULT

!---------------------------------------------------------
      SUBROUTINE FOCARRAY_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFOCARRAY):: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call FOCARRAY_INP_R(OBJ,trim(ARG%ID),NUM,LR)
      end select
      NARG=LR
      END SUBROUTINE FOCARRAY_INP

!---------------------------------------------------------
      SUBROUTINE FOCARRAY_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFOCARRAY) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,NA
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FOCARRAY_OUT_R(OBJ,trim(ARG%ID),NUM,NA)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
      end select
      NARG=LR
      END SUBROUTINE FOCARRAY_OUT

!---------------------------------------------------------
      SUBROUTINE FOCARRAY_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFOCARRAY),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('FOC1');   OBJ%FOC1=ARG(1)
        CASE('FOC2');   OBJ%FOC2=ARG(1)
        CASE('GAP');    OBJ%GAP=ARG(1)
        CASE('RHO');    OBJ%RHO=ARG(1)*1.D-3
        CASE('NSEG');   OBJ%NSEG=max(1,NINT(ARG(1)))
        CASE('AUTO');   OBJ%AUTO=(ARG(1).eq.1.D0)
        CASE('BENT');   OBJ%BENT=(ARG(1).eq.1.D0)
        CASE('STACK');  OBJ%STACK=(ARG(1).eq.1.D0)
      END SELECT
      NARG=LR
      END SUBROUTINE FOCARRAY_INP_R

!---------------------------------------------------------
      SUBROUTINE FOCARRAY_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TFOCARRAY),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('FOC1');   ARG(1)=OBJ%FOC1
        CASE('FOC2');   ARG(1)=OBJ%FOC2
        CASE('GAP');    ARG(1)=OBJ%GAP
        CASE('RHO');    ARG(1)=OBJ%RHO*1.D3
        CASE('NSEG');   ARG(1)=OBJ%NSEG
        CASE('AUTO');   ARG(1)=LOG2INT(OBJ%AUTO)
        CASE('BENT');   ARG(1)=LOG2INT(OBJ%BENT)
        CASE('STACK');  ARG(1)=LOG2INT(OBJ%STACK)
      END SELECT
      NARG=LR
      END SUBROUTINE FOCARRAY_OUT_R

      end MODULE FOCARRAYS