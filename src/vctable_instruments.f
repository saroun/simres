!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable_instruments.f,v 1.10 2009/12/07 22:32:49 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.10 $
!////     $Date: 2009/12/07 22:32:49 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table for INSTRUMENTS
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE_INSTRUMENTS
      use CLASSES
      use CLASSDEF,ONLY:GetICLS
      use SPECTROMETER
      implicit none

! TINSOBJ - abstract type with pointer to a component
      type TINSOBJ
        integer :: ICLS   ! class ID from CLS_IDS (all defined in classes.xml)
        integer :: DCLS   ! class ID from DCLS_NAMES (device interfaces)
        type(TSPEC),pointer :: P_SPEC
    ! to be added ...
    !   type(TPOWDER),pointer :: P_POWDER
    !   type(TTAS),pointer :: P_TAS
      end type TINSOBJ

! overloaded PSPECOBJ
      interface PINSOBJ
        module procedure pins_spec_
      !  module procedure pins_powder_
      !  module procedure pins_tas_
      end interface

! global instance of TINSTOBJ
      type(TINSOBJ) :: INSOBJ

      contains

!--------------------------------------------------------
! nullifies all pointers in TOBJ
!--------------------------------------------------------
      subroutine pinsobj_clear(OBJ)
      type(TINSOBJ) :: OBJ
        OBJ%ICLS=0
        OBJ%DCLS=0
        NULLIFY(OBJ%P_SPEC)
    !    NULLIFY(OBJ%P_POWDER)
    !    NULLIFY(OBJ%P_TAS)
      end subroutine pinsobj_clear

!-------------------------------------------
! overloaded PINSOBJ: define for each object type
!-------------------------------------------
      type(TINSOBJ) function pins_spec_(OBJ)
        type(TSPEC),target :: OBJ
        type(TINSOBJ) :: RES
        call PINSOBJ_CLEAR(RES)
        RES%DCLS=DCLS_SPEC
        RES%ICLS=getICLS(RES%DCLS,DCLS_NAMES)
        RES%P_SPEC => OBJ
        pins_spec_=RES
      end function pins_spec_

!---------------------------------------------------------
       SUBROUTINE PINSOBJ_ASSIGN(OBJ,SRC)
! assign object data from SRC to OBJ
!---------------------------------------------------------
      TYPE(TINSOBJ) :: OBJ,SRC
      if (associated(SRC%P_SPEC)) then
        call pinsobj_clear(OBJ)
        select case(SRC%DCLS)
        case(DCLS_SPEC)
          OBJ=PINSOBJ(SPEC)
          OBJ%P_SPEC=SRC%P_SPEC
        case(DCLS_PWDIF)
!          OBJ=PINSOBJ(PWDIF)
!          OBJ%P_PWDIF=SRC%P_PWDIF
!          OBJ%P_SPEC => OBJ%P_POWDER%SPEC
        case(DCLS_TAS)
!          OBJ=PINSOBJ(TAS)
!          OBJ%P_TAS=SRC%P_TAS
!          OBJ%P_SPEC => OBJ%P_TAS%SPEC
        case default
          return
        end select
        OBJ%ICLS=SRC%ICLS
        OBJ%DCLS=SRC%DCLS
      endif
      end SUBROUTINE PINSOBJ_ASSIGN

!---------------------------------------------------------
       SUBROUTINE PINSOBJ_INP(CLS,PNAME,ARG,NARG)
! dispatch calls to *_INP procedures for each object
!---------------------------------------------------------
      TYPE(TINSOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%DCLS)
        case(DCLS_SPEC)
          call SPEC_INP(CLS%P_SPEC,PNAME,ARG,NARG)
        case(DCLS_PWDIF)
  !        call POWDER_INP(CLS%P_POWDER,PNAME,ARG,NARG)
        case(DCLS_TAS)
  !        call TAS_INP(CLS%P_TAS,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PINSOBJ_INP

!---------------------------------------------------------
       SUBROUTINE PINSOBJ_OUT(CLS,PNAME,ARG,NARG)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TINSOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%DCLS)
        case(DCLS_SPEC)
          call SPEC_OUT(CLS%P_SPEC,PNAME,ARG,NARG)
        case(DCLS_PWDIF)
  !        call POWDER_OUT(CLS%P_POWDER,PNAME,ARG,NARG)
        case(DCLS_TAS)
  !        call TAS_OUT(CLS%P_TAS,PNAME,ARG,NARG)
        end select
      end SUBROUTINE PINSOBJ_OUT

!---------------------------------------------------------
       SUBROUTINE PINSOBJ_DEFAULT(CLS)
! dispatch calls to *_DEFAULT procedures for each object
!---------------------------------------------------------
      TYPE(TINSOBJ) :: CLS
        select case(CLS%DCLS)
        case(DCLS_SPEC)
          call SPEC_DEFAULT(CLS%P_SPEC)
        case(DCLS_PWDIF)
  !        call PWD_DEFAULT(CLS%P_POWDER)
        case(DCLS_TAS)
  !        call TAS_DEFAULT(CLS%P_TAS)
        end select
      end SUBROUTINE PINSOBJ_DEFAULT


!---------------------------------------------------------
       SUBROUTINE PINSOBJ_INIT(CLS)
! dispatch calls to *_DEFAULT procedures for each object
!---------------------------------------------------------
      TYPE(TINSOBJ) :: CLS
        select case(CLS%DCLS)
        case(DCLS_SPEC)
          call SPEC_INIT(CLS%P_SPEC)
        case(DCLS_PWDIF)
  !        call PWD_INIT(CLS%P_POWDER)
        case(DCLS_TAS)
  !        call TAS_INIT(CLS%P_TAS)
        end select


      end SUBROUTINE PINSOBJ_INIT

      end MODULE VCTABLE_INSTRUMENTS