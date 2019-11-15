!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable_options.f,v 1.8 2009/12/07 22:32:49 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.8 $
!////     $Date: 2009/12/07 22:32:49 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table for COMPONENTS
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE_OPTIONS
      USE FIELDDEF
      use CLASSDEF,ONLY:getICLS
      use FIELDDATA
      use TRACINGOPT
      use REPORTSOPT
      implicit none


! TOPTOBJ - abstract type with pointers to option types
      type TOPTOBJ
        integer :: ICLS   ! class ID from CLS_IDS (all defined in classes.xml)
        integer :: OCLS
        type(TTROPT),pointer :: P_TRA
        type(TREPOPT),pointer :: P_REP
      end type TOPTOBJ

      private popt_topt_

! overloaded POPTOBJ
      interface POPTOBJ
        module procedure popt_topt_
        module procedure popt_trep_
      end interface

! global instances of TOPTOBJ
      type(TOPTOBJ) :: TRAOBJ   ! tracing
      type(TOPTOBJ) :: REPOBJ   ! reports

      contains

C-------------------------------------------------------------
! get options pointer by ID string
C-------------------------------------------------------------
      type(TOPTOBJ) function getOptByID(IDSTR)
        character(*),intent(in) :: IDSTR
        type(TOPTOBJ) :: RES
        integer :: IDCOMPARE
        call poptobj_clear(RES)
        if (IDCOMPARE(IDSTR,TROPT%ID).eq.0) then
          RES=POPTOBJ(TROPT)
        else if (IDCOMPARE(IDSTR,REPOPT%ID).eq.0) then
          RES=POPTOBJ(REPOPT)
        endif
        getOptByID=RES
      end function getOptByID

!--------------------------------------------------------
! nullifies all pointers in TCOMOBJ
!--------------------------------------------------------
      subroutine poptobj_clear(OBJ)
      type(TOPTOBJ) :: OBJ
        OBJ%ICLS=0
        OBJ%OCLS=0
        NULLIFY(OBJ%P_TRA)
        NULLIFY(OBJ%P_REP)
      end subroutine poptobj_clear

!-------------------------------------------
! overloaded POBJ: define for each object type
!-------------------------------------------
      type(TOPTOBJ) function popt_topt_(OBJ)
        type(TTROPT),target :: OBJ
        type(TOPTOBJ) :: RES
        call poptobj_clear(RES)
        RES%OCLS=OCLS_TRACING
        RES%ICLS=getICLS(RES%OCLS,OCLS_NAMES)
        RES%P_TRA => OBJ
        popt_topt_=RES
      end function popt_topt_
      type(TOPTOBJ) function popt_trep_(OBJ)
        type(TREPOPT),target :: OBJ
        type(TOPTOBJ) :: RES
        call poptobj_clear(RES)
        RES%OCLS=OCLS_REPORTS
        RES%ICLS=getICLS(RES%OCLS,OCLS_NAMES)
        RES%P_REP => OBJ
        popt_trep_=RES
      end function popt_trep_

!---------------------------------------------------------
       SUBROUTINE POPTOBJ_INP(CLS,PNAME,ARG,NARG)
! dispatch calls to *_INP procedures for each object
!---------------------------------------------------------
      TYPE(TOPTOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%OCLS)
        case(OCLS_TRACING)
          call TROPT_INP(CLS%P_TRA,PNAME,ARG,NARG)
        case(OCLS_REPORTS)
          call REPOPT_INP(CLS%P_REP,PNAME,ARG,NARG)
        end select
      end SUBROUTINE POPTOBJ_INP

!---------------------------------------------------------
       SUBROUTINE POPTOBJ_OUT(CLS,PNAME,ARG,NARG)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TOPTOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        select case(CLS%OCLS)
        case(OCLS_TRACING)
          call TROPT_OUT(CLS%P_TRA,PNAME,ARG,NARG)
        case(OCLS_REPORTS)
          call REPOPT_OUT(CLS%P_REP,PNAME,ARG,NARG)
        end select
      end SUBROUTINE POPTOBJ_OUT

!---------------------------------------------------------
       SUBROUTINE POPTOBJ_DEFAULT(CLS)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TOPTOBJ) :: CLS
        select case(CLS%OCLS)
        case(OCLS_TRACING)
          call TROPT_DEFAULT(CLS%P_TRA)
        case(OCLS_REPORTS)
          call REPOPT_DEFAULT(CLS%P_REP)
        end select
      end SUBROUTINE POPTOBJ_DEFAULT


C-------------------------------------------------------------
      SUBROUTINE GetOptAttrib(CLS,CLSID,OBJID,OBJNAME)
! Get attributes (class ID, object ID, name) by object pointer
! INPUT
!   CLS    .. TOPTOBJ reference
! OUTPUT
!   CLSID     .. class ID
!   OBJID     .. object ID
!   OBJNAME   .. descriptive name
C-------------------------------------------------------------
      TYPE(TOPTOBJ),intent(in) :: CLS
      character(*),intent(out) :: CLSID,OBJID,OBJNAME
      integer :: IS,IL
      CLSID=' '
      OBJID=' '
      OBJNAME=' '
      if ((CLS%ICLS.le.0).or.(CLS%OCLS.le.0)) return
      call FINDSTRPAR(OCLS_NAMES,':',CLS%OCLS,IS,IL)
      if (IL.le.0) return
      CLSID=OCLS_NAMES(IS:IS+IL-1)
      select case(CLS%OCLS)
      case(OCLS_TRACING)
        if (associated(CLS%P_TRA)) then
          OBJID=trim(CLS%P_TRA%ID)
          OBJNAME=trim(CLS%P_TRA%NAME)
        endif
      case(OCLS_REPORTS)
        if (associated(CLS%P_REP)) then
          OBJID=trim(CLS%P_REP%ID)
          OBJNAME=trim(CLS%P_REP%NAME)
        endif
      end select
      end SUBROUTINE GetOptAttrib

C-------------------------------------------------------------
      SUBROUTINE SetOptAttrib(CLS,OBJID,OBJNAME)
! Set class attributes (object ID, name) to the object given as pointer
! INPUT
!   OBJID     .. object ID
!   OBJNAME   .. descriptive name
! OUTPUT
!   CLS    .. TOPTOBJ reference, for which we change ID,NAME and CLASS
C-------------------------------------------------------------
      TYPE(TOPTOBJ),intent(inout) :: CLS
      character(*),intent(in) :: OBJID,OBJNAME
      select case(CLS%OCLS)
      case(OCLS_TRACING)
        CLS%P_TRA%ID=trim(OBJID)
        CLS%P_TRA%NAME=trim(OBJNAME)
      case(OCLS_REPORTS)
        CLS%P_REP%ID=trim(OBJID)
        CLS%P_REP%NAME=trim(OBJNAME)
      end select
      end SUBROUTINE SetOptAttrib

      end MODULE VCTABLE_OPTIONS

