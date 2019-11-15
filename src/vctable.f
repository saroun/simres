!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable.f,v 1.34 2019/01/10 20:39:00 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.34 $
!////     $Date: 2019/01/10 20:39:00 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table
!////  - defines virtual classes - common ancestor of all classes
!////  - implements common input/output procedures
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE
      use CLASSES
    !  use FIELDDEF
      use CLASSDEF
      use COMMANDS
      use VCTABLE_COMPONENTS
      use VCTABLE_SAMPLES
      use VCTABLE_INSTRUMENTS
      use VCTABLE_OPTIONS
      use VCTABLE_COMMANDS

      implicit none
      save
      private

! TCLASS - base type for all abstract classes
      type TCLASS
        integer :: ICLS   ! class ID from CLS_IDS (all defined in classes.xml)
        integer :: TCLS ! TCLS_* value
        type(TCLSDEF),pointer :: DEF
        character(LEN_ID) :: ID
        character(LEN_NAME) :: NAME
        type(TCOMOBJ) :: P_COM
        type(TSAMOBJ) :: P_SAM
        type(TINSOBJ) :: P_INS
        type(TOPTOBJ) :: P_OPT
        type(TCMDOBJ) :: P_CMD
        integer :: P_IGP
      end type TCLASS

! overloaded PCLASS
      interface PCLASS
        module procedure pclass_com_
        module procedure pclass_sam_
        module procedure pclass_ins_
        module procedure pclass_opt_
        module procedure pclass_cmd_
      end interface

      private pclass_com_
      private pclass_sam_
      private pclass_ins_
      private pclass_opt_
      private pclass_cmd_


! local char array for full component namespaces, shared by _INP procedures
      character(2*DIM_ID),private :: SGLOB

      public TCLASS,PCLASS,PCLASS_CLEAR,PCLASS_DISPOSE
      public GetClassAttrib,SetClassAttrib
      public PCLASS_DEFAULT,PCLASS_INIT,PCLASS_INP,PCLASS_OUT

      contains

!--------------------------------------------------------
! nullifies all pointers in TOBJ
!--------------------------------------------------------
      subroutine PCLASS_CLEAR(OBJ)
      type(TCLASS) :: OBJ
        OBJ%TCLS=0
        OBJ%ICLS=0
        OBJ%ID=' '
        OBJ%NAME=' '
        nullify(OBJ%DEF)
        call PCOM_CLEAR(OBJ%P_COM)
        call PSAMOBJ_CLEAR(OBJ%P_SAM)
        call PINSOBJ_CLEAR(OBJ%P_INS)
        call POPTOBJ_CLEAR(OBJ%P_OPT)
        call PCMDOBJ_CLEAR(OBJ%P_CMD)
        OBJ%P_IGP=0
      end subroutine PCLASS_CLEAR

!--------------------------------------------------------
! disposes memory and nullifies all pointers in TOBJ
!--------------------------------------------------------
      subroutine PCLASS_DISPOSE(OBJ)
      type(TCLASS) :: OBJ
        select case(OBJ%TCLS)
        case(TCLS_COM)
          call PCOM_DISPOSE(OBJ%P_COM)
      !  case(TCLS_SAM)
      !  case(TCLS_INS)
        case(TCLS_CMD)
          call PCMDOBJ_DISPOSE(OBJ%P_CMD)
        end select
        call PCLASS_CLEAR(OBJ)
      end subroutine PCLASS_DISPOSE

!-------------------------------------------
! overloaded PCLASS:
!-------------------------------------------

      type(TCLASS) function pclass_com_(OBJ)
        type(TCOMOBJ) :: OBJ
        type(TCLASS) :: RES
        type(PCLSDEF) :: CDEF
        call pclass_clear(RES)
        call GetClassDef(OBJ%ICLS,CDEF)
        if (associated(CDEF%C).and.(OBJ%ICLS.gt.0)) then
          RES%DEF=>CDEF%C
          RES%ICLS=OBJ%ICLS
          RES%TCLS=TCLS_COM
          RES%P_COM = OBJ
          RES%ID=OBJ%ID
          RES%NAME='not needed' !OBJ%P_FRAME%NAME
        endif
        pclass_com_=RES
      end function pclass_com_
      type(TCLASS) function pclass_ins_(OBJ)
        type(TINSOBJ) :: OBJ
        type(TCLASS) :: RES
        type(PCLSDEF) :: CDEF
        call pclass_clear(RES)
        call GetClassDef(OBJ%ICLS,CDEF)
        if (associated(CDEF%C).and.(OBJ%ICLS.gt.0)) then
          RES%DEF=>CDEF%C
          RES%ICLS=OBJ%ICLS
          RES%TCLS=TCLS_INS
          RES%P_INS = OBJ
          RES%ID=OBJ%P_SPEC%FRAME%ID
          RES%NAME=OBJ%P_SPEC%FRAME%NAME
        endif
        pclass_ins_=RES
      end function pclass_ins_
      type(TCLASS) function pclass_sam_(OBJ)
        type(TSAMOBJ) :: OBJ
        type(TCLASS) :: RES
        type(PCLSDEF) :: CDEF
        call pclass_clear(RES)
        call GetClassDef(OBJ%ICLS,CDEF)
        if (associated(CDEF%C).and.(OBJ%ICLS.gt.0)) then
          RES%DEF=>CDEF%C
          RES%ICLS=OBJ%ICLS
          RES%TCLS=TCLS_SAM
          RES%P_SAM = OBJ
          RES%ID=OBJ%P_FRAME%ID
          RES%NAME=OBJ%P_FRAME%NAME
        endif
        pclass_sam_=RES
      end function pclass_sam_
      type(TCLASS) function pclass_opt_(OBJ)
        type(TOPTOBJ) :: OBJ
        type(TCLASS) :: RES
        type(PCLSDEF) :: CDEF
        character(LEN_ID) :: CLSID,OBJID
        character(LEN_NAME) :: OBJNAME
        call pclass_clear(RES)
        call GetClassDef(OBJ%ICLS,CDEF)
        if (associated(CDEF%C).and.(OBJ%ICLS.gt.0)) then
          RES%DEF=>CDEF%C
          RES%ICLS=OBJ%ICLS
          RES%TCLS=TCLS_OPT
          RES%P_OPT = OBJ
          call GetOptAttrib(OBJ,CLSID,OBJID,OBJNAME)
          RES%ID=OBJID
          RES%NAME=OBJNAME
        endif
        pclass_opt_=RES
      end function pclass_opt_
      type(TCLASS) function pclass_cmd_(OBJ)
        type(TCMDOBJ) :: OBJ
        type(TCLASS) :: RES
        type(PCLSDEF) :: CDEF
        call pclass_clear(RES)
        call GetClassDef(OBJ%ICLS,CDEF)
        if (associated(CDEF%C).and.(OBJ%ICLS.gt.0)) then
          RES%DEF=>CDEF%C
          RES%ICLS=OBJ%ICLS
          RES%TCLS=TCLS_CMD
          RES%P_CMD = OBJ
          RES%ID=OBJ%P_CMD%ID
          RES%NAME=OBJ%P_CMD%NAME
        endif
        pclass_cmd_=RES
      end function pclass_cmd_



!******************** DISPATCH CLASS INPUT / OUTPUT ******************

C---------------------------------------------------------
      SUBROUTINE PCLASS_INP(CLS,PNAME,ARG,NARG)
C input data from ARG to CLS for parameter name PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... TYPE(TFIELD) parameter data
C OUTPUT
C    CLS     .... TOBJ (abstract type)
C    NARG    .... number of values read from ARG
C---------------------------------------------------------
      TYPE(TCLASS) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
    ! OBJ must have class definition associated
      NARG=0
      if (associated(ARG%DEF)) THEN
        select case (CLS%TCLS)
        case(TCLS_COM)
          call PCOMOBJ_INP(CLS%P_COM,PNAME,ARG,NARG)
        case(TCLS_SAM)
          call PSAMOBJ_INP(CLS%P_SAM,PNAME,ARG,NARG)
        case(TCLS_INS)
          call PINSOBJ_INP(CLS%P_INS,PNAME,ARG,NARG)
        case(TCLS_OPT)
          call POPTOBJ_INP(CLS%P_OPT,PNAME,ARG,NARG)
        case(TCLS_CMD)
          call PCMDOBJ_INP(CLS%P_CMD,PNAME,ARG,NARG)
        end select
      endif
      END SUBROUTINE PCLASS_INP

C---------------------------------------------------------
      SUBROUTINE PCLASS_OUT(CLS,PNAME,ARG,NARG)
C input data from ARG to CLS for parameter name PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... TYPE(TFIELD) parameter data
C OUTPUT
C    CLS     .... TOBJ (abstract type)
C    NARG    .... number of values read from ARG
C---------------------------------------------------------
      TYPE(TCLASS) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
    ! ARG must have field definition associated
      NARG=0
      if (associated(ARG%DEF)) THEN
        select case (CLS%TCLS)
        case(TCLS_COM)
          call PCOMOBJ_OUT(CLS%P_COM,PNAME,ARG,NARG)
        case(TCLS_SAM)
          call PSAMOBJ_OUT(CLS%P_SAM,PNAME,ARG,NARG)
        case(TCLS_INS)
          call PINSOBJ_OUT(CLS%P_INS,PNAME,ARG,NARG)
        case(TCLS_OPT)
          call POPTOBJ_OUT(CLS%P_OPT,PNAME,ARG,NARG)
        case(TCLS_CMD)
          call PCMDOBJ_OUT(CLS%P_CMD,PNAME,ARG,NARG)
        end select
      endif
      END SUBROUTINE PCLASS_OUT

C-------------------------------------------------------------
      SUBROUTINE GetClassAttrib(CLS,CLSID,OBJID,OBJNAME)
! Get class attributes (class ID, object ID, name) by object pointer
! INPUT
!   OBJ    .. TOBJ reference
! OUTPUT
!   CLSID     .. class ID
!   OBJID     .. object ID
!   OBJNAME   .. descriptive name
C-------------------------------------------------------------
      TYPE(TCLASS),intent(in) :: CLS
      character(*),intent(out) :: CLSID,OBJID,OBJNAME
      integer :: IS,IL
      TYPE(PFRAME) :: F
      CLSID=' '
      OBJID=' '
      OBJNAME=' '
      select case(CLS%TCLS)
      case(TCLS_COM)
        if (CLS%P_COM%ICLS.le.0) return
        call PCOMOBJ_GET(CLS%P_COM,F)
        if (associated(F%X)) then
          call FINDSTRPAR(CCLS_NAMES,':',F%X%CLASS,IS,IL)
          if (IL.gt.0) then
            CLSID=CCLS_NAMES(IS:IS+IL-1)
            OBJID=trim(F%X%ID)
            OBJNAME=trim(F%X%NAME)
          endif
        endif
      case(TCLS_SAM)
        if (CLS%P_SAM%ICLS.le.0) return
        if (associated(CLS%P_SAM%P_FRAME)) then
          call FINDSTRPAR(SCLS_NAMES,':',CLS%P_SAM%P_FRAME%CLASS,IS,IL)
          if (IL.gt.0) then
            CLSID=SCLS_NAMES(IS:IS+IL-1)
            OBJID=trim(CLS%P_SAM%P_FRAME%ID)
            OBJNAME=trim(CLS%P_SAM%P_FRAME%NAME)
          endif
        endif
      case(TCLS_INS)
        if (CLS%P_INS%ICLS.le.0) return
        if (associated(CLS%P_INS%P_SPEC)) then
          call FINDSTRPAR(DCLS_NAMES,':',CLS%P_INS%P_SPEC%CLASS,IS,IL)
          if (IL.gt.0) then
            CLSID=DCLS_NAMES(IS:IS+IL-1)
            OBJID=trim(CLS%P_INS%P_SPEC%FRAME%ID)
            OBJNAME=trim(CLS%P_INS%P_SPEC%FRAME%NAME)
          endif
        endif
      case(TCLS_OPT)
        call GetOptAttrib(CLS%P_OPT,CLSID,OBJID,OBJNAME)
      case(TCLS_CMD)
        if (CLS%P_CMD%ICLS.le.0) return
        if (associated(CLS%P_CMD%P_CMD)) then
            CLSID=CLS%ID
            OBJID=CLS%ID
            OBJNAME=trim(CLS%P_CMD%P_CMD%NAME)
        endif
      end select
      end SUBROUTINE GetClassAttrib

C-------------------------------------------------------------
      SUBROUTINE SetClassAttrib(CLS,OBJID,OBJNAME)
! Set class attributes (object ID, name) to the object given as pointer
! INPUT
!   OBJID     .. object ID
!   OBJNAME   .. descriptive name
! OUTPUT
!   CLS    .. TCLASS reference, for which we change ID,NAME and CLASS
C-------------------------------------------------------------
      TYPE(TCLASS),intent(inout) :: CLS
      character(*),intent(in) :: OBJID,OBJNAME
        select case(CLS%TCLS)
        case(TCLS_COM)
        !  CLS%P_COM%P_FRAME%ID=trim(OBJID)
        !  CLS%P_COM%P_FRAME%NAME=trim(OBJNAME)
        case(TCLS_SAM)
          CLS%P_SAM%P_FRAME%ID=trim(OBJID)
          CLS%P_SAM%P_FRAME%NAME=trim(OBJNAME)
        case(TCLS_INS)
          CLS%P_INS%P_SPEC%FRAME%ID=trim(OBJID)
          CLS%P_INS%P_SPEC%FRAME%NAME=trim(OBJNAME)
        case(TCLS_OPT)
          call SetOptAttrib(CLS%P_OPT,OBJID,OBJNAME)
        end select
        CLS%ID=trim(OBJID)
        CLS%NAME=trim(OBJNAME)
      end SUBROUTINE SetClassAttrib

!---------------------------------------------------------
       SUBROUTINE PCLASS_DEFAULT(CLS)
! dispatch calls to *_DEFAULT procedures for each object
!---------------------------------------------------------
      TYPE(TCLASS) :: CLS
        select case(CLS%TCLS)
        case(TCLS_COM)
          call PCOMOBJ_DEFAULT(CLS%P_COM)
        case(TCLS_SAM)
          call PSAMOBJ_DEFAULT(CLS%P_SAM)
        case(TCLS_INS)
          call PINSOBJ_DEFAULT(CLS%P_INS)
        case(TCLS_OPT)
          call POPTOBJ_DEFAULT(CLS%P_OPT)
      !    write(*,*) 'PCLASS_DEFAULT ',TCLS_OPT
        end select
      end SUBROUTINE PCLASS_DEFAULT


!---------------------------------------------------------
       SUBROUTINE PCLASS_INIT(CLS)
! dispatch calls to *_INIT procedures for each object
! not used (2/6/2011)
!---------------------------------------------------------
      TYPE(TCLASS) :: CLS
        select case(CLS%TCLS)
        case(TCLS_COM)
          call PCOMOBJ_INIT(CLS%P_COM)
        case(TCLS_SAM)
          call PSAMOBJ_INIT(CLS%P_SAM)
        case(TCLS_INS)
          call PINSOBJ_INIT(CLS%P_INS)
        end select
      end SUBROUTINE PCLASS_INIT


!---------------------------------------------------------

      end MODULE VCTABLE