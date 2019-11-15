!//////////////////////////////////////////////////////////////////////
!////  $Id: vctable_commands.f,v 1.9 2018/10/29 17:50:48 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2018/10/29 17:50:48 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Virtual classes table for COMMANDS
!////  - serves to emulation of polymorphism in Fortran
!//////////////////////////////////////////////////////////////////////
      MODULE VCTABLE_COMMANDS
      use FIELDDATA
      USE COMMANDS
      implicit none


! TCMDOBJ - abstract type with pointer to a commmand
      type TCMDOBJ
        integer :: ICLS   ! class ID from CLS_IDS (all defined in classes.xml)
        integer :: ICMD   ! class ID from CCLS_NAMES
        type(TCOMMAND),pointer :: P_CMD
      end type TCMDOBJ

      contains


!--------------------------------------------------------
! clear data in TCOMOBJ, do not dispose allocated fielads
!--------------------------------------------------------
      subroutine PCMDOBJ_CLEAR(OBJ)
      type(TCMDOBJ) :: OBJ
        OBJ%ICLS=0
        OBJ%ICMD=0
        NULLIFY(OBJ%P_CMD)
      end subroutine PCMDOBJ_CLEAR

!--------------------------------------------------------
! nullifies all pointers in TCOMOBJ
!--------------------------------------------------------
      subroutine PCMDOBJ_DISPOSE(OBJ)
      type(TCMDOBJ) :: OBJ
        if (associated(OBJ%P_CMD)) call DISPOSE_COMMAND(OBJ%P_CMD)
        call PCMDOBJ_CLEAR(OBJ)
      end subroutine PCMDOBJ_DISPOSE

!-------------------------------------------
! overloaded PSAMOBJ: define for each object type
!-------------------------------------------
      type(TCMDOBJ) function PCMDOBJ(OBJ)
        type(TCOMMAND),target :: OBJ
        type(TCMDOBJ) :: RES
        call PCMDOBJ_CLEAR(RES)
        RES%ICMD=OBJ%ICMD
        RES%ICLS=OBJ%ICLS
        RES%P_CMD => OBJ
        PCMDOBJ=RES
      end function PCMDOBJ

C---------------------------------------------------------
      type(TCMDOBJ) function getCmdObjByID(IDSTR)
! return command pointer for given ID string
C---------------------------------------------------------
      character(*),intent(in) :: IDSTR
      integer :: ic
      type(TCMDOBJ) :: RES
      call PCMDOBJ_CLEAR(RES)
      ic=getCmdIndex(IDSTR)
      if (ic.gt.0) then
        RES=PCMDOBJ(myCOMMANDS(ic))
  !      write(*,*) 'getCmdObjByID ',trim(IDSTR),' ic=',ic,' ICLS=',RES%ICLS,' ICMD=',RES%ICMD
      endif
      getCmdObjByID=RES
      end function getCmdObjByID

C---------------------------------------------------------
      type(TCMDOBJ) function getCmdObj(ic)
! return command pointer for given ID string
C---------------------------------------------------------
      integer,intent(in) :: ic
      type(TCMDOBJ) :: RES
      call PCMDOBJ_CLEAR(RES)
      if ((ic.gt.0).and.(ic.le.COMMANDS_NUM)) then
        RES=PCMDOBJ(myCOMMANDS(ic))
      endif
      getCmdObj=RES
      end function getCmdObj

!---------------------------------------------------------
       SUBROUTINE PCMDOBJ_INP(CLS,PNAME,ARG,NARG)
! dispatch calls to *_INP procedures for each object
!---------------------------------------------------------
      TYPE(TCMDOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
  !      write(*,*) 'PCMDOBJ_INP_R [',trim(CLS%P_CMD%ID),'] ',trim(PNAME)
        call COMMANDS_INP(CLS%P_CMD,PNAME,ARG,NARG)
      end SUBROUTINE PCMDOBJ_INP

!---------------------------------------------------------
       SUBROUTINE PCMDOBJ_OUT(CLS,PNAME,ARG,NARG)
! dispatch calls to *_OUT procedures for each object
!---------------------------------------------------------
      TYPE(TCMDOBJ) :: CLS
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
        call COMMANDS_OUT(CLS%P_CMD,PNAME,ARG,NARG)
      end SUBROUTINE PCMDOBJ_OUT

      end MODULE VCTABLE_COMMANDS
