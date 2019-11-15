!//////////////////////////////////////////////////////////////////////
!////  $Id: classdata.f,v 1.2 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.2 $
!////     $Date: 2019/08/15 15:02:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Defines general class data
!////
!//////////////////////////////////////////////////////////////////////
      MODULE CLASSDATA
      use XMLINFO
      USE XMLUTILS
      use FIELDDEF
      use CLASSDEF
      USE FIELDDATA
      implicit none

      private
      save

! CLSDATA definition
      integer,parameter :: CLSD_MAXARG=16
      integer,parameter :: CLSDATA_DIM=32

! common type for all CLSDATA
      type TCLSDATA
        INTEGER :: ICLS=0 ! class index (see global definitions in CLASSES module)
      !  INTEGER :: ICMD=0 ! CLSDATA index (pointer to ACLSDATA array in this module)
        character(LEN_ID)     :: ID=''       ! ID string
        type(TCLSDEF),pointer :: CDEF
        integer :: NARG=0 ! number of numerical arguments
        TYPE(TFIELD) :: FIELDS(CLSD_MAXARG)
      end TYPE TCLSDATA

      TYPE PCLSDATA
        integer :: INST
        type(TCLSDATA),pointer :: X
      end TYPE PCLSDATA

  ! number of CLSDATA defined
      integer :: CLSDATA_NC
  ! instances of command classes
      type(PCLSDATA) :: ACLSDATA(CLSDATA_DIM)

      public CLSDATA_CREATE,CLSDATA_DISPOSE,CLSDATA_DISPOSE_ALL
      public CLSDATA_INP,CLSDATA_OUT

      contains

!-------------------------------------------------------------
      logical function CLSDATA_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      CLSDATA_isValid= ((INST.gt.0).and.(INST.le.CLSDATA_DIM).and.associated(ACLSDATA(INST)%X))
      end function CLSDATA_isValid

!-------------------------------------------------------------------------
! Creator for CLSDATA, return instance number
! Memory is allocated on the first free element of ACLSDATAS array
!-------------------------------------------------------------------------
      integer function CLSDATA_CREATE(ID,CDEF)
      character(*) :: ID
      type(TCLSDEF),pointer :: CDEF
      integer :: ierr,i
      character(32) :: S
      logical :: error
        ierr=1
        i=1
        do while ((i.lt.CLSDATA_DIM).and.(ACLSDATA(i)%INST.gt.0))
          i=i+1
        enddo
        if (ACLSDATA(i)%INST.le.0) then
          allocate(ACLSDATA(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            CLSDATA_NC=CLSDATA_NC+1
            ACLSDATA(i)%INST=i
            ACLSDATA(i)%X%ID=trim(ID)
            ACLSDATA(i)%X%ICLS=CDEF%IDX
            ACLSDATA(i)%X%CDEF=>CDEF
            call CLSDATA_ALLOC_FIELDS(ACLSDATA(i)%X,error)
            if (error) ierr=2
          endif
        endif
        if (ierr.eq.0) then
          CLSDATA_CREATE=i
        else
          call INT2STR(CLSDATA_DIM,S)
          call MSG_ERROR('CLSDATA_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          CLSDATA_CREATE=0
        endif
      end function CLSDATA_CREATE

!---------------------------------------------------------
      subroutine CLSDATA_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if (CLSDATA_isValid(INST)) then
          call CLSDATA_DISPOSE_FIELDS(INST)
          DEallocate(ACLSDATA(inst)%X)
          nullify(ACLSDATA(inst)%X)
          CLSDATA_NC=CLSDATA_NC-1
          ACLSDATA(inst)%INST=0
        endif
      end subroutine CLSDATA_DISPOSE

!--------------------------------------------------------
      subroutine CLSDATA_DISPOSE_ALL
!--------------------------------------------------------
      integer :: i
      do i=1,CLSDATA_DIM
        call CLSDATA_DISPOSE(i)
      enddo
      CLSDATA_NC=0
      end subroutine CLSDATA_DISPOSE_ALL

!--------------------------------------------------------
      subroutine CLSDATA_ALLOC_FIELDS(OBJ,error)
! create new command instance from definition object TCLSDEF !
!--------------------------------------------------------
      type(TCLSDATA) :: OBJ
      logical,intent(out)           :: error
      integer :: i
        error=.true.
        OBJ%NARG=0
        do i=1,OBJ%CDEF%NF
          OBJ%FIELDS(i)%DEF=>OBJ%CDEF%FIELDS(i)%F
          if (ALLOC_FIELD(OBJ%FIELDS(i),OBJ%CDEF%FIELDS(i)%F%DIM)) then
            OBJ%NARG=OBJ%NARG+1
          else
            call MSG_ERROR('CLSDATA','cannot allocate field '//trim(OBJ%CDEF%ID),0,1)
            exit
          endif
        enddo
        error=(OBJ%CDEF%NF.ne.OBJ%NARG)
        if (error) call MSG_WARN('Some fields of '//trim(OBJ%CDEF%ID)//' where not allocated.',1)
      end subroutine CLSDATA_ALLOC_FIELDS

C---------------------------------------------------------
      SUBROUTINE CLSDATA_DISPOSE_FIELDS(INST)
C deallocate field buffers of the object
C---------------------------------------------------------
      integer,intent(in) :: INST
      integer :: i
      if (CLSDATA_isValid(INST)) then
        do i=1,ACLSDATA(INST)%X%NARG
          call DISPOSE_FIELD(ACLSDATA(INST)%X%FIELDS(i))
        enddo
      endif
      end SUBROUTINE CLSDATA_DISPOSE_FIELDS

C---------------------------------------------------------
      integer function getParIndex(CLSD,IDSTR)
! return command index for given ID string
! unique substring is handled as an exact match => abbreviations are allowed
C---------------------------------------------------------
      type(TCLSDATA) :: CLSD
      character(*) :: IDSTR
      integer i,IDCOMPARE
      logical :: found
      found=.false.
      i=0
      do WHILE((.not.found).and.(i.lt.CLSD%NARG))
        i=i+1
        found=(IDCOMPARE(trim(IDSTR),trim(CLSD%FIELDS(i)%ID)).eq.0)
      enddo
      if (found) then
        getParIndex=i
      else
        getParIndex=0
      endif
      end function getParIndex

!---------------------------------------------------------
      SUBROUTINE CLSDATA_INP(INST,PNAME,ARG,NARG)
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(TCLSDATA),pointer  :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: IP
      logical :: error
      NARG=0
  !    write(*,*) 'CLSDATA_INP  PNAME=',trim(PNAME)
      OBJ => ACLSDATA(INST)%X
      if (trim(PNAME).ne.trim(ARG%ID)) then
        write(*,*) 'CLSDATA_INP: unequal argument IDs: ['//trim(PNAME)//','//trim(ARG%ID)//']'
        return
      endif
      IP=getParIndex(OBJ,PNAME)
      if (IP.gt.0) then
        call FIELD_COPY(ARG,OBJ%FIELDS(IP),error)
        if (.not.error) NARG=OBJ%FIELDS(IP)%NP
      else
        write(*,*) 'CLSDATA_INP: undefined parameter : ['//trim(PNAME)//']'
      endif
      END SUBROUTINE CLSDATA_INP

!---------------------------------------------------------
      SUBROUTINE CLSDATA_OUT(INST,PNAME,ARG,NARG)
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(TCLSDATA),pointer  :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: IP
      logical :: error
      NARG=0
      if (trim(PNAME).ne.trim(ARG%ID)) then
        write(*,*) 'CLSDATA_OUT: unequal argument IDs: ['//trim(PNAME)//','//trim(ARG%ID)//']'
        return
      endif
      OBJ => ACLSDATA(INST)%X
      IP=getParIndex(OBJ,PNAME)
      if (IP.gt.0) then
        call FIELD_COPY(OBJ%FIELDS(IP),ARG,error)
        if (.not.error) NARG=ARG%NP
      else
        write(*,*) 'CLSDATA_OUT: undefined parameter : ['//trim(PNAME)//']'
      endif
      END SUBROUTINE CLSDATA_OUT

!--------------------------------------------------------
      subroutine send_message(what,ID,TID)
!--------------------------------------------------------
      character(*),intent(in) :: what,ID,TID
      select case(what)
      case('NR')
        call MSG_WARN('cannot read numerical value from: '//trim(ID)//':'//trim(tid),1)
      case('NA')
        call MSG_WARN('cannot assign numerical value to: '//trim(ID)//':'//trim(tid),1)
      case('SR')
        call MSG_WARN('cannot read string value from: '//trim(ID)//':'//trim(tid),1)
      case('SA')
        call MSG_WARN('cannot assign string value to: '//trim(ID)//':'//trim(tid),1)
      case('P')
        call MSG_WARN('unknown command parameter: '//trim(ID),1)
      case('C')
        call MSG_WARN('unknown command ID: '//trim(ID),1)
      case('ASSIGN_ICLS')
        call MSG_WARN('CLSDATA_ASSIGN: wrong ICLS for '//trim(ID)//':'//trim(tid),1)
      case('ASSIGN_ICMD')
        call MSG_WARN('CLSDATA_ASSIGN: wrong ICMD for '//trim(ID)//':'//trim(tid),1)
      case('ADD_FIELD')
        call MSG_WARN('CLSDATA_ASSIGN: cannot add field '//trim(ID)//':'//trim(tid),1)
      case default
        call MSG_WARN('undefined error: '//trim(what),1)
      end select

      !READ(*,*)
      end subroutine send_message

      end module CLASSDATA
