!//////////////////////////////////////////////////////////////////////
!////  $Id: igroups.f,v 1.6 2012/02/01 13:50:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.6 $
!////     $Date: 2012/02/01 13:50:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes IGROUP, a group of interacting components
!////
!////////////////////////////////////////////////////////////////////////
      MODULE IGROUPS
      use CLASSES
      use FIELDDEF
      use FIELDDATA
      use VCTABLE_COMPONENTS
      implicit none


! static storage of IGROUP instances
      integer, parameter :: IGROUPS_DIM=5
      integer,parameter :: MAX_IGROUP_ITEMS=128
      TYPE TIGROUP
        TYPE(TFRAME) :: FRAME
        integer :: NCOM  ! nuber of member components
        TYPE(TCOMOBJ) :: ITEMS(MAX_IGROUP_ITEMS) ! instance indices of member components (see VCTABLE_COMPONENTS)
      END TYPE TIGROUP

      type PIGROUP
        integer :: IDX  ! instance index
        TYPE(TIGROUP),pointer :: X
      end type PIGROUP

      integer :: IGROUPS_NC
      type(PIGROUP) :: AIGROUPS(IGROUPS_DIM)

      contains

!-------------------------------------------------------------------------
! Creator for IGROUP, return instance number
! Memory is allocated on the first free element of AIGROUPS array
!-------------------------------------------------------------------------
      integer function IGROUP_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.IGROUPS_DIM).and.(AIGROUPS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (AIGROUPS(i)%IDX.le.0) then
          allocate(AIGROUPS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            IGROUPS_NC=IGROUPS_NC+1
            call IGROUP_CLEAR(i)
            AIGROUPS(i)%IDX=i
            AIGROUPS(i)%X%FRAME%ID=trim(ID)
            AIGROUPS(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          IGROUP_CREATE=i
        else
          call INT2STR(IGROUPS_DIM,S)
          call MSG_ERROR('IGROUP_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          IGROUP_CREATE=0
        endif
      end function IGROUP_CREATE

!---------------------------------------------------------
      subroutine IGROUP_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.IGROUPS_DIM)) then
          if (associated(AIGROUPS(inst)%X)) then
            DEallocate(AIGROUPS(inst)%X)
            nullify(AIGROUPS(inst)%X)
            IGROUPS_NC=IGROUPS_NC-1
          endif
          AIGROUPS(inst)%IDX=0
        endif
      end subroutine IGROUP_DISPOSE

!---------------------------------------------------------
      subroutine IGROUP_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,IGROUPS_DIM
          call IGROUP_DISPOSE(i)
        enddo
      end subroutine IGROUP_DISPOSE_ALL

!-------------------------------------------------------------
      SUBROUTINE IGROUP_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PIGROUP),intent(out) :: OBJ
        if (IGROUP_isValid(INST)) then
          OBJ%IDX=AIGROUPS(INST)%IDX
          OBJ%X => AIGROUPS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE IGROUP_GET

!-------------------------------------------------------------
      logical function IGROUP_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      IGROUP_isValid= ((INST.gt.0).and.(INST.le.IGROUPS_DIM).and.associated(AIGROUPS(INST)%X))
      end function IGROUP_isValid

!-------------------------------------------
      SUBROUTINE IGROUP_ADD_MEMBER(INST,OBJ)
! set default parameters
!-------------------------------------------
      INTEGER,intent(in) :: INST
      TYPE(TCOMOBJ),intent(in) :: OBJ
      integer :: N
      if (IGROUP_isValid(INST)) then
        N=AIGROUPS(INST)%X%NCOM
        if (N.lt.MAX_IGROUP_ITEMS) then
          AIGROUPS(INST)%X%ITEMS(N+1)=OBJ
          AIGROUPS(INST)%X%NCOM=N+1
        endif
      endif
      END SUBROUTINE IGROUP_ADD_MEMBER

!------------------------------------
      SUBROUTINE IGROUP_CLEAR(INST)
! clear contents of the group
!------------------------------------
      integer,intent(in) :: INST
      if (IGROUP_isValid(INST)) then
        call FRAME_CLEAR(AIGROUPS(INST)%X%FRAME)
        AIGROUPS(INST)%X%NCOM=0
      endif
      END SUBROUTINE IGROUP_CLEAR

!---------------------------------------------------------
      SUBROUTINE IGROUP_RESET(INST)
! clear event counter of a component
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      integer :: i
      if (IGROUP_isValid(INST)) then
        AIGROUPS(INST)%X%FRAME%COUNT=0
        do i=1,AIGROUPS(INST)%X%NCOM
          call PCOMOBJ_RESET(AIGROUPS(INST)%X%ITEMS(i))
        enddo
      endif
      end SUBROUTINE IGROUP_RESET


!----------------------------------------------------
      SUBROUTINE IGROUP_INIT(INST)
! call INIT procedure on all objects in the group
!----------------------------------------------------
      INTEGER,intent(in) :: INST
      integer :: i
      if (IGROUP_isValid(INST)) then
        call FRAME_INIT(AIGROUPS(INST)%X%FRAME)
        do i=1,AIGROUPS(INST)%X%NCOM
          call PCOMOBJ_INIT(AIGROUPS(INST)%X%ITEMS(i))
        enddo
      endif
      END SUBROUTINE IGROUP_INIT

!---------------------------------------------------------
       SUBROUTINE IGROUP_FINALIZE(INST)
! dispatch calls to *_FINALIZE procedures to all group members
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      integer :: i
      if (IGROUP_isValid(INST)) then
        do i=1,AIGROUPS(INST)%X%NCOM
          call PCOMOBJ_FINALIZE(AIGROUPS(INST)%X%ITEMS(i))
        enddo
      endif
      end SUBROUTINE IGROUP_FINALIZE

!---------------------------------------------------------
       SUBROUTINE IGROUP_ENDTRACE(INST,NEU)
! dispatch calls to *_FINALIZE procedures to all group members
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      TYPE(NEUTRON) :: NEU
      integer :: i
      if (IGROUP_isValid(INST)) then
        do i=1,AIGROUPS(INST)%X%NCOM
          call PCOMOBJ_ENDTRACE(AIGROUPS(INST)%X%ITEMS(i),NEU)
        enddo
      endif
      end SUBROUTINE IGROUP_ENDTRACE

!---------------------------------------------------------
       SUBROUTINE IGROUP_ADJUST(INST,IERR)
! dispatch calls to *_ADJUST procedures for each object
! IERR<>0 if there is a problem with positioning
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      integer,intent(out) :: IERR
      integer :: i,ie
      ierr=0
      if (IGROUP_isValid(INST)) then
        do i=1,AIGROUPS(INST)%X%NCOM
          call PCOMOBJ_ADJUST(AIGROUPS(INST)%X%ITEMS(i),IE)
          if (ie.ne.0) ierr=1
        enddo
      endif
      end SUBROUTINE IGROUP_ADJUST


!---------------------------------------------------------
      logical function  IGROUP_GO(INST)
! dispatch calls to *_FINALIZE procedures to all group members
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      integer :: i
      logical :: RES
      RES=.true.
      if (IGROUP_isValid(INST)) then
        do i=1,AIGROUPS(INST)%X%NCOM
          RES=PCOM_GO(AIGROUPS(INST)%X%ITEMS(i))
          if (.not.RES) exit
        enddo
      endif
      IGROUP_GO=RES
      end function IGROUP_GO

!---------------------------------------------------------
      SUBROUTINE IGROUP_INP(INST,PNAME,ARG,NARG)
! call INP on group FRAME object
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      if (IGROUP_isValid(INST)) then
        call FRAME_INP(AIGROUPS(INST)%X%FRAME,PNAME,ARG,NARG)
      endif
      END SUBROUTINE IGROUP_INP

!---------------------------------------------------------
      SUBROUTINE IGROUP_OUT(INST,PNAME,ARG,NARG)
! call OUT on group FRAME object
!---------------------------------------------------------
      INTEGER,intent(in) :: INST
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      if (IGROUP_isValid(INST)) then
        call FRAME_OUT(AIGROUPS(INST)%X%FRAME,PNAME,ARG,NARG)
      endif
      END SUBROUTINE IGROUP_OUT

      end module IGROUPS