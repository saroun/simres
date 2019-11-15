!//////////////////////////////////////////////////////////////////////
!////  $Id: monitors.f,v 1.10 2015/04/13 18:38:51 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights AREAerved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.10 $
!////     $Dat$
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: MONITOR
!////
!////////////////////////////////////////////////////////////////////////
      MODULE MONITORS
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      implicit none
      save
      private

      integer,parameter,public :: monit_mode_add=0
      integer,parameter,public :: monit_mode_inner=1
      integer,parameter,public :: monit_mode_outer=2
      integer,parameter,public :: monit_mode_sub=3 ! subtract inner-outer

      integer :: LREC=4 ! record length for events
      integer :: MREC=10000 ! segment size in records for memory allocation

      TYPE  MONITOR; SEQUENCE
        TYPE(TFRAME) :: FRAME
        INTEGER :: CNTMODE ! counting mode
        INTEGER :: IX,IY ! variable type (phase space variable as in BEAM1D command)
        REAL(kind(1.D0)) :: X0,DX   ! horizontal clip range
        REAL(kind(1.D0)) :: Y0,DY   ! vertical clip range
        ! calculated fields
        REAL(kind(1.D0)) :: IIN,IOUT   ! sums outside and inside clip area
        REAL(kind(1.D0)) :: K0   ! nominal |K|
        logical :: BLOCKAREA ! if true, block events passing outside the clip area
      END TYPE MONITOR

      integer, parameter :: MONITORS_DIM=6
! pointer type for MONITOR
      type PMONITOR
        integer :: IDX  ! instance index
        TYPE(MONITOR),pointer :: X
      end type PMONITOR
! instances of MONITOR. AMONITORS(0) is always unallocated
      integer :: MONITORS_NC
      type(PMONITOR) :: AMONITORS(1:MONITORS_DIM)

      public MONITOR_CREATE,MONITOR_DISPOSE,MONITOR_DISPOSE_ALL,AddMONITOR
      public MONITOR,PMONITOR,AMONITORS,MONITOR_GET,MONITOR_PREPARE
      public MONITOR_DEFAULT,MONITOR_INP,MONITOR_OUT,MONITOR_GET_FRAME,MONITOR_isValid

      public MONITORS_DIM,MONITORS_NC

      contains


!-------------------------------------------------------------------------
! Creator for MONITOR, return instance number
! Memory is allocated on the first free element of AMONITORS array
!-------------------------------------------------------------------------
      integer function MONITOR_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.MONITORS_DIM).and.(AMONITORS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (AMONITORS(i)%IDX.le.0) then
          allocate(AMONITORS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            MONITORS_NC=MONITORS_NC+1
            call MONITOR_DEFAULT(i)
            AMONITORS(i)%IDX=i
            AMONITORS(i)%X%FRAME%ID=trim(ID)
            AMONITORS(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          MONITOR_CREATE=i
        else
          call INT2STR(MONITORS_DIM,S)
          call MSG_ERROR('MONITOR_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          MONITOR_CREATE=0
        endif
      end function MONITOR_CREATE


!---------------------------------------------------------
      subroutine MONITOR_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.MONITORS_DIM)) then
          if (associated(AMONITORS(inst)%X)) then
            DEallocate(AMONITORS(inst)%X)
            nullify(AMONITORS(inst)%X)
            MONITORS_NC=MONITORS_NC-1
          endif
          AMONITORS(inst)%IDX=0
        endif
      end subroutine MONITOR_DISPOSE

!---------------------------------------------------------
      subroutine MONITOR_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,MONITORS_DIM
          call MONITOR_DISPOSE(i)
        enddo
      end subroutine MONITOR_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddMONITOR(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (MONITOR_isValid(INST)) then
        i=MONITOR_CREATE(AMONITORS(INST)%X%FRAME%ID,AMONITORS(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          AMONITORS(i)%X=AMONITORS(INST)%X
          AMONITORS(i)%IDX=i
        !  write(*,*) 'AddMONITOR ',trim(AMONITORS(i)%X%FRAME%ID),i
        endif
      endif
      AddMONITOR=i
      end function AddMONITOR


!-------------------------------------------------------------
      logical function MONITOR_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      MONITOR_isValid= ((INST.gt.0).and.(INST.le.MONITORS_DIM).and.associated(AMONITORS(INST)%X))
      end function MONITOR_isValid


!-------------------------------------------------------------
      SUBROUTINE MONITOR_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PMONITOR) :: OBJ
        if (MONITOR_isValid(INST)) then
          OBJ%IDX=AMONITORS(INST)%IDX
          OBJ%IDX=2
          OBJ%X => AMONITORS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE MONITOR_GET

!------------------------------------
      SUBROUTINE MONITOR_PREPARE(INST,IERR)
! prepare dependent fields after input
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(MONITOR),POINTER :: OBJ
        IERR=0
        if (.not.MONITOR_isValid(INST)) return
        OBJ => AMONITORS(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE MONITOR_PREPARE

!-------------------------------------------------------------
      SUBROUTINE MONITOR_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (MONITOR_isValid(INST)) then
          OBJ%IDX=AMONITORS(INST)%IDX
          OBJ%X => AMONITORS(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE MONITOR_GET_FRAME

!------------------------------------
      SUBROUTINE MONITOR_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(MONITOR),POINTER :: OBJ
        if (.not.MONITOR_isValid(INST)) return
        OBJ => AMONITORS(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%CLASS=CCLS_MONITOR
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%FRAME%SIZE=(/64.D0,64.D0,1.D0/)
        OBJ%CNTMODE=0
        OBJ%IIN=0.D0
        OBJ%IOUT=0.D0
        OBJ%BLOCKAREA=.false.
      END SUBROUTINE MONITOR_DEFAULT

!---------------------------------------------------------
      SUBROUTINE MONITOR_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(MONITOR),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.MONITOR_isValid(INST)) return
      OBJ => AMONITORS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call MONITOR_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'MONITOR_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE MONITOR_INP

!---------------------------------------------------------
      SUBROUTINE MONITOR_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(MONITOR),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.MONITOR_isValid(INST)) return
      OBJ => AMONITORS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call MONITOR_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'MONITOR_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE MONITOR_OUT

!---------------------------------------------------------
      SUBROUTINE MONITOR_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(MONITOR),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! counting mode
        CASE('CMODE')
          OBJ%CNTMODE=MAX(0,MIN(NINT(ARG(1)),3))
! variable type
          CASE('IX'); OBJ%IX=NINT(ARG(1))+1
          CASE('IY'); OBJ%IY=NINT(ARG(1))+1
! inner area [arb. units]
          CASE('X0');  OBJ%X0=ARG(1)
          CASE('DX');  OBJ%DX=ARG(1)
          CASE('Y0');  OBJ%Y0=ARG(1)
          CASE('DY');  OBJ%DY=ARG(1)
          CASE('BLOCK');  OBJ%BLOCKAREA=(NINT(ARG(1)).eq.1)
! try FRAME parameters
        CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE MONITOR_INP_R

!---------------------------------------------------------
      SUBROUTINE MONITOR_OUT_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(MONITOR),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! counting mode
          CASE('CMODE'); ARG(1)=OBJ%CNTMODE
! variable type
          CASE('IX'); ARG(1)=OBJ%IX-1
          CASE('IY'); ARG(1)=OBJ%IY-1
! inner area [arb. units]
          CASE('X0');   ARG(1)=OBJ%X0
          CASE('DX');   ARG(1)=OBJ%DX
          CASE('Y0');   ARG(1)=OBJ%Y0
          CASE('DY');   ARG(1)=OBJ%DY
          CASE('BLOCK');   ARG(1)=LOG2INT(OBJ%BLOCKAREA)
! try FRAME parameters
          CASE DEFAULT
            CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE MONITOR_OUT_R

      end MODULE MONITORS