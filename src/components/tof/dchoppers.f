!//////////////////////////////////////////////////////////////////////
!////  $Id $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2016/03/02 11:30:36 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: DCHOPPER
!////
!////////////////////////////////////////////////////////////////////////
      MODULE DCHOPPERS
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      implicit none
      private
      save


      integer,parameter :: MAXWINS=64
      TYPE  DCHOPPER; SEQUENCE
        TYPE(TFRAME) :: FRAME
        real(kind(1.D0)) :: FRQ          ! frequency [MHz]
        real(kind(1.D0)) :: PHI          ! phase [2PI]
        real(kind(1.D0)) :: RAD          ! radius [mm]
        real(kind(1.D0)) :: WIN          ! disc window size (sector angle) [2PI]
        INTEGER :: NW      ! number of windows
        INTEGER :: PROF    ! window profile (0) V-shape, (1) U-shape
        INTEGER :: ORI     ! orientation (0) vertical (window above the axis) (1) horizontal
        LOGICAL :: AUTOADJ ! automatically adjust phase for given wavelength
        LOGICAL :: isT0    ! is T0 chopper
        LOGICAL :: lockT0    ! is T0 chopper
        LOGICAL :: OVERLAP ! allow frame overlaps
        character(1024) ::  FRAMES ! multi-frame data,  "phase1|width1|phase2|width2|..."
      ! calculated fields
        REAL (kind(1.D0)) :: T0 ! nominal TOF to the chopper
		LOGICAL :: ACTIVE  ! Ignore the chopper if not
      ! data for multi-frame chopper
        REAL (kind(1.D0)) :: PHASES(MAXWINS)
        REAL (kind(1.D0)) :: WIDTHS(MAXWINS) ! widths as fractions of WIN
        real(kind(1.D0)) ::  TANWS(MAXWINS) ! tan(WIN/2)
      END TYPE DCHOPPER

      integer, parameter :: DCHOPPERS_DIM=32
! pointer type for DCHOPPER
      type PDCHOPPER
        integer :: IDX  ! instance index
        TYPE(DCHOPPER),pointer :: X
      end type PDCHOPPER
! instances of DCHOPPER. ADCHOPPERS(0) is always unallocated
      integer :: DCHOPPERS_NC
      type(PDCHOPPER) :: ADCHOPPERS(1:DCHOPPERS_DIM)

      public DCHOPPER,PDCHOPPER,ADCHOPPERS,DCHOPPER_PREPARE
      public DCHOPPER_CREATE,DCHOPPER_DISPOSE,DCHOPPER_DISPOSE_ALL,AddDCHOPPER
      public DCHOPPER_isValid,DCHOPPER_DEFAULT,DCHOPPER_GET,DCHOPPER_GET_FRAME
      public DCHOPPER_INP,DCHOPPER_OUT,DCHOPPER_READ_FRAMES

      contains

!-------------------------------------------------------------------------
! Creator for DCHOPPER, return instance number
! Memory is allocated on the first free element of ADCHOPPERS array
!-------------------------------------------------------------------------
      integer function DCHOPPER_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.DCHOPPERS_DIM).and.(ADCHOPPERS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (ADCHOPPERS(i)%IDX.le.0) then
          allocate(ADCHOPPERS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            DCHOPPERS_NC=DCHOPPERS_NC+1
            call DCHOPPER_DEFAULT(i)
            ADCHOPPERS(i)%IDX=i
            ADCHOPPERS(i)%X%FRAME%ID=trim(ID)
            ADCHOPPERS(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          DCHOPPER_CREATE=i
        else
          call INT2STR(DCHOPPERS_DIM,S)
          call MSG_ERROR('DCHOPPER_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          DCHOPPER_CREATE=0
        endif
      end function DCHOPPER_CREATE

!---------------------------------------------------------
      subroutine DCHOPPER_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.DCHOPPERS_DIM)) then
          if (associated(ADCHOPPERS(inst)%X)) then
            DEallocate(ADCHOPPERS(inst)%X)
            nullify(ADCHOPPERS(inst)%X)
            DCHOPPERS_NC=DCHOPPERS_NC-1
          endif
          ADCHOPPERS(inst)%IDX=0
        endif
      end subroutine DCHOPPER_DISPOSE

!---------------------------------------------------------
      subroutine DCHOPPER_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,DCHOPPERS_DIM
          call DCHOPPER_DISPOSE(i)
        enddo
      end subroutine DCHOPPER_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddDCHOPPER(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (DCHOPPER_isValid(INST)) then
        i=DCHOPPER_CREATE(ADCHOPPERS(INST)%X%FRAME%ID,ADCHOPPERS(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          ADCHOPPERS(i)%X=ADCHOPPERS(INST)%X
          ADCHOPPERS(i)%IDX=i
        !  write(*,*) 'AddDCHOPPER ',trim(ADCHOPPERS(i)%X%FRAME%ID),i
        endif
      endif
      AddDCHOPPER=i
      end function AddDCHOPPER

!------------------------------------
      SUBROUTINE DCHOPPER_PREPARE(INST,IERR)
! prepare dependent fields after input
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(DCHOPPER),POINTER :: OBJ
        IERR=0
        if (.not.DCHOPPER_isValid(INST)) return
        OBJ => ADCHOPPERS(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE DCHOPPER_PREPARE

!-------------------------------------------------------------
      logical function DCHOPPER_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      DCHOPPER_isValid= ((INST.gt.0).and.(INST.le.DCHOPPERS_DIM).and.associated(ADCHOPPERS(INST)%X))
      end function DCHOPPER_isValid

!------------------------------------
      SUBROUTINE DCHOPPER_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
        if (.not.DCHOPPER_isValid(INST)) return
        OBJ => ADCHOPPERS(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%SIZE(3)=1.D0
        OBJ%FRAME%CLASS=CCLS_DCHOPPER
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%FRQ=1.D-4 ! = 100 Hz
        OBJ%PHI=0.D0
        OBJ%RAD=150.D0
        OBJ%WIN=0.1 ! = 10% of the period
        OBJ%NW=1
        OBJ%PROF=0
        OBJ%ORI=0
        OBJ%AUTOADJ=.true.
        OBJ%isT0=.false.
        OBJ%lockT0=.false.
        OBJ%OVERLAP=.false.
        OBJ%T0=0.D0
		OBJ%ACTIVE=.true.
      END SUBROUTINE DCHOPPER_DEFAULT

!-------------------------------------------------------------
      SUBROUTINE DCHOPPER_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PDCHOPPER),intent(out) :: OBJ
        if (DCHOPPER_isValid(INST)) then
          OBJ%IDX=ADCHOPPERS(INST)%IDX
          OBJ%X => ADCHOPPERS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE DCHOPPER_GET

!-------------------------------------------------------------
      SUBROUTINE DCHOPPER_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (DCHOPPER_isValid(INST)) then
          OBJ%IDX=ADCHOPPERS(INST)%IDX
          OBJ%X => ADCHOPPERS(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE DCHOPPER_GET_FRAME

!---------------------------------------------------------
      SUBROUTINE DCHOPPER_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      if (.not.DCHOPPER_isValid(INST)) return
      OBJ => ADCHOPPERS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call DCHOPPER_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case(ID_FIELD_STR)
          call FIELD2STR(ARG,0,SARG,LR)
          call DCHOPPER_INP_S(OBJ,trim(ARG%ID),SARG,LR)
        case default
          write(*,*) 'DCHOPPER_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE DCHOPPER_INP


!---------------------------------------------------------
      SUBROUTINE DCHOPPER_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DCHOPPER),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,NA
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      if (.not.DCHOPPER_isValid(INST)) return
      OBJ => ADCHOPPERS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call DCHOPPER_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case(ID_FIELD_STR)
          call DCHOPPER_OUT_S(OBJ,trim(ARG%ID),SARG,NA)
          call STR2FIELD(ARG,0,SARG,LR)
        case default
          write(*,*) 'DCHOPPER_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE DCHOPPER_OUT

!---------------------------------------------------------
      SUBROUTINE DCHOPPER_READ_FRAMES(OBJ)
! read multi-frame data from the OBJ%FRAMES string
! syntax is "phase1|width1|phase2|width2|..."
!---------------------------------------------------------
      TYPE(DCHOPPER),intent(out) :: OBJ
      real :: AUX(MAXWINS*2+1)
      integer :: N,N1,i
      character*32 :: CNUM1,CNUM2
      !write(*,*) 'DCHOPPER_RF [',trim(OBJ%FRAMES),']'
      call STR2ARRAY4(OBJ%FRAMES,':',AUX,MAXWINS*2+1,N)
1     format(a,6(G12.5))
      N=MIN(MAXWINS*2,N)
      !write(*,1) 'DCHOPPER_RF ' ,N
      N1=N/2
      if (N>1) then
        DO i=0,N1-1
          OBJ%PHASES(i+1)=AUX(2*i+1)
          OBJ%WIDTHS(i+1)=AUX(2*i+2)
          !write(*,1) 'DCHOPPER_RF ' ,i,OBJ%PHASES(i+1)
        enddo
        OBJ%NW=N1
    ! if no valid description found, use equidistant frames
      else
        OBJ%FRAMES=''
        DO i=0,OBJ%NW-1
          OBJ%PHASES(i+1)=i*1.D0/OBJ%NW
          OBJ%WIDTHS(i+1)=1.D0
        enddo
      endif
      ! rewrite FRAMES string to the canonical form
      OBJ%FRAMES=''
      DO i=1,OBJ%NW
        call FLOAT2STR(OBJ%PHASES(i),CNUM1)
        call FLOAT2STR(OBJ%WIDTHS(i),CNUM2)
        call APPENDPAR(OBJ%FRAMES,':',trim(CNUM1))
        if (i==OBJ%NW) then
          call APPENDPAR(OBJ%FRAMES,' ',trim(CNUM2))
        else
          call APPENDPAR(OBJ%FRAMES,':',trim(CNUM2))
        endif
      enddo
      end SUBROUTINE DCHOPPER_READ_FRAMES

!---------------------------------------------------------
      SUBROUTINE DCHOPPER_INP_S(OBJ,PNAME,SARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DCHOPPER),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(in) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('TIMING')
          OBJ%FRAMES=trim(SARG)
      END SELECT
      NARG=LR
      END SUBROUTINE DCHOPPER_INP_S


!---------------------------------------------------------
      SUBROUTINE DCHOPPER_OUT_S(OBJ,PNAME,SARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DCHOPPER),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*), intent(out) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=0
      SELECT CASE (trim(PNAME))
        CASE('TIMING')
          SARG=' '
          if (len_trim(OBJ%FRAMES).gt.0) SARG=trim(OBJ%FRAMES)
          LR=1
      END SELECT
      NARG=LR
      END SUBROUTINE DCHOPPER_OUT_S


!---------------------------------------------------------
      SUBROUTINE DCHOPPER_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DCHOPPER),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
      case('FRQ');   OBJ%FRQ=ARG(1)*1.D-6
      case('PHASE'); OBJ%PHI=ARG(1)
      case('RAD');   OBJ%RAD=ARG(1)
      case('WIN');   OBJ%WIN=ARG(1)
      case('NWIN');  OBJ%NW=MIN(MAXWINS,MAX(NINT(ARG(1)),1))
      case('PROF');  OBJ%PROF=MIN(NINT(abs(ARG(1))),1)
      case('ORI');   OBJ%ORI=MIN(NINT(abs(ARG(1))),1)
      case('ADJ');   OBJ%AUTOADJ=(NINT(ARG(1)).eq.1)
      case('T0');    OBJ%isT0=(NINT(ARG(1)).eq.1)
      case('LOCKT0');    OBJ%lockT0=(NINT(ARG(1)).eq.1)
      case('OVERLAP'); OBJ%OVERLAP=(NINT(ARG(1)).eq.1)
! try FRAME parameters
      CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE DCHOPPER_INP_R

!---------------------------------------------------------
      SUBROUTINE DCHOPPER_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DCHOPPER),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
      case('FRQ');  ARG(1)=OBJ%FRQ*1.D6
      case('PHASE'); ARG(1)=OBJ%PHI
      case('RAD');   ARG(1)=OBJ%RAD
      case('WIN');   ARG(1)=OBJ%WIN
      case('NWIN');  ARG(1)=OBJ%NW
      case('PROF');  ARG(1)=OBJ%PROF
      case('ORI');   ARG(1)=OBJ%ORI
      case('ADJ');   ARG(1)=LOG2INT(OBJ%AUTOADJ)
      case('T0');    ARG(1)=LOG2INT(OBJ%isT0)
      case('LOCKT0');    ARG(1)=LOG2INT(OBJ%lockT0)
      case('OVERLAP');   ARG(1)=LOG2INT(OBJ%OVERLAP)
! try FRAME parameters
          CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE DCHOPPER_OUT_R

      end MODULE DCHOPPERS
