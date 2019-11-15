!//////////////////////////////////////////////////////////////////////
!////  $Id: sources.f,v 1.29 2016/11/13 15:20:31 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.29 $
!////     $Date: 2016/11/13 15:20:31 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: SOURCE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SOURCES
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      use SOURCES_TABLE
      use SOURCES_PULSE
      implicit none
      save
      private

! constants
      integer,parameter,PUBLIC :: src_steady=0
      integer,parameter,PUBLIC :: src_pulsed=1

! SOURCE
      TYPE SOURCE; SEQUENCE
        TYPE(TFRAME) :: FRAME
        real(kind(1.D0)) :: FLUX ! integrated isotropic flux in 10^14/cm2/s
        real(kind(1.D0))  :: TEMP
        real(kind(1.D0)) :: FRQ    ! frequency [MHz]
        real(kind(1.D0)) :: DELAY   ! time delay in [ms]
        real(kind(1.D0)) :: PULSW ! pulse width [ms]
        real(kind(1.D0)) :: LAMW ! wavelength band (d_lambda/lambda)
        logical :: OVERLAP ! allow frame overlap
        logical :: is_pulsed ! set to true by SOURCE_ADJUST if the source is time limited
        integer :: TYP    ! steady-state (0) or pulsed (1)
        CHARACTER(256) :: FLUXTAB  ! filename with flux distribution table
      END TYPE SOURCE
      integer, parameter :: SOURCES_DIM=2
! pointer type for SOURCE
      type PSOURCE
        integer :: IDX  ! instance index
        TYPE(SOURCE),pointer :: X
      end type PSOURCE
! instances of SOURCE. ASOURCES(0) is always unallocated
      integer :: SOURCES_NC
      type(PSOURCE) :: ASOURCES(1:SOURCES_DIM)

      public SOURCE,PSOURCE,ASOURCES,SOURCE_PREPARE
      public SOURCES_NC
      public SOURCE_DEFAULT,SOURCE_INP,SOURCE_OUT,SOURCE_GET_FRAME,SOURCE_isValid
      public SOURCE_CREATE,SOURCE_DISPOSE,SOURCE_DISPOSE_ALL,AddSOURCE,SOURCE_GET

      contains

!-------------------------------------------------------------------------
! Creator for SOURCE, return instance number
! Memory is allocated on the first free element of ASOURCES array
!-------------------------------------------------------------------------
      integer function SOURCE_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.SOURCES_DIM).and.(ASOURCES(i)%IDX.gt.0))
          i=i+1
        enddo
        if (ASOURCES(i)%IDX.le.0) then
          allocate(ASOURCES(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            SOURCES_NC=SOURCES_NC+1
            call SOURCE_DEFAULT(i)
            ASOURCES(i)%IDX=i
            ASOURCES(i)%X%FRAME%ID=trim(ID)
            ASOURCES(i)%X%FRAME%NAME=trim(NAMESTR)
          !  write(*,*) 'SOURCE_CREATE ',i,' NC=',SOURCES_NC
          endif
        endif
        if (ierr.eq.0) then
          SOURCE_CREATE=i
        else
          call INT2STR(SOURCES_DIM,S)
          call MSG_ERROR('SOURCE_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          SOURCE_CREATE=0
        endif
      end function SOURCE_CREATE

!---------------------------------------------------------
      subroutine SOURCE_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.SOURCES_DIM)) then
          if (associated(ASOURCES(inst)%X)) then
            DEallocate(ASOURCES(inst)%X)
            nullify(ASOURCES(inst)%X)
            SOURCES_NC=SOURCES_NC-1
          !  write(*,*) 'SOURCE_DISPOSE ',INST,' NC=',SOURCES_NC
          endif
          ASOURCES(inst)%IDX=0
        endif
      end subroutine SOURCE_DISPOSE

!---------------------------------------------------------
      subroutine SOURCE_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,SOURCES_DIM
          call SOURCE_DISPOSE(i)
        enddo
        SOURCES_NC=0
      end subroutine SOURCE_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddSOURCE(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
    !  write(*,*) 'AddSOURCE ',INST
      if (SOURCE_isValid(INST)) then
        i=SOURCE_CREATE(ASOURCES(INST)%X%FRAME%ID,ASOURCES(INST)%X%FRAME%NAME)
      !  write(*,*) 'AddSOURCE created',i
        if (i.gt.0) then
          ASOURCES(i)%X=ASOURCES(INST)%X
          ASOURCES(i)%IDX=i
        !  write(*,*) 'AddSOURCE ',trim(ASOURCES(i)%X%FRAME%ID),i
        endif
      endif
      AddSOURCE=i
      end function AddSOURCE

!------------------------------------
      SUBROUTINE SOURCE_PREPARE(INST,IERR)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(SOURCE),POINTER :: OBJ
        IERR=0
        if (.not.SOURCE_isValid(INST)) return
        OBJ => ASOURCES(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE SOURCE_PREPARE

!-------------------------------------------------------------
      logical function SOURCE_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      SOURCE_isValid= ((INST.gt.0).and.(INST.le.SOURCES_DIM).and.associated(ASOURCES(INST)%X))
      end function SOURCE_isValid

!-------------------------------------------------------------
      SUBROUTINE SOURCE_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PSOURCE),intent(out) :: OBJ
        if (SOURCE_isValid(INST)) then
          OBJ%IDX=ASOURCES(INST)%IDX
          OBJ%X => ASOURCES(INST)%X
        else
          write(*,*) 'SOURCE_GET not valid: ',INST
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE SOURCE_GET

!-------------------------------------------------------------
      SUBROUTINE SOURCE_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (SOURCE_isValid(INST)) then
          OBJ%IDX=ASOURCES(INST)%IDX
          OBJ%X => ASOURCES(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE SOURCE_GET_FRAME

!------------------------------------
      SUBROUTINE SOURCE_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(SOURCE),POINTER :: OBJ
        if (.not.SOURCE_isValid(INST)) return
        OBJ => ASOURCES(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%CLASS=CCLS_SOURCE
        OBJ%FRAME%SHAPE=FRAME_SHAPE_DISC
        OBJ%FRAME%SIZE=(/100.D0,100.D0,1.D-1/)
        OBJ%FLUX=1.D0 ! x 10^14 n/s/cm^2
        OBJ%TEMP=310 ! temperature in [K]
        OBJ%FLUXTAB='none' ! lookup table name, no by default
        OBJ%PULSW=0.D0
        OBJ%DELAY=0.D0
        OBJ%LAMW=0.01
        OBJ%TYP=0
        OBJ%OVERLAP=.false.
        OBJ%FRQ=1.D-6
        OBJ%is_pulsed=.false.
      END SUBROUTINE SOURCE_DEFAULT

!---------------------------------------------------------
      SUBROUTINE SOURCE_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SOURCE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      if (.not.SOURCE_isValid(INST)) return
      OBJ => ASOURCES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call SOURCE_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case(ID_FIELD_STR)
          call FIELD2STR(ARG,0,SARG,LR)
          call SOURCE_INP_S(OBJ,trim(ARG%ID),SARG,LR)
        case default
          write(*,*) 'SOURCE_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SOURCE_INP


!---------------------------------------------------------
      SUBROUTINE SOURCE_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(SOURCE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,NA
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      if (.not.SOURCE_isValid(INST)) return
      OBJ => ASOURCES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call SOURCE_OUT_R(OBJ,trim(ARG%ID),NUM,NA)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case(ID_FIELD_STR)
          call SOURCE_OUT_S(OBJ,trim(ARG%ID),SARG,NA)
          call STR2FIELD(ARG,0,SARG,LR)
        case default
          write(*,*) 'SOURCE_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
  !    write(*,*) 'SOURCE_OUT ',trim(PNAME),' LR= ',LR,' NA=',NA
      END SUBROUTINE SOURCE_OUT


!---------------------------------------------------------
      SUBROUTINE SOURCE_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SOURCE),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
    !  if (INDEX(PNAME,'FLUX').gt.0) write(*,*) '  SOURCE_INP [',trim(PNAME),']'
      SELECT CASE (trim(PNAME))
        CASE('FLUX');   OBJ%FLUX=ARG(1)
        CASE('TEMP');   OBJ%TEMP=ARG(1)
        case('FRQ');   OBJ%FRQ=ARG(1)*1.D-6
        CASE('DELAY');  OBJ%DELAY=ARG(1)*1.D3
        CASE('LAMW');   OBJ%LAMW=ARG(1)
        CASE('PULSW');  OBJ%PULSW=ARG(1)*1.D3
        case('TYPE');    OBJ%TYP=MIN(NINT(abs(ARG(1))),1)
        case('OVERLAP'); OBJ%OVERLAP=(NINT(ARG(1)).eq.1)
        CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SOURCE_INP_R

!---------------------------------------------------------
      SUBROUTINE SOURCE_INP_S(OBJ,PNAME,SARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SOURCE),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(in) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('TABNAME','FLUXTAB')
          call READ_FLUX(OBJ,trim(SARG))
      END SELECT
      NARG=LR
      END SUBROUTINE SOURCE_INP_S

!---------------------------------------------------------
      SUBROUTINE SOURCE_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SOURCE),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('FLUX');   ARG(1)=OBJ%FLUX
        CASE('TEMP');   ARG(1)=OBJ%TEMP
        case('FRQ');  ARG(1)=OBJ%FRQ*1.D6
        CASE('DELAY');  ARG(1)=OBJ%DELAY/1.D3
        CASE('LAMW');   ARG(1)=OBJ%LAMW
        CASE('PULSW');  ARG(1)=OBJ%PULSW/1.D3
        case('TYPE');    ARG(1)=OBJ%TYP
        case('OVERLAP');   ARG(1)=LOG2INT(OBJ%OVERLAP)
        CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SOURCE_OUT_R

!---------------------------------------------------------
      SUBROUTINE SOURCE_OUT_S(OBJ,PNAME,SARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(SOURCE),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*), intent(out) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=0
      SELECT CASE (trim(PNAME))
      ! backward compatibility requires TABNAME ...
        CASE('TABNAME','FLUXTAB')
          if (len_trim(OBJ%FLUXTAB).gt.0) SARG=trim(OBJ%FLUXTAB)
          LR=1
      END SELECT
      NARG=LR
      END SUBROUTINE SOURCE_OUT_S

!==============================================================
! Component specific code
!==============================================================
!---------------------------------------------------------------
      SUBROUTINE READ_FLUX(OBJ,FNAME)
! read source lookup table
!---------------------------------------------------------------
      TYPE(SOURCE),intent(out) :: OBJ
      CHARACTER*(*),intent(in) :: FNAME
      integer :: ilin
      call READ_FLUX_TABLE(trim(FNAME),ilin)
      if (ilin.gt.0) call REPORT_FLUX_TABLE(ilin)
      OBJ%FLUXTAB=trim(FLUX_TABLE_NAMES)
      end SUBROUTINE READ_FLUX

      end MODULE SOURCES