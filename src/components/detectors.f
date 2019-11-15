!//////////////////////////////////////////////////////////////////////
!////  $Id: detectors.f,v 1.17 2015/04/13 18:38:51 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.17 $
!////     $Date: 2015/04/13 18:38:51 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: DETECTOR
!////
!////////////////////////////////////////////////////////////////////////
      MODULE DETECTORS
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      implicit none
      save
      private

      integer,parameter,public :: det_type_monit=0
      integer,parameter,public :: det_type_array=1
      integer,parameter,public :: det_type_psd=2
      integer,parameter,public :: det_type_cyl=3
      integer,parameter,public :: det_type_sph=4

! DETECTOR
! TYP = 0        .. area, inf. thin, 100% efficiency, 0 mm resolution
! TYP = 1        .. tube array of ND elements, efficiency=ALPHA
! TYP = 2        .. PSD with spatial resolution = RESX,RESY
! TYP = 3        .. curved cylindric
! TYP = 4        .. curved spherical
      TYPE  DETECTOR; SEQUENCE
        TYPE(TFRAME) :: FRAME
        INTEGER :: TYP ! detector type
        REAL(kind(1.D0)) :: RESX,RESY ! spatial resolution [mm], gaussian
        REAL(kind(1.D0)) :: ALPHA     ! efficiency coeff. [1/A/cm]
        REAL(kind(1.D0)) :: RAD       ! radius for curved layout [mm]
        REAL(kind(1.D0)) :: THMIN     ! lower angular limit for data analysis [deg]
        REAL(kind(1.D0)) :: THMAX     ! upper angular limit for data analysis [deg]
        REAL(kind(1.D0)) :: SPACE     ! distance between adjacent tubes [mm]
        INTEGER :: ND                 ! number of tubes in an array
      ! calculated fields during ADJUST
        REAL(kind(1.D0)) :: COSX,SINX,COSY,SINY ! cos and sin of the sector angles defined by the left/right and top/bottom walls
        REAL(kind(1.D0)) :: KI(3), KI0 ! nominal KI vector in detector local coordinates
        REAL(kind(1.D0)) :: TOF        ! nominal time of flight [m/hbar units]
        REAL(kind(1.D0)) :: LAMBDA     ! nominal wavelength
      END TYPE DETECTOR

      integer, parameter :: DETECTORS_DIM=128
! pointer type for DETECTOR
      type PDETECTOR
        integer :: IDX  ! instance index
        TYPE(DETECTOR),pointer :: X
      end type PDETECTOR
! instances of DETECTOR. ADETECTORS(0) is always unallocated
      integer :: DETECTORS_NC
      type(PDETECTOR) :: ADETECTORS(1:DETECTORS_DIM)

      public DETECTOR,PDETECTOR,ADETECTORS,DETECTOR_GET,DETECTOR_PREPARE,DETECTORS_DIM
      public DETECTOR_DEFAULT,DETECTOR_INP,DETECTOR_OUT,DETECTOR_GET_FRAME,DETECTOR_isValid
      public DETECTOR_CREATE,DETECTOR_DISPOSE,DETECTOR_DISPOSE_ALL,AddDETECTOR

      contains

!-------------------------------------------------------------------------
! Creator for DETECTOR, return instance number
! Memory is allocated on the first free element of ADETECTORS array
!-------------------------------------------------------------------------
      integer function DETECTOR_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.DETECTORS_DIM).and.(ADETECTORS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (ADETECTORS(i)%IDX.le.0) then
          allocate(ADETECTORS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            DETECTORS_NC=DETECTORS_NC+1
            call DETECTOR_DEFAULT(i)
            ADETECTORS(i)%IDX=i
            ADETECTORS(i)%X%FRAME%ID=trim(ID)
            ADETECTORS(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          DETECTOR_CREATE=i
        else
          call INT2STR(DETECTORS_DIM,S)
          call MSG_ERROR('DETECTOR_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          DETECTOR_CREATE=0
        endif
      end function DETECTOR_CREATE

!---------------------------------------------------------
      subroutine DETECTOR_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.DETECTORS_DIM)) then
          if (associated(ADETECTORS(inst)%X)) then
            DEallocate(ADETECTORS(inst)%X)
            nullify(ADETECTORS(inst)%X)
            DETECTORS_NC=DETECTORS_NC-1
          endif
          ADETECTORS(inst)%IDX=0
        endif
      end subroutine DETECTOR_DISPOSE

!---------------------------------------------------------
      subroutine DETECTOR_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,DETECTORS_DIM
          call DETECTOR_DISPOSE(i)
        enddo
      end subroutine DETECTOR_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddDETECTOR(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (DETECTOR_isValid(INST)) then
        i=DETECTOR_CREATE(ADETECTORS(INST)%X%FRAME%ID,ADETECTORS(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          ADETECTORS(i)%X=ADETECTORS(INST)%X
          ADETECTORS(i)%IDX=i
        !  write(*,*) 'AddDETECTOR ',trim(ADETECTORS(i)%X%FRAME%ID),i
        endif
      endif
      AddDETECTOR=i
      end function AddDETECTOR

!------------------------------------
      SUBROUTINE DETECTOR_PREPARE(INST,IERR)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(DETECTOR),POINTER :: OBJ
        IERR=0
        if (.not.DETECTOR_isValid(INST)) return
        OBJ => ADETECTORS(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE DETECTOR_PREPARE


!-------------------------------------------------------------
      logical function DETECTOR_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      DETECTOR_isValid= ((INST.gt.0).and.(INST.le.DETECTORS_DIM).and.associated(ADETECTORS(INST)%X))
      end function DETECTOR_isValid


!-------------------------------------------------------------
      SUBROUTINE DETECTOR_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PDETECTOR) :: OBJ
        if (DETECTOR_isValid(INST)) then
          OBJ%IDX=ADETECTORS(INST)%IDX
          OBJ%X => ADETECTORS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE DETECTOR_GET

!-------------------------------------------------------------
      SUBROUTINE DETECTOR_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (DETECTOR_isValid(INST)) then
          OBJ%IDX=ADETECTORS(INST)%IDX
          OBJ%X => ADETECTORS(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE DETECTOR_GET_FRAME

!------------------------------------
      SUBROUTINE DETECTOR_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(DETECTOR),POINTER :: OBJ
        if (.not.DETECTOR_isValid(INST)) return
        OBJ => ADETECTORS(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%CLASS=CCLS_DETECTOR
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%FRAME%SIZE=(/64.D0,64.D0,1.D0/)
        OBJ%ND=1
        OBJ%SPACE=0.D0
        OBJ%ALPHA=5.D-1
        OBJ%RESX=1.D0
        OBJ%RESY=1.D0
        OBJ%RAD=0.D0
        OBJ%THMIN=0.D0
        OBJ%THMAX=180.D0*deg
      END SUBROUTINE DETECTOR_DEFAULT


!---------------------------------------------------------
      SUBROUTINE DETECTOR_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DETECTOR),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.DETECTOR_isValid(INST)) return
      OBJ => ADETECTORS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call DETECTOR_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'DETECTOR_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE DETECTOR_INP

!---------------------------------------------------------
      SUBROUTINE DETECTOR_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(DETECTOR),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.DETECTOR_isValid(INST)) return
      OBJ => ADETECTORS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call DETECTOR_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'DETECTOR_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE DETECTOR_OUT

!---------------------------------------------------------
      SUBROUTINE DETECTOR_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DETECTOR),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! type
        CASE('TYPE')
          OBJ%TYP=MAX(0,MIN(NINT(ARG(1)),4))
! radius [mm]
        CASE('RAD');  OBJ%RAD=ARG(1)
! resolution [mm]
        CASE('RES');   OBJ%RESX=ARG(1);OBJ%RESY=ARG(2);LR=2
        CASE('RESX');  OBJ%RESX=ARG(1)
        CASE('RESY');  OBJ%RESY=ARG(1)
! efficiency coeff. [1/A/cm]
        CASE('ALPHA');   OBJ%ALPHA=ARG(1)/10.D0
! number of segments
        CASE('ND');  OBJ%ND=max(1,min(64,NINT(ARG(1))))
! spacing between segments [mm]
        CASE('SPACE'); OBJ%SPACE=ARG(1)
! lower angular limit for data analysis [deg]
        CASE('THMIN');   OBJ%THMIN=ARG(1)*deg
! upper angular limit for data analysis [deg]
        CASE('THMAX');   OBJ%THMAX=ARG(1)*deg
! try FRAME parameters
        CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE DETECTOR_INP_R

!---------------------------------------------------------
      SUBROUTINE DETECTOR_OUT_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(DETECTOR),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! type
          CASE('TYPE'); ARG(1)=OBJ%TYP
! resolution [mm]
          CASE('RES');   ARG(1)=OBJ%RESX;ARG(2)=OBJ%RESY;LR=2
          CASE('RESX');  ARG(1)=OBJ%RESX
          CASE('RESY');  ARG(1)=OBJ%RESY
! radius [mm]
          CASE('RAD');  ARG(1)=OBJ%RAD
! efficiency coeff. [1/A/cm]
          CASE('ALPHA');   ARG(1)=OBJ%ALPHA*10.D0
! number of segments
          CASE('ND');  ARG(1)=OBJ%ND
! spacing between segments [mm]
          CASE('SPACE'); ARG(1)=OBJ%SPACE
! lower angular limit for data analysis [deg]
          CASE('THMIN');   ARG(1)=OBJ%THMIN/deg
! upper angular limit for data analysis [deg]
          CASE('THMAX');   ARG(1)=OBJ%THMAX/deg
! try FRAME parameters
          CASE DEFAULT
            CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE DETECTOR_OUT_R

      end MODULE DETECTORS