!//////////////////////////////////////////////////////////////////////
!////  $Id: guides.f,v 1.32 2019/08/15 15:02:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2016, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.32 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: GUIDE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE GUIDES
      use CONSTANTS
      use CLASSES
      use FIELDDATA
      use FRAMES
      use MIRRLOG
      use RNDGEN
      implicit none

! GUIDE:
! TYP < 0        .. ignore
! TYP = 0        .. normal course or soller collimator
! TYP = 1        .. guide or bender (reflecting walls)
! TYP = 2        .. parabolic guide
! TYP = 3        .. parabolic guide with optimum lamella lengths
! TYP = 4        .. elliptic guide, parallel lamellae at the wider end
      TYPE  GUIDE; SEQUENCE
        TYPE(TFRAME) :: FRAME
        real(kind(1.D0)) :: W2,H2           ! width & height at the exit
        real(kind(1.D0)) :: GHLU,GHLD,GHRU,GHRD,GVT,GVB  ! horizontal & vertical (H&V) critical angles
        real(kind(1.D0)) :: RHLU,RHLD,RHRU,RHRD,RVT,RVB  ! horizontal& vertical reflectivities
        real(kind(1.D0)) :: DLH,DLV         ! thickness of H&V lamellae
        real(kind(1.D0)) :: ROH,ROV         ! curvatures
        real(kind(1.D0)) :: MU              ! absorption coeff. (1/mm) for lamellae
        !real(kind(1.D0)) :: LAMBDA          ! nominal wavelength (used to calculate mean critical angle)
        real(kind(1.D0)) :: WAV  ! waviness (sigma, [rad])
        real(kind(1.D0)) :: MISALIGN(2)   ! segments misalignment in mm (gaussian sigma)
        INTEGER :: NLH,NLV      ! number of slits (H&V), max.=127
        INTEGER :: TYP          ! (-1) ignored, normal(0), guide (1), parabolic (2,3), elliptic (4)
        INTEGER :: NHLU,NHLD,NHRU,NHRD,NVT,NVB ! indexes pointing to the field containg reflectivity data
        LOGICAL :: OSCILLATE     ! oscillating collimator, if > 0
        LOGICAL :: NODIR        ! no direct view (ballistic guides only)
        logical :: LAMONLY      ! only lamellae are used, no outer walls
        logical :: MONITOR      ! monitor neutrons captured by the wall
        LOGICAL :: CLOSED       ! if true, no neutrons can path through
        integer :: LOGBNC    ! log bounces at reflecting walls - produce map in BNCREC
        INTEGER :: ONESIDE      ! reflection only from inner sides of lamellae, if >0
        INTEGER :: TRFRONT      ! transparent front surface of lamellae, if >0
        INTEGER :: MATER        ! blade material: (0) absorbing (1) MU*lambda (2) Si (3) Al2O3
! calculated fields
        integer :: MLOG   ! reference to a mirror logger  (see MIRRLOG unit)
        real(kind(1.D0)) ::  CH,CV           ! calculated curvatures, parabolic, elliptic, etc guides
        REAL :: AH(0:255),AV(0:255) ! parameters of lamellae (for parabolic profile)
        REAL :: LH(0:255),LV(0:255) ! lengths of lamellae (for parabolic profile)
      END TYPE GUIDE

      integer, parameter :: GUIDES_DIM=128
! pointer type for GUIDE
      type PGUIDE
        integer :: IDX  ! instance index
        TYPE(GUIDE),pointer :: X
        ! record bounces for walls (left,right,top,bootom)
      end type PGUIDE
! instances of GUIDE. AGUIDES(0) is always unallocated
      integer :: GUIDES_NC
      type(PGUIDE) :: AGUIDES(1:GUIDES_DIM)

      logical :: GUIDES_dbg
      integer :: GUIDES_idbg

      contains

!---------------------------------------------------
      subroutine GUIDE_GET_MTAB(ID,OBJ,NR,MC)
! return mirror table index for given surface
! return 0 if no table is defined
! ID=1  left
! ID=2  right
! ID=3  top
! ID=4  bottom
! NOTE: current version recognizes only two tables
! All horizontal/vertical tables have the same m
! reflectivity is ignored (it is calculaed analytically later if needed)
!---------------------------------------------------
      INTEGER,intent(in) :: ID
      TYPE(GUIDE) :: OBJ
      integer,intent(out) :: NR
      real(kind(1.D0)),intent(out) :: MC
      select case(ID)
      case(1,2)
        NR=OBJ%NHLU
        MC=OBJ%GHLU/gammaNi
        !R0=OBJ%RHLU
      case(3,4)
        NR=OBJ%NVT
        MC=OBJ%GVT/gammaNi
        !R0=OBJ%RVT
      case DEFAULT
        NR=0
        MC=0.D0
        !R0=0.D0
      END select
      END subroutine GUIDE_GET_MTAB


!---------------------------------------------------------
      SUBROUTINE GUIDE_MISALIGN(OBJ)
! generate random misalignment for each segment
!---------------------------------------------------------
      TYPE(GUIDE) :: OBJ
      integer :: i
      do i=1,2
        if (OBJ%MISALIGN(i).ne.0.D0) then
          OBJ%FRAME%DM(i)=OBJ%MISALIGN(i)*GASDEV1(0.D0,3.D0)
        endif
      enddo
      END SUBROUTINE GUIDE_MISALIGN

!---------------------------------------------------------
      subroutine GUIDE_REGMLOG(INST)
! register logger for mirror in MIRRLOG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
1     FORMAT(a,1x,6(1x,G11.5))
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
      if (OBJ%LOGBNC>0) then
        OBJ%MLOG=MLOGS_REGISTER(OBJ%FRAME%ID,CCLS_GUIDE,INST,OBJ%LOGBNC,OBJ%FRAME%SIZE(3))
      else
        OBJ%MLOG=0
      endif
      end subroutine GUIDE_REGMLOG

!---------------------------------------------------------
      subroutine GUIDE_PROFILE(INST,IU)
! write guide profile in the given file unit
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
      integer,intent(in) :: IU
      integer :: i,NP
      real(kind(1.D0)) :: Z,DZ,X(4),LN,R(3),R1(3),R2(3),W,M(4)
      real(kind(1.D0)) :: GUIDE_LAM,PGUIDE_LAM,EGUIDE_LAM
      character(LEN_LINE) :: LINE,LINE2
      character(32) :: CNUM
1     FORMAT(a,1x,6(1x,G11.5))
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
      select case(OBJ%TYP)
        case(1)
          NP=64
        case(2,3)
          NP=64
        case(4)
          NP=64
        case DEFAULT
          NP=2
      END select
      LN=OBJ%FRAME%SIZE(3)
      M(1)=OBJ%GHLU/GammaNi
      M(2)=OBJ%GHRU/GammaNi
      M(3)=OBJ%GVT/GammaNi
      M(4)=OBJ%GVB/GammaNi
      ! include random misalignment if set
      call GUIDE_MISALIGN(OBJ)
      call ARRAY2STR(M,4,LINE2)
      DZ=LN/NP
      ! write(*,1) 'GUIDE_PROFILE typ=',OBJ%TYP
      do i=0,NP
        Z=i*DZ
        select case(OBJ%TYP)
        case(0)
          W=OBJ%FRAME%SIZE(1)
          X(1)=0.5D0*(W+(OBJ%W2-W)*Z/LN)
          X(2)=-0.5D0*(W+(OBJ%W2-W)*Z/LN)
          W=OBJ%FRAME%SIZE(2)
          X(3)=0.5D0*(W+(OBJ%H2-W)*Z/LN)
          X(4)=-0.5D0*(W+(OBJ%H2-W)*Z/LN)
        case(1)
          X(2)=GUIDE_LAM(OBJ,0,0,Z,1)
          X(1)=GUIDE_LAM(OBJ,0,OBJ%NLH,Z,1)
          X(4)=GUIDE_LAM(OBJ,0,0,Z,2)
          X(3)=GUIDE_LAM(OBJ,0,OBJ%NLV,Z,2)
        case(2,3)
          X(2)=PGUIDE_LAM(OBJ,0,0,Z,1)
          X(1)=PGUIDE_LAM(OBJ,0,OBJ%NLH,Z,1)
          X(4)=PGUIDE_LAM(OBJ,0,0,Z,2)
          X(3)=PGUIDE_LAM(OBJ,0,OBJ%NLV,Z,2)
        case(4)
          X(2)=EGUIDE_LAM(OBJ,0,0,Z,1)
          X(1)=EGUIDE_LAM(OBJ,0,OBJ%NLH,Z,1)
          X(4)=EGUIDE_LAM(OBJ,0,0,Z,2)
          X(3)=EGUIDE_LAM(OBJ,0,OBJ%NLV,Z,2)
        case DEFAULT
          X=0.D0
        END select
        ! write(*,1) 'GUIDE_PROFILE ', X(1:4)
! convert to lab. coordinates:
        R=(/X(1),X(3),Z/)
        call Loc2GlobR(OBJ%FRAME,R,R1)
        R=(/X(2),X(4),Z/)
        call Loc2GlobR(OBJ%FRAME,R,R2)
        X=(/R1(1),R2(1),R1(2),R2(2)/)
        R=(/0.D0,0.D0,Z/)
        call Loc2GlobR(OBJ%FRAME,R,R1)
        Z=R1(3)
        call ARRAY2STR(X,4,LINE)
        call FLOAT2STR(Z,CNUM)
        write(IU,*) i,' ',trim(CNUM),' ',trim(LINE),' ',trim(LINE2)
      enddo
      end subroutine GUIDE_PROFILE

!-------------------------------------------------------------------------
! Creator for GUIDE, return instance number
! Memory is allocated on the first free element of AGUIDES array
!-------------------------------------------------------------------------
      integer function GUIDE_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.GUIDES_DIM).and.(AGUIDES(i)%IDX.gt.0))
          i=i+1
        enddo
        if (AGUIDES(i)%IDX.le.0) then
          allocate(AGUIDES(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            GUIDES_NC=GUIDES_NC+1
            call GUIDE_DEFAULT(i)
            AGUIDES(i)%IDX=i
            AGUIDES(i)%X%FRAME%ID=trim(ID)
            AGUIDES(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          GUIDE_CREATE=i
        else
          call INT2STR(GUIDES_DIM,S)
          call MSG_ERROR('GUIDE_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          GUIDE_CREATE=0
        endif
      end function GUIDE_CREATE

!---------------------------------------------------------
      subroutine GUIDE_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.GUIDES_DIM)) then
          if (associated(AGUIDES(inst)%X)) then
            DEallocate(AGUIDES(inst)%X)
            nullify(AGUIDES(inst)%X)
            GUIDES_NC=GUIDES_NC-1
          endif
          AGUIDES(inst)%IDX=0
        endif
      end subroutine GUIDE_DISPOSE

!---------------------------------------------------------
      subroutine GUIDE_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,GUIDES_DIM
          call GUIDE_DISPOSE(i)
        enddo
      end subroutine GUIDE_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddGUIDE(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (GUIDE_isValid(INST)) then
        i=GUIDE_CREATE(AGUIDES(INST)%X%FRAME%ID,AGUIDES(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          AGUIDES(i)%X=AGUIDES(INST)%X
          AGUIDES(i)%IDX=i
        !  write(*,*) 'AddGUIDE ',trim(AGUIDES(i)%X%FRAME%ID),i
        endif
      endif
      AddGUIDE=i
      end function AddGUIDE

!------------------------------------
      SUBROUTINE GUIDE_PREPARE(INST,IERR)
! prepare dependent fields after input
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      integer :: SS
      REAL(KIND(1.D0)) :: ROPT
      REAL(KIND(1.D0)) :: V(3), DV(3), R(3,3), ALPHA, BETA, DX, DY
      TYPE(GUIDE),POINTER :: OBJ
1     format(a,': ',7(1x,G12.6))
      IERR=0
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
! adjust curvature of a curved guide if requiered
      if ((OBJ%TYP == 1).and.(OBJ%NODIR).and.(OBJ%FRAME%SIZE(3)>0.D0)) then
        if (abs(OBJ%ROH)>0.D0) then
          SS=NINT(sign(1.D0,OBJ%ROH))
          ROPT=4.D0*(OBJ%FRAME%SIZE(1)+OBJ%W2)/OBJ%FRAME%SIZE(3)**2
          if (abs(OBJ%ROH)<ROPT) OBJ%ROH=SS*ROPT
        else if (abs(OBJ%ROV)>0.D0) then
          SS=NINT(sign(1.D0,OBJ%ROV))
          ROPT=4.D0*(OBJ%FRAME%SIZE(2)+OBJ%H2)/OBJ%FRAME%SIZE(3)**2
          if (abs(OBJ%ROV)<ROPT) OBJ%ROV=SS*ROPT
        endif
      endif
      if (OBJ%TYP<=0) OBJ%LOGBNC=0
      ! translucent material is allowed only for straight guides
      if (OBJ%TYP.ne.1) OBJ%MATER=0
      if (OBJ%MONITOR) then
        write(*,1) trim(OBJ%FRAME%ID)//' is monitor'
      endif
! initialize transformation matrices
      call FRAME_INIT_MAT(OBJ%FRAME)
! for curved guides, add deflection
      OBJ%FRAME%ISDEFL=(ABS(OBJ%ROH)+ABS(OBJ%ROV)>1.D-9)
      if (OBJ%FRAME%ISDEFL) call FRAME_DEFLECTION(OBJ%FRAME, (/OBJ%ROH, OBJ%ROV/))
      END SUBROUTINE GUIDE_PREPARE

!-------------------------------------------------------------
      logical function GUIDE_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      GUIDE_isValid= ((INST.gt.0).and.(INST.le.GUIDES_DIM).and.associated(AGUIDES(INST)%X))
      end function GUIDE_isValid

!-------------------------------------------------------------
      SUBROUTINE GUIDE_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PGUIDE),intent(out) :: OBJ
        if (GUIDE_isValid(INST)) then
          OBJ%IDX=AGUIDES(INST)%IDX
          OBJ%X => AGUIDES(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE GUIDE_GET

!-------------------------------------------------------------
      SUBROUTINE GUIDE_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (GUIDE_isValid(INST)) then
          OBJ%IDX=AGUIDES(INST)%IDX
          OBJ%X => AGUIDES(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE GUIDE_GET_FRAME


!------------------------------------
      SUBROUTINE GUIDE_DEFAULT(INST)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
        if (.not.GUIDE_isValid(INST)) return
        OBJ => AGUIDES(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%SIZE(3)=1.D0
        OBJ%FRAME%CLASS=CCLS_GUIDE
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%TYP=0
        OBJ%W2=100.D0
        OBJ%H2=100.D0
        OBJ%ROH=0.D0
        OBJ%ROV=0.D0
        OBJ%NLH=1
        OBJ%NLV=1
        OBJ%DLH=8.D-2
        OBJ%DLV=8.D-2
        OBJ%MU=1.D9
        OBJ%WAV=2.D-4
        OBJ%OSCILLATE=.FALSE.
        OBJ%NODIR=.false.
        OBJ%LAMONLY=.false.
        OBJ%MONITOR=.false.
        OBJ%CLOSED=.false.
        OBJ%LOGBNC=0
        OBJ%MLOG=0
        OBJ%MATER=0
        OBJ%MISALIGN=(/0.D0,0.D0/)
      END SUBROUTINE GUIDE_DEFAULT

!---------------------------------------------------------
      SUBROUTINE GUIDE_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call GUIDE_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'GUIDE_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE GUIDE_INP


!---------------------------------------------------------
      SUBROUTINE GUIDE_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(GUIDE),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.GUIDE_isValid(INST)) return
      OBJ => AGUIDES(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call GUIDE_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'GUIDE_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE GUIDE_OUT

!---------------------------------------------------------
      SUBROUTINE GUIDE_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(GUIDE),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! type
          CASE('TYPE')
  !          write(*,*) 'GUIDE_INP ARG=',NINT(ARG(1))
            OBJ%TYP=MIN(NINT(ARG(1)),4)
            OBJ%TYP=MAX(OBJ%TYP,-1)
! exit size (H,V) [mm]
          CASE('EXIT');    OBJ%W2=ARG(1);OBJ%H2=ARG(1+1);LR=2
          CASE('EXITH');   OBJ%W2=ARG(1)
          CASE('EXITV');   OBJ%H2=ARG(1)
! curvatures [1/m] - convert to [1/mm]
          CASE('RO');  OBJ%ROH=ARG(1)/1000.;OBJ%ROV=ARG(1+1)/1000.;LR=2
          CASE('ROH'); OBJ%ROH=ARG(1)/1000.
          CASE('ROV'); OBJ%ROV=ARG(1)/1000.
! number of lamellae [mm]
          CASE('N');  OBJ%NLH=NINT(ARG(1));OBJ%NLV=NINT(ARG(1+1));LR=2
          CASE('NH');  OBJ%NLH=NINT(ARG(1))
          CASE('NV');  OBJ%NLV=NINT(ARG(1))
! lamellae thickness [mm]
          CASE('DL');   OBJ%DLH=ARG(1);OBJ%DLV=ARG(1+1);LR=2
          CASE('DLH');  OBJ%DLH=ARG(1)
          CASE('DLV');  OBJ%DLV=ARG(1)
! critical angles [Ni nat.]
          CASE('M'); call BENDER_SETM(OBJ,ARG(1),1); call BENDER_SETM(OBJ,ARG(1+1),2); LR=2
          CASE('MH'); call BENDER_SETM(OBJ,ARG(1),1)
          CASE('MV'); call BENDER_SETM(OBJ,ARG(1),2)
! reflectivity
          CASE('REF'); call BENDER_SETREF(OBJ,ARG(1),1); call BENDER_SETREF(OBJ,ARG(1+1),2); LR=2
          CASE('REFH'); call BENDER_SETREF(OBJ,ARG(1),1)
          CASE('REFV'); call BENDER_SETREF(OBJ,ARG(1),2)
! absorption of lamellae [1/cm] - convert to [1/mm]
          CASE('MU'); OBJ%MU=ARG(1)/10.
! material index
          CASE('MATER'); OBJ%MATER=NINT(ARG(1))
! misalignment [mm]
          CASE('MISALIGN');  OBJ%MISALIGN(1)=ARG(1);OBJ%MISALIGN(2)=ARG(2);LR=2
! waviness [mrad]
          CASE('WAV'); OBJ%WAV=ARG(1)/1.D3
! onesided
          CASE('ONESIDE');  OBJ%ONESIDE=NINT(ARG(1))
! transmission front
          CASE('TRFRONT');  OBJ%TRFRONT=NINT(ARG(1))
! only lamnellae
          CASE('LAMONLY');  OBJ%LAMONLY=(NINT(ARG(1)).eq.1)
! monitor (bit mask)
          case('MONITOR'); OBJ%MONITOR=(NINT(ARG(1)).eq.1)
! closed guide
          case('CLOSED'); OBJ%CLOSED=(NINT(ARG(1)).eq.1)
! log bounces
          case('LOGBNC'); OBJ%LOGBNC=MAX(MIN(NINT(ARG(1)),7),0)
! automatically adjust in reflection position
          CASE('OSC');   OBJ%OSCILLATE=(NINT(ARG(1)).eq.1)
! no direct view (ballistic guides only)
          CASE('NODIR');   OBJ%NODIR=(NINT(ARG(1)).eq.1)
! try FRAME parameters
          CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE GUIDE_INP_R

!---------------------------------------------------------
      SUBROUTINE GUIDE_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(GUIDE),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! type
          CASE('TYPE')
            ARG(1)=OBJ%TYP ! ;write(*,*) 'GUIDE_OUT OBJ%TYP=',OBJ%TYP
! exit size (H,V) [mm]
          CASE('EXIT');    ARG(1)=OBJ%W2;ARG(1+1)=OBJ%H2;LR=2
          CASE('EXITH');   ARG(1)=OBJ%W2
          CASE('EXITV');   ARG(1)=OBJ%H2
! curvatures [1/m] - convert to [1/mm]
          CASE('RO');  ARG(1)=OBJ%ROH*1000.;ARG(1+1)=OBJ%ROV*1000.;LR=2
          CASE('ROH'); ARG(1)=OBJ%ROH*1000.
          CASE('ROV'); ARG(1)=OBJ%ROV*1000.
! number of lamellae [mm]
          CASE('N');   ARG(1)=OBJ%NLH;ARG(1+1)=OBJ%NLV;LR=2
          CASE('NH');  ARG(1)=OBJ%NLH
          CASE('NV');  ARG(1)=OBJ%NLV
! lamellae thickness [mm]
          CASE('DL');   ARG(1)=OBJ%DLH;ARG(1+1)=OBJ%DLV;LR=2
          CASE('DLH');  ARG(1)=OBJ%DLH
          CASE('DLV');  ARG(1)=OBJ%DLV
! critical angles [Ni nat.]
          CASE('M');  ARG(1)=OBJ%GHLU/GammaNi; ARG(1+1)=OBJ%GVT/GammaNi; LR=2
          CASE('MH'); ARG(1)=OBJ%GHLU/GammaNi
          CASE('MV'); ARG(1)=OBJ%GVT/GammaNi
! reflectivity
          CASE('REF');  ARG(1)=OBJ%RHLU; ARG(1+1)=OBJ%RVT; LR=2
          CASE('REFH'); ARG(1)=OBJ%RHLU
          CASE('REFV'); ARG(1)=OBJ%RVT
! absorption of lamellae [1/cm] - convert to [1/mm]
          CASE('MU');   ARG(1)=OBJ%MU*10.
! material index
          CASE('MATER');   ARG(1)=OBJ%MATER
! misalignment [mm]
          CASE('MISALIGN');    ARG(1)=OBJ%MISALIGN(1);ARG(2)=OBJ%MISALIGN(2);LR=2
! waviness [mrad]
          CASE('WAV'); ARG(1)=OBJ%WAV*1.D3
! onesided
          CASE('ONESIDE');  ARG(1)=OBJ%ONESIDE
! transmission front
          CASE('TRFRONT');  ARG(1)=OBJ%TRFRONT
! only lamnellae
          CASE('LAMONLY');  ARG(1)=LOG2INT(OBJ%LAMONLY)
! monitor (bit mask)
          case('MONITOR'); ARG(1)=LOG2INT(OBJ%MONITOR)
! closed guide
          case('CLOSED'); ARG(1)=LOG2INT(OBJ%CLOSED)
! log bounces
          case('LOGBNC'); ARG(1)=OBJ%LOGBNC
! automatically adjust in reflection position
          CASE('OSC');   ARG(1)=LOG2INT(OBJ%OSCILLATE)
! no direct view (ballistic guides only)
          CASE('NODIR');   ARG(1)=LOG2INT(OBJ%NODIR)
! try FRAME parameters
          CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE GUIDE_OUT_R

      end MODULE GUIDES