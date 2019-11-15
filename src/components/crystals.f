!//////////////////////////////////////////////////////////////////////
!////  $Id: crystals.f,v 1.51 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.51 $
!////     $Date: 2019/08/15 17:24:08 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: CRYSTAL
!////
!////////////////////////////////////////////////////////////////////////
      MODULE CRYSTALS
      use CONSTANTS
      use CLASSES
!     use FIELDS
      use FIELDDATA
      use FRAMES
      ! use TRACINGDATA
      use TABLES
      use LATTICE
      USE TABLE_CRYSTALS
      use ARRAY3D
      implicit none
      SAVE

      integer,PARAMETER :: ctyp_filter=0
      integer,PARAMETER :: ctyp_mosaic=1
      integer,PARAMETER :: ctyp_bent=2
      integer,PARAMETER :: ctyp_gradient=3
      integer,PARAMETER :: cmos_gauss=0
      integer,PARAMETER :: cmos_lorentz=1
      integer,PARAMETER :: cmos_voigt=2
      integer,PARAMETER :: cmos_box=3
      INTEGER :: MDIST
      REAL(KIND(1.D0)) ::RDIST

! CRYSTAL
      TYPE CRYSTAL
        TYPE(TFRAME) :: FRAME
        character*8 REFNAME
        integer*4 TYP        ! (-1) ignored, (0) mosaic, (1) gradient, (2) mosaic + gradient
        integer*4 MODEL ! tracing method: simple (0), random walk (1)
        integer*4 SGN   ! sign of the take-off angle
        REAL*8 DHKL ! dhkl
        REAL*8 CHI,PSI ! cutting angles (H,V)
        REAL*8 POI,ANIZ ! dhkl,cutting angle, Poisson const., mos. anisotropy
        REAL*8 HMOS,VMOS,QML,VOL,MI,DW ! mosaicity (H,V), M.L. reflectivity, cell volume, absorption, Debye-Waller
        REAL*8 DEPTH !  thickness along K=(0,0,1) through R=(0,0,0)
        REAL*8 P0    !  reflection probability along K=(0,0,1) through R=(0,0,0) (without absorption)
        REAL*8 MAG !  magnetization
        REAL*8 RH,RV,RB,DH,DV,DB ! curvatures, segemnts spacing (H,V,z)
        REAL*8 SIGMAB, SIGMAA, A, THETAD,C2 ! Freund parameters for absorption (see crystals.lib)
        REAL*8 DLAM,DEXT,EXT1,DELTA         ! extinction parameters: lamellae thickness, ext. length, primary ext., Darwin box
        REAL*8 LAMBDA,QHKL,REF ! wavelength, kin. refl., ...
        REAL*8 DG_DR(3,3),G(3),GAMA(3),GTOT,DGR,DGA ! deformation gradient, diff. vector, d-gradient, gradient orientation
        LOGICAL*4 MAPG(3)
        INTEGER*4 IRNDV ! indices to random numbers array (vertical mosaicity)
        INTEGER*4 NH,NV,NB    ! number of segments (H,V,z)
        LOGICAL FOCH,FOCV     ! automatic focusing flags
        REAL*8 FH1,FH2,FV1,FV2 ! focal distances
        REAL*8 ASTACK ! stacking angle and its tangent
        LOGICAL :: STACKH,STACKV    ! smooth segment stacking on a curved surface
        LOGICAL :: FCONE  ! force flat-cone orientation
        logical :: AUTOADJ   ! automatically adjust in reflection position
        TYPE(TARRAY3D) :: A3D  ! used for tracing in 3D segments array
      END TYPE CRYSTAL

      integer, parameter :: CRYSTALS_DIM=128
! pointer type for CRYSTAL
      type PCRYSTAL
        integer :: IDX  ! instance index
        TYPE(CRYSTAL),pointer :: X
      end type PCRYSTAL
! instances of CRYSTAL.
      integer :: CRYSTALS_NC
      type(PCRYSTAL) :: ACRYSTALS(1:CRYSTALS_DIM)

      contains

!-------------------------------------------------------------------------
! Creator for CRYSTAL, return instance number
! Memory is allocated on the first free element of ACRYSTALS array
!-------------------------------------------------------------------------
      integer function CRYSTAL_CREATE(ID,NAMESTR)
      character(*) :: ID,NAMESTR
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i.lt.CRYSTALS_DIM).and.(ACRYSTALS(i)%IDX.gt.0))
          i=i+1
        enddo
        if (ACRYSTALS(i)%IDX.le.0) then
          allocate(ACRYSTALS(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            CRYSTALS_NC=CRYSTALS_NC+1
            call CRYSTAL_DEFAULT(i)
            ACRYSTALS(i)%IDX=i
            ACRYSTALS(i)%X%FRAME%ID=trim(ID)
            ACRYSTALS(i)%X%FRAME%NAME=trim(NAMESTR)
          endif
        endif
        if (ierr.eq.0) then
          CRYSTAL_CREATE=i
        else
          call INT2STR(CRYSTALS_DIM,S)
          call MSG_ERROR('CRYSTAL_CREATE','Reached maximum number or instances ('//trim(S)//')',1,1)
          CRYSTAL_CREATE=0
        endif
      end function CRYSTAL_CREATE

!---------------------------------------------------------
      subroutine CRYSTAL_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((inst.gt.0).and.(inst.le.CRYSTALS_DIM)) then
          if (associated(ACRYSTALS(inst)%X)) then
            DEallocate(ACRYSTALS(inst)%X)
            nullify(ACRYSTALS(inst)%X)
            CRYSTALS_NC=CRYSTALS_NC-1
          endif
          ACRYSTALS(inst)%IDX=0
        endif
      end subroutine CRYSTAL_DISPOSE

!---------------------------------------------------------
      subroutine CRYSTAL_DISPOSE_ALL
!---------------------------------------------------------
        integer :: i
        do i=1,CRYSTALS_DIM
          call CRYSTAL_DISPOSE(i)
        enddo
      end subroutine CRYSTAL_DISPOSE_ALL

!-------------------------------------------------------------------------
! Add a component
!-------------------------------------------------------------------------
      INTEGER function AddCRYSTAL(INST)
      integer,intent(in) :: INST
      integer :: i
      i=0
      if (CRYSTAL_isValid(INST)) then
        i=CRYSTAL_CREATE(ACRYSTALS(INST)%X%FRAME%ID,ACRYSTALS(INST)%X%FRAME%NAME)
        if (i.gt.0) then
          ACRYSTALS(i)%X=ACRYSTALS(INST)%X
          ACRYSTALS(i)%IDX=i
        !  write(*,*) 'AddCRYSTAL ',trim(ACRYSTALS(i)%X%FRAME%ID),i
        endif
      endif
      AddCRYSTAL=i
      end function AddCRYSTAL

!------------------------------------
      SUBROUTINE CRYSTAL_PREPARE(INST,IERR)
! set default parameters
!------------------------------------
      integer,intent(in) :: INST
      integer,intent(out) :: IERR
      TYPE(CRYSTAL),POINTER :: OBJ
        IERR=0
        if (.not.CRYSTAL_isValid(INST)) return
        OBJ => ACRYSTALS(INST)%X
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE CRYSTAL_PREPARE

!-------------------------------------------------------------
      logical function CRYSTAL_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      CRYSTAL_isValid= ((INST.gt.0).and.(INST.le.CRYSTALS_DIM).and.associated(ACRYSTALS(INST)%X))
      end function CRYSTAL_isValid


!-------------------------------------------------------------
      SUBROUTINE CRYSTAL_GET(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PCRYSTAL),intent(out) :: OBJ
        if (CRYSTAL_isValid(INST)) then
          OBJ%IDX=ACRYSTALS(INST)%IDX
          OBJ%X => ACRYSTALS(INST)%X
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE CRYSTAL_GET

!-------------------------------------------------------------
      SUBROUTINE CRYSTAL_GET_FRAME(INST,OBJ)
! return reference to given instance
!-------------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(PFRAME),intent(out) :: OBJ
        if (CRYSTAL_isValid(INST)) then
          OBJ%IDX=ACRYSTALS(INST)%IDX
          OBJ%X => ACRYSTALS(INST)%X%FRAME
        else
          NULLIFY(OBJ%X)
          OBJ%IDX=0
        endif
      end SUBROUTINE CRYSTAL_GET_FRAME

!------------------------------------
      SUBROUTINE CRYSTAL_DEFAULT(INST)
!------------------------------------
      integer,intent(in) :: INST
      TYPE(CRYSTAL),POINTER :: OBJ
        if (.not.CRYSTAL_isValid(INST)) return
        OBJ => ACRYSTALS(INST)%X
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%CLASS=CCLS_CRYSTAL
        OBJ%FRAME%SHAPE=FRAME_SHAPE_BOX
        OBJ%FRAME%SIZE=(/100.D0,50.D0,5.D-1/)
        OBJ%SGN=1
        OBJ%CHI=0.D0*deg
        OBJ%PSI=0.D0*deg
        OBJ%DLAM=10.D0
        OBJ%HMOS=20.D0*minute
        OBJ%VMOS=20.D0*minute
        OBJ%ANIZ=1.D0
        OBJ%POI=0.3D0
        OBJ%NH=1
        OBJ%NV=1
        OBJ%NB=1
        OBJ%MODEL=1
        OBJ%FOCH=.false.
        OBJ%FOCV=.false.
        OBJ%FH1=0.D0
        OBJ%FH2=0.D0
        OBJ%FV1=0.D0
        OBJ%FV2=0.D0
        OBJ%FCONE=.false.
        OBJ%MAG=0.D0
        OBJ%ASTACK=0.D0
        OBJ%STACKH=.false.
        OBJ%STACKV=.false.
        call SET_CRYST_PARAM(OBJ,'Ge 111',.true.)
      END SUBROUTINE CRYSTAL_DEFAULT

!-----------------------------------------
      SUBROUTINE LIST_REFLECTIONS(OBJ)
!-----------------------------------------
      TYPE(CRYSTAL) :: OBJ
      TYPE(TCR_ITEM) :: CR
      character*(2) :: ID
      integer :: HKL(3),IX
      real(kind(1.D0)) :: Q,D
      character*(128) :: STR
1     format(a,a,2(1x,G11.5))
      ID=OBJ%REFNAME(1:2)
      IX=GET_CR_INDEX(ID)
      write(*,*) 'LIST_REFLECTIONS ',ID,IX
      if (IX.gt.0) then
        call GET_CR_ITEM(IX,CR)
        HKL=(/1,1,1/)
        call IARRAY2STR(HKL,3,STR)
        Q=GET_QML(CR%LAT,HKL,298.D0)
        D=GET_DHKL(CR%LAT,HKL,298.D0)
        write(*,1) ID,trim(STR),D,Q
        HKL=(/2,2,0/)
        call IARRAY2STR(HKL,3,STR)
        Q=GET_QML(CR%LAT,HKL,298.D0)
        D=GET_DHKL(CR%LAT,HKL,298.D0)
        write(*,1) ID,trim(STR),D,Q
        HKL=(/1,0,0/)
        call IARRAY2STR(HKL,3,STR)
        Q=GET_QML(CR%LAT,HKL,298.D0)
        D=GET_DHKL(CR%LAT,HKL,298.D0)
        write(*,1) ID,trim(STR),D,Q
      endif

      end SUBROUTINE LIST_REFLECTIONS


!------------------------------------------------------------
      SUBROUTINE CRYSTAL_AUTOFOCUS(CR)
! Adjust crystal focal distances
!------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      real(kind(1.D0)) :: K0,STH,LL,TH,TH1,TH2
      real(kind(1.D0)) :: DS(3),SZ(3),R(3)
      integer :: NS(3),i
      logical :: BENTH,BENTV
        K0=2.D0*PI/CR%LAMBDA
        STH=PI/(K0*CR%DHKL) ! sin(thetaB)
        if (ABS(STH).ge.1.D0) return
        TH=ASIN(STH)
        TH1=SIN(TH-CR%SGN*CR%CHI) ! angle on the side of FH1
        TH2=SIN(TH+CR%SGN*CR%CHI) ! angle on the side of FH2
        if (CR%FOCH) THEN
          LL=CR%FH1*CR%FH2
          if (LL.gt.0.D0) CR%RH=(CR%FH1*TH2+CR%FH2*TH1)/LL*0.5D-3
        endif
        if (CR%FOCV) THEN
          LL=CR%FV1*CR%FV2*COS(CR%CHI)
          if (LL.gt.0.D0) CR%RV=(CR%FV1+CR%FV2)/LL/STH*0.5D-3
        endif
    ! initialize A3D data
      ! get gap vector
        DS=(/CR%DH,CR%DV,CR%DB/)
      ! get segment number vector
        NS=(/CR%NH,CR%NV,CR%NB/)
      ! get curvatures vector
        R=(/CR%RH,CR%RV,CR%RB/)
      ! slab size
        do i=1,3
          SZ(i)=CR%FRAME%SIZE(i)/NS(i)-DS(i)
        enddo
        BENTH=(R(1).ne.0.D0)
        BENTV=.false.
        call ARRAY3D_INIT(CR%A3D,SZ,DS,NS,R,tan(CR%ASTACK),CR%STACKH,CR%STACKV,BENTH,BENTV)
      end SUBROUTINE CRYSTAL_AUTOFOCUS

!------------------------------------------------------------
      SUBROUTINE CRYSTAL_ORIENT(CR,DTHETA,ierr)
! Adjust crystal in the position with G || (Kf-Ki)
! Lambda and Kf axis (REXI) must already be defined, Ki || (0 0 1)
! DTHETA is additional rocking angle (misalignment in GON(3))
! NOTE: Scattering angle (THETA) is defined by the direction of Kf axis
!       Therefore, this procedure ensures Bragg position
!       only if DTHETA=0 and LAMBDA=2*DHKL*sin(THETA/2)
!------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      integer,intent(out) :: IERR
      real(kind(1.D0)),intent(in) :: DTHETA
      real(kind(1.D0)) :: K0,KI(3),KF(3),Q(3),GON(3),THETA
      real(kind(1.D0)),parameter :: EPS=1.D-10
      ierr=1
      if (CR%TYP.eq.ctyp_filter) then
        call SetAxisAngles(CR%FRAME,(/0.D0,0.D0,0.D0/))
        call SetGonioAngles(CR%FRAME,(/0.D0,0.D0,0.D0/))
        ierr=0
        return
      else if ((CR%LAMBDA.gt.0.D0).and.(CR%LAMBDA.le.2.D0*CR%DHKL)) then
        THETA=2.D0*ASIN(CR%LAMBDA/2.D0/CR%DHKL)*CR%SGN
        if (CR%FCONE) then
          call SetAxisAngles(CR%FRAME,(/0.D0,-ABS(THETA),0.D0/))
        else
          call SetAxisAngles(CR%FRAME,(/THETA,0.D0,0.D0/))
        endif
        ierr=0
      else
        return
      endif
! adjust goniometer in reflection position + DTHETA
      K0=2.D0*PI/CR%LAMBDA
      KI=(/0.D0,0.D0,K0/)
      call M3xV3(-CR%FRAME%MEXI,CR%FRAME%REXI,KI,KF)
      Q=KF-KI
      call CalGonioAngles(Q,CR%G,.false.,GON)
      !write(*,*) 'CRYSTAL_ORIENT ',trim(CR%FRAME%ID),' GON=',GON(1:3)
      GON(3)=GON(3)+DTHETA
      call SetGonioAngles(CR%FRAME,GON)
      end SUBROUTINE CRYSTAL_ORIENT

!---------------------------------------------------------------
      SUBROUTINE SET_CRYST_PARAM(CR,CRNAME,SIL)
!---------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      CHARACTER*(*) CRNAME
      logical,intent(in) :: SIL
      TYPE(TCR_PARAM) :: PARAM
1     FORMAT(a,' set to ',a)
2     FORMAT(a,' is ',a)
      if (len_trim(CRNAME).gt.0) then
        if (FINDCRYSTPARAM(trim(CRNAME),PARAM)) then
          CR%REFNAME=PARAM%ID
          CR%dhkl=PARAM%dhkl
          CR%QML= PARAM%QML
          CR%sigmab=PARAM%sigmab
          CR%sigmaa=PARAM%sigmaa
          CR%VOL=PARAM%VOL
          CR%A=PARAM%A
          CR%thetaD=PARAM%thetaD
          CR%C2=PARAM%C2
          CR%poi=PARAM%poi
          if (.NOT.SIL) write(*,1) trim(CR%FRAME%NAME),trim(CR%REFNAME)
        else
          call MSG_WARN(trim(CRNAME)//' not found in crystal library',1)
        endif
      else
        write(*,2) trim(CR%FRAME%NAME),trim(CR%REFNAME)
      endif
      end SUBROUTINE SET_CRYST_PARAM

!-------------------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION CRYST_GET_LAMBDA(CR)
! get wavelength setting for a crystal
! calculated from vector Bragg condition using current crystal orientation
! and direct beam along (0,0,1) in incident coordinate system
!--------------------------------------------------------------------------
      TYPE (CRYSTAL) :: CR
      real(kind(1.D0)) :: K(3),VK(3),Z,K0
      K=(/0.D0,0.D0,1.D0/)
  ! transform (0,0,1) to local coordinates
      CALL M3XV3(CR%FRAME%MLOC,CR%FRAME%RLOC,K,VK)
      CALL V3AV3(-1,VK,CR%FRAME%VEL,K) ! subtract m/h*velocity
      Z=2.D0*v3xv3(K,CR%G)
      if (abs(Z).ge.1.D-10) then
        K0=-CR%GTOT**2/Z
        CRYST_GET_LAMBDA=2.D0*PI/K0
      else
        CRYST_GET_LAMBDA=0.D0
      endif
      end FUNCTION CRYST_GET_LAMBDA

!---------------------------------------------------------
      SUBROUTINE CRYSTAL_INP(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(CRYSTAL),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.CRYSTAL_isValid(INST)) return
      OBJ => ACRYSTALS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call CRYSTAL_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'CRYSTAL_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE CRYSTAL_INP


!---------------------------------------------------------
      SUBROUTINE CRYSTAL_OUT(INST,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      integer,intent(in) :: INST
      TYPE(CRYSTAL),POINTER :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      if (.not.CRYSTAL_isValid(INST)) return
      OBJ => ACRYSTALS(INST)%X
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call CRYSTAL_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'SAMPLE_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE CRYSTAL_OUT

!---------------------------------------------------------
      SUBROUTINE CRYSTAL_INP_R(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(CRYSTAL),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      character(LEN_NAME) :: CRNAME
      LR=1
      SELECT CASE (trim(PNAME))
! reflection name as integer !! = ord number on the list from crystal.lib
          case('TYPE'); OBJ%TYP=NINT(ARG(1))
          CASE('REFNAME')
        ! get enum type corresponding to this variable
            call GetEnumString('CRYSTAL_REFNAME',NINT(ARG(1)),CRNAME)
            call SET_CRYST_PARAM(OBJ,trim(CRNAME),.false.)

      !      call GETPARINDX(ENUM_TYPES,':','CRYSTAL_REFNAME',IE)
      !      call FINDSTRPAR(ENUM_TYPE_LIST(IE),':',NINT(ARG(1))+1,JS,JL)
      !      JL=MIN(JL,8)
      !      if (JL.gt.0) then
      !        call SET_CRYST_PARAM(OBJ,ENUM_TYPE_LIST(IE)(JS:JS+JL-1),.false.)
      !      else
      !        OBJ%REFNAME='undef'
      !      endif

! d-spacing
          CASE('DHKL');  OBJ%DHKL=ARG(1)
! sign of the take-off angle
! NOTE: ARG contains ordinal value, i.e. ARG=0 for SGN=-1
          CASE('SGN')
            if (ARG(1).le.0.D0) then
              OBJ%SGN=-1
            else
              OBJ%SGN=1
            endif
          !  write(*,*) 'CRYSTAL_INP SGN=',OBJ%SGN
! cutting angle [deg]
          CASE('CHI');   OBJ%CHI=ARG(1)*deg
! mosaicity [min]
          CASE('MOS')
!   write(*,*) 'CRYSTAL_INP ', ARG(1),minute,R8LN2
            OBJ%HMOS=ARG(1)*minute/R8LN2
            OBJ%VMOS=OBJ%ANIZ*OBJ%HMOS
! anisotropy of mosaicity (vertical/horizontal)
          CASE('ANIZ')
            OBJ%ANIZ=ARG(1)
            OBJ%VMOS=OBJ%ANIZ*OBJ%HMOS
! poisson number
          CASE('POISS');   OBJ%POI=ARG(1)
! curvatures
          case('RO') ; OBJ%RH=ARG(1)/1.D3;OBJ%RV=ARG(2)/1.D3;OBJ%RB=ARG(3)/1.D3; LR=3
! number of segments
          CASE('N')
            OBJ%NH=MAX(1,NINT(ARG(1)))
            OBJ%NV=MAX(1,NINT(ARG(2)))
            OBJ%NB=MAX(1,NINT(ARG(3)))
            LR=3
          CASE('NH');  OBJ%NH=MAX(1,NINT(ARG(1)))
          CASE('NV');  OBJ%NV=MAX(1,NINT(ARG(1)))
          CASE('NB');  OBJ%NB=MAX(1,NINT(ARG(1)))
! spacing between segments [mm]
          CASE('D');   OBJ%DH=ARG(1);OBJ%DV=ARG(1+1);OBJ%DB=ARG(1+2);LR=3
          CASE('DH');  OBJ%DH=ARG(1)
          CASE('DV');  OBJ%DV=ARG(1)
          CASE('DB');  OBJ%DB=ARG(1)
! d-gradient: magnitude (DGR) in [0.001/cm] and direction angle (DGA) in [deg]
          CASE('DGR');  OBJ%DGR=ARG(1)
          CASE('DGA');  OBJ%DGA=ARG(1)*deg
! extinction thickness [um]
          CASE('DLAM');  OBJ%DLAM=ARG(1)
! stacking angle [deg]
          CASE('ASTACK');  OBJ%ASTACK=ARG(1)*deg
          CASE('STACKH');   OBJ%STACKH=(ARG(1).eq.1.D0)
          CASE('STACKV');   OBJ%STACKV=(ARG(1).eq.1.D0)
! model type: (-1) ignored, (0) mosaic, (1) gradient, (2) mosaic + gradient
          CASE('MODEL');  OBJ%MODEL=NINT(ARG(1))
! settings for autofocusing
          case('FOCH'); OBJ%FOCH=(NINT(ARG(1)).eq.1)
          case('FOCV'); OBJ%FOCV=(NINT(ARG(1)).eq.1)
          CASE('FH1');  OBJ%FH1=ARG(1)
          CASE('FH2');  OBJ%FH2=ARG(1)
          CASE('FV1');  OBJ%FV1=ARG(1)
          CASE('FV2');  OBJ%FV2=ARG(1)
! automatically adjust in reflection position
          CASE('AUTOADJ');   OBJ%AUTOADJ=(ARG(1).eq.1.D0)
          CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
! adjust crystal type
  ! filter

      if (OBJ%NH.LE.0) then
        OBJ%TYP=ctyp_filter
  ! elastically bent
      else if (OBJ%HMOS.LT.sec) then
        OBJ%TYP=ctyp_bent
      else
  ! mosaic with gradient
        if (OBJ%DGR.NE.0.D0) then
          OBJ%TYP=ctyp_gradient
        else
  ! mosaic
          OBJ%TYP=ctyp_mosaic
        endif
      endif
      NARG=LR
      END SUBROUTINE CRYSTAL_INP_R

!---------------------------------------------------------
      SUBROUTINE CRYSTAL_OUT_R(OBJ,PNAME,ARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(CRYSTAL),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR,ITYPE,IORD
      LR=1
      SELECT CASE (trim(PNAME))
          CASE('TYPE');  ARG(1)=OBJ%TYP
! reflection name as integer !! = ord number on the list from crystal.lib
          CASE('REFNAME')
         ! get enum type corresponding to this variable
            call GetEnumOrd('CRYSTAL_REFNAME',trim(OBJ%REFNAME),ITYPE,IORD)
            ARG(1)=IORD
          !  JL=-1
          !  call GETPARINDX(ENUM_TYPES,':','CRYSTAL_REFNAME',JS)
          !  if (JS.gt.0) call GETPARINDX(ENUM_TYPE_LIST(JS),':',OBJ%REFNAME,JL)
           ! ARG(1)=JL-1
! d-spacing
          CASE('DHKL');  ARG(1)=OBJ%DHKL
! sign of the take-off angle
! ARG must return ORD number !!
          CASE('SGN')
            if (OBJ%SGN<0) then
              ARG(1)=0.D0
            else
              ARG(1)=1.D0
            endif
          !  write(*,*) 'CRYSTAL_OUT SGN=',OBJ%SGN
! cutting angle [deg]
          CASE('CHI');   ARG(1)=OBJ%CHI/deg
! mosaicity [min]
          CASE('MOS')
        ! write(*,*) 'CRYSTAL_OUT ', OBJ%HMOS,minute,R8LN2
            ARG(1)=OBJ%HMOS/minute*R8LN2
! anisotropy of mosaicity (vertical/horizontal)
          CASE('ANIZ')
            ARG(1)=OBJ%ANIZ
! poisson number
          CASE('POISS');   ARG(1)=OBJ%POI
! curvatures
          case('RO') ; ARG(1)=OBJ%RH*1.D3;ARG(2)=OBJ%RV*1.D3;ARG(3)=OBJ%RB*1.D3; LR=3
! number of segments
          CASE('N');   ARG(1)=OBJ%NH;ARG(1+1)=OBJ%NV;ARG(1+2)=OBJ%NB; LR=3
          CASE('NH');  ARG(1)=OBJ%NH
          CASE('NV');  ARG(1)=OBJ%NV
          CASE('NB');  ARG(1)=OBJ%NB
! spacing between segments [mm]
          CASE('D');   ARG(1)=OBJ%DH;ARG(1+1)=OBJ%DV;ARG(1+2)=OBJ%DB;LR=3
          CASE('DH');  ARG(1)=OBJ%DH
          CASE('DV');  ARG(1)=OBJ%DV
          CASE('DB');  ARG(1)=OBJ%DB
! d-gradient: magnitude (DGR) in [0.001/cm] and direction angle (DGA) in [deg]
          CASE('DGR');  ARG(1)=OBJ%DGR
          CASE('DGA');  ARG(1)=OBJ%DGA/deg
! extinction thickness [um]
          CASE('DLAM');  ARG(1)=OBJ%DLAM
! stacking angle [deg]
          CASE('ASTACK');  ARG(1)=OBJ%ASTACK/deg
          CASE('STACKH');   ARG(1)=LOG2INT(OBJ%STACKH)
          CASE('STACKV');   ARG(1)=LOG2INT(OBJ%STACKV)
! model type: (-1) ignored, (0) mosaic, (1) gradient, (2) mosaic + gradient
          CASE('MODEL');  ARG(1)=OBJ%MODEL
! settings for autofocusing
          case('FOCH');ARG(1)=0; if (OBJ%FOCH) ARG(1)=1
          case('FOCV');ARG(1)=0; if (OBJ%FOCV) ARG(1)=1
          CASE('FH1');  ARG(1)=OBJ%FH1
          CASE('FH2');  ARG(1)=OBJ%FH2
          CASE('FV1');  ARG(1)=OBJ%FV1
          CASE('FV2');  ARG(1)=OBJ%FV2
! automatically adjust in reflection position
          CASE('AUTOADJ');   ARG(1)=LOG2INT(OBJ%AUTOADJ)
          CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE CRYSTAL_OUT_R

      end MODULE CRYSTALS