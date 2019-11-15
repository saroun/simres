!//////////////////////////////////////////////////////////////////////
!////  $Id: samples.f,v 1.47 2019/06/22 20:48:05 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.47 $
!////     $Date: 2019/06/22 20:48:05 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: SAMPLE
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SAMPLES
      use CONSTANTS
      use FIELDDATA
      use CLASSES
      use FRAMES
      implicit none

! SAMPLE
! values of TYP:
      integer,parameter :: isam_INELASTIC=0
      integer,parameter :: isam_ELASTIC=1
      integer,parameter :: isam_POWDER=2
      integer,parameter :: isam_VANAD=3
      integer,parameter :: isam_PHONON=4

! values of COORD:
      integer,parameter :: csam_PWD=0
      integer,parameter :: csam_SANS=1

      TYPE TSAMPLE; SEQUENCE
        TYPE(TFRAME) :: FRAME
        real(KIND(1.D0)) :: THICK     ! wall thickness for hollow shapes [mm]
        real(KIND(1.D0)) :: SIGSC ! total scattering  cross-section [1/mm] (! TODO: modify for phonons !!)
        real(KIND(1.D0)) :: SIGA      ! total absorption cross-section [1/mm/A]
        real(KIND(1.D0)) :: SIGI      ! total incoherent cross-section [1/mm]
        real(KIND(1.D0)) :: TEMP      ! temperature [K]
        real(KIND(1.D0)) :: MAG(3)    ! magnetic field [gauss]
        real(KIND(1.D0)) :: STRAIN(3) ! true strain []
        real(KIND(1.D0)) :: PRESS     ! isotropic preasure [MPa]
        INTEGER :: TYP    ! sample type
        integer :: SGN    ! sign of the scattering angle
        logical :: TRANS  ! true if transmission is allowed
        logical :: SCATT  ! true if scattering is allowed
! auxilliary parameters, not accessible by user input
        real(KIND(1.D0)) :: QSC(3)  ! nominal scattering vector in local coordinates
        real(KIND(1.D0)) :: ESC     ! nominal energy transfer [meV]
        real(KIND(1.D0)) :: EXC(3)  ! planar excitation in local coordinates (ref. to QSC,ESC) = h_bar*group_velocity [meV*Ang]
        real(KIND(1.D0)) :: KI0,KF0 ! nominal magnitudes of KI,KF
        real(KIND(1.D0)) :: Q0      ! =|QSC|
        real(KIND(1.D0)) :: QMAT(3,3) ! conversion matrix local -> Q-coordinates
        real(KIND(1.D0)) :: MCL(3,3) ! conversion matrix entry (=Lab) -> Cooper & Nathans
        real(KIND(1.D0)) :: QN(3) ! nominal Q vector in Q-coordinates
        real(KIND(1.D0)) :: KF(3) ! nominal KF vector in Q-coordinates
        integer :: COORD ! choice of polar axis for scattering angles, (0) z-axis, like in powder dif., (1) y-axis, like in SANS
        INTEGER :: IRND(3)   ! indices to random numbers array (theta,phi,dE)
      END TYPE TSAMPLE

! overloaded SAMPLE_ADJAX
      interface SAMPLE_ADJAX
        module procedure SAMPLE_ADJAX_1
        module procedure SAMPLE_ADJAX_2
      end interface

      private SAMPLE_ADJAX_1
      private SAMPLE_ADJAX_2

! global instance of SAMPLE
      type(TSAMPLE),TARGET :: SAMPLE

      logicAL,PRIVATE :: dbg=.false.

      contains

C------------------------------------
      SUBROUTINE SAMPLE_DEFAULT(OBJ)
C------------------------------------
        type(TSAMPLE) :: OBJ
        call FRAME_CLEAR(OBJ%FRAME)
        OBJ%FRAME%CLASS=SCLS_SAMPLE
        OBJ%FRAME%NAME='default'
        OBJ%FRAME%SHAPE=FRAME_SHAPE_CYLINDER
        OBJ%FRAME%SIZE=(/10.D0,20.D0,10.D0/)
        OBJ%SGN=1
        OBJ%THICK=0.D0
        OBJ%SIGSC=1.D-2
        OBJ%SIGA=0.D0
        OBJ%SIGI=0.D0
        OBJ%TEMP=300.D0
        OBJ%MAG=0.D0
        OBJ%STRAIN=0.D0
        OBJ%PRESS=0.D0
        OBJ%QSC=(/0.D0,0.D0,1.D0/)
        OBJ%ESC=0.D0
        OBJ%EXC=0.D0
        OBJ%IRND=0
        OBJ%TRANS=.false.
        OBJ%SCATT=.true.
        call UNIMAT(OBJ%QMAT,3,3)
      END SUBROUTINE SAMPLE_DEFAULT


!---------------------------------------------------------------------
      SUBROUTINE SAMPLE_INIT(OBJ)
!----------------------------------------------------------------------
      type(TSAMPLE) :: OBJ
      CALL FRAME_INIT(OBJ%FRAME)
      OBJ%IRND=0
      END SUBROUTINE SAMPLE_INIT

!------------------------------------
      SUBROUTINE SAMPLE_PREPARE(OBJ,IERR)
! prepare dependent fields after input
!------------------------------------
      type(TSAMPLE) :: OBJ
      integer,intent(out) :: IERR
        IERR=0
        call FRAME_INIT_MAT(OBJ%FRAME)
      END SUBROUTINE SAMPLE_PREPARE

!-----------------------------------------------------------------
      SUBROUTINE SAMPLE_ADJAX_1(OBJ,KI0,KF0,Q0,PHI,KFMODE,IERR)
! Adjust sample for given KI0,KF0,Q0,PHI,KFMODE
! 1) Rotates AX to match the scattering angle
! 2) Does not move gonio, only updates transformation matrices
! Following fields of OBJ are updated here:
! FRAME%AX,KI0,KF0,KI0,ESC + dependences in SETQSC
!------------------------------------------- ----------------------
      TYPE (TSAMPLE) :: OBJ
      real(kind(1.D0)),intent(in) :: KI0,KF0,Q0,PHI
      integer,intent(in) :: KFMODE
      integer,intent(out) :: ierr
      real(kind(1.D0)) :: AX(3), AXIS(3)
      call GetKfAngles(KI0,KF0,Q0,OBJ%SGN,PHI,KFMODE,AX,IERR)
      if (IERR==0) then
        select case(KFMODE)
        ! PHI is sagital angle (X-craddle)
        case(KFMODE_FLAT)
          AXIS=(/AX(1),AX(2),OBJ%FRAME%AX(3)/)
        ! PHI is azimuthal angle (rotation of Kf around Ki)
        case(KFMODE_TOF)
          AXIS=AX
        end select
        call SetAxisAngles(OBJ%FRAME,AXIS)
      endif
      end SUBROUTINE SAMPLE_ADJAX_1

!-----------------------------------------------------------------
      SUBROUTINE SAMPLE_ADJAX_2(OBJ,KI0,THETA,PHI,KFMODE)
! Adjust sample for given KI0,THETA,PHI,KFMODE
! version for elastic scattering from a polycrystal
! 1) Rotates AX to match the scattering angle
! 2) Does not move gonio, only updates transformation matrices and QSC
! Following fields of OBJ are updated here:
! FRAME%AX,KI0,KF0=KI0,ESC=0 + dependences in SETQSC
!------------------------------------------- ----------------------
      TYPE (TSAMPLE) :: OBJ
      real(kind(1.D0)),intent(in) :: KI0,THETA,PHI
      integer,intent(in) :: KFMODE
      real(kind(1.D0)) :: AX(3)
      call CalAxisAngles(THETA,PHI,KFMODE,AX)
      select case(KFMODE)
    ! PHI is sagital angle (X-craddle)
      case(KFMODE_FLAT)
        ! preserve AX(3) as the scanning angle
          call SetAxisAngles(OBJ%FRAME,(/AX(1),AX(2),OBJ%FRAME%AX(3)/))
    ! PHI is azimuthal angle (rotation of Kf around Ki)
      case(KFMODE_TOF)
          call SetAxisAngles(OBJ%FRAME,AX)
      end select
      end SUBROUTINE SAMPLE_ADJAX_2

!------------------------------------------------------------
      SUBROUTINE SAMPLE_ORIENT(OBJ,XStayHorizontal)
! Orient sample by setting goniometer so that QSC is aligned with Q=KF-KI
! Assume:
! 1) transformation betwen incident and exit coordinates (REXI) is already defined
! 2) KI0,KF0 and QSC are already defined
! 3) KI is parallel to (0,0,1) of the incident coordinates
! 4) KF is parallel to (0,0,1) of the exit coordinates
!------------------------------------------------------------
      type(TSAMPLE) :: OBJ
      logical :: XStayHorizontal
      real(kind(1.D0)) :: KI(3),KF(3),VQ(3),GON(3)
      integer :: i
      KI=(/0.D0,0.D0,1.D0/)
    ! get KF in incident coordinates
      call M3xV3(-OBJ%FRAME%MEXI,OBJ%FRAME%REXI,KI,KF)
      do i=1,3
        KI(i)=KI(i)*OBJ%KI0
        KF(i)=KF(i)*OBJ%KF0
      enddo
      VQ=KF-KI
      !call PRN_VECTOR('SAMPLE_ORIENT, VQ',VQ)
      !call PRN_VECTOR('SAMPLE_ORIENT, QSC',OBJ%QSC)

      call CalGonioAngles(VQ,OBJ%QSC,XStayHorizontal,GON)
      call SetGonioAngles(OBJ%FRAME,GON)
      end SUBROUTINE SAMPLE_ORIENT

!------------------------------------------------------------
      SUBROUTINE SAMPLE_SETQSC(OBJ)
! Set scattering vector QSC [local coordinates], does not move gonio
! Assume:
! 1) transformation betwen incident and exit coordinates (REXI) is already defined
!    done by calling SetAxisAngles
! 2) transformation betwen incident and local coordinates (RLOC) is already defined
!    done by calling SetGonioAngles
! 3) KI0,KF0 are already defined
! 4) KI is parallel to (0,0,1) of the incident coordinates
! 5) KF is parallel to (0,0,1) of the exit coordinates
! Following fields of OBJ are updated here:
! MCL,QSC,QMAT,KF,QN
!------------------------------------------------------------
      type(TSAMPLE) :: OBJ
      real(kind(1.D0)) :: KI(3),KF(3),VQ(3),U(3,3),R(3,3),R1(3,3),VV(3)
      integer :: i
  ! get KI,KF,Q in Lab (entry) coordinates
      KI=(/0.D0,0.D0,1.D0/)
      call M3xV3(-OBJ%FRAME%MEXI,OBJ%FRAME%REXI,KI,KF)
      do i=1,3
        KI(i)=KI(i)*OBJ%KI0
        KF(i)=KF(i)*OBJ%KF0
      enddo
      VQ=KF-KI
      OBJ%Q0=SQRT(VQ(1)**2+VQ(2)**2+VQ(3)**2)
  ! get unit vector along KI x VQ in Lab coord.
      call V3cV3(KI,VQ,VV)
  ! get trans. matrix from Lab to C&N coordinates
  ! NOTE: we adopt the convention Qx=KF-KI, Qz vertical for C&N coordinates
      call getRotZX(OBJ%SGN*VV,VQ,OBJ%MCL)
  ! convert Q to local coordinates
      call  M3xV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,VQ,OBJ%QSC)

  ! get QMAT: conversion from local to Q coordinates
  ! Q coordinates are defined by z'//KI, y'//KIxQSC
      call getRotZX(KI,VQ,R)  ! R is QMAT for entry coordinates
      call UNIMAT(U,3,3)
      CALL M3XM3(-1,OBJ%FRAME%RLOC,U,R1)  ! get RLOC^T
      CALL M3XM3(1,R,R1,OBJ%QMAT)         ! QMAT=R.RLOC^T

    ! get KF in Q-coordinates
      call  M3xV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,KF,VQ)
      CALL M3XV3(1,OBJ%QMAT,VQ,OBJ%KF)
    ! get QSC in Q-coordinates
      CALL M3XV3(1,OBJ%QMAT,OBJ%QSC,OBJ%QN)

      end SUBROUTINE SAMPLE_SETQSC


!---------------------------------------------------------
      subroutine SAMPLE_GETQE(OBJ,Q,E)
! return nominal Q,E values
!---------------------------------------------------------
      TYPE (TSAMPLE) :: OBJ
      real(kind(1.D0)), intent(out) :: Q,E
      Q=OBJ%Q0
      E=OBJ%ESC
      end subroutine SAMPLE_GETQE


!-------------------------------------------
      subroutine SAMPLE_SETHOLLOW(OBJ, t)
! Sets wall thickness for hollow shape
!-------------------------------------------
      type(TFRAME) :: OBJ
      real(kind(1.D0)),intent(in) :: t
      integer :: i
      if (t.le.0.D0) then
        OBJ%INSIZE = 0.D0
        OBJ%ISHOLLOW = .false.
        return
      endif
      select case (OBJ%SHAPE)
      case(FRAME_SHAPE_BOX)
        OBJ%INSIZE = OBJ%SIZE - 2*(/t, 0.D0, t/)
        OBJ%ISHOLLOW = .true.
      case(FRAME_SHAPE_DISC)
        OBJ%INSIZE = OBJ%SIZE - 2*(/t, t, 0.D0/)
        OBJ%ISHOLLOW = .true.
      case(FRAME_SHAPE_CYLINDER)
        OBJ%INSIZE = OBJ%SIZE - 2*(/t, 0.D0, t/)
        OBJ%ISHOLLOW = .true.
      case(FRAME_SHAPE_ELLIPSOID)
        OBJ%INSIZE = OBJ%SIZE - 2*(/t, t, t/)
        OBJ%ISHOLLOW = .true.
      case DEFAULT
        OBJ%INSIZE = 0.D0
        OBJ%ISHOLLOW = .false.
      END select
      if (OBJ%ISHOLLOW) then
        do i=1,3
          OBJ%INSIZE(i) = max(0.D0,OBJ%INSIZE(i))
        enddo
      endif
      END subroutine SAMPLE_SETHOLLOW



!---------------------------------------------------------
      SUBROUTINE SAMPLE_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSAMPLE) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call SAMPLE_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'SAMPLE_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SAMPLE_INP

!---------------------------------------------------------
      SUBROUTINE SAMPLE_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSAMPLE) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call SAMPLE_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'SAMPLE_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SAMPLE_OUT


C---------------------------------------------------------
      SUBROUTINE SAMPLE_INP_R(OBJ,PNAME,ARG,NARG)
C input data from ARG to OBJ for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      type(TSAMPLE) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
  !    write(*,*) 'SAMPLE_INP_R ',trim(PNAME),' ARG=',ARG(1)
      SELECT CASE (trim(PNAME))
! sample type
        CASE('TYPE');  OBJ%TYP=NINT(ARG(1)); !write(*,*) 'SAMPLE_INP TYPE=',OBJ%TYP
! sign of the scattering angle
        CASE('SGN')
            if (ARG(1).le.0.D0) then
              OBJ%SGN=-1
            else
              OBJ%SGN=1
            endif
! total coherent scattering cross-section [1/cm]
        CASE('SIGSC');  OBJ%SIGSC=abs(ARG(1))/10
! total absorption cross-section [1/cm/A]
        CASE('SIGA');   OBJ%SIGA=abs(ARG(1))/10
! total incoherent cross-section [1/cm]
        CASE('SIGI');   OBJ%SIGI=abs(ARG(1))/10
! temperature [K]
        CASE('TEMP');   OBJ%TEMP=abs(ARG(1))
! magnetic field [gauss]
        CASE('MAG');   OBJ%MAG(1:3)=ARG(1:3);LR=3
! strain [], direction in local coord.
        CASE('STRAIN');   OBJ%STRAIN(1:3)=ARG(1:3);LR=3
! isotropic preasure [MPa]
        CASE('PRESS');   OBJ%PRESS=ARG(1)
! transmission allowed
        CASE('TRANS');   OBJ%TRANS=(NINT(ARG(1))==1)
! scattering allowed
        CASE('SCATT');   OBJ%SCATT=(NINT(ARG(1))==1)
! thickness
        CASE('THICK')
          call SAMPLE_SETHOLLOW(OBJ%FRAME, ARG(1))
          if (OBJ%FRAME%ISHOLLOW) then
            OBJ%THICK=ARG(1)
          else
            OBJ%THICK=0.D0
          endif

! try FRAME parameters
        CASE DEFAULT; CALL FRAME_INP_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SAMPLE_INP_R

C---------------------------------------------------------
      SUBROUTINE SAMPLE_OUT_R(OBJ,PNAME,ARG,NARG)
C output data from OBJ to ARG for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      type(TSAMPLE) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! sample type
        CASE('TYPE');  ARG(1)=OBJ%TYP
! sign of the scattering angle
        CASE('SGN')
            if (OBJ%SGN<0) then
              ARG(1)=0.D0
            else
              ARG(1)=1.D0
            endif
! total coherent scattering cross-section [1/cm]
        CASE('SIGSC');  ARG(1)=abs(OBJ%SIGSC)*10
! total absorption cross-section [1/cm]
        CASE('SIGA');   ARG(1)=abs(OBJ%SIGA)*10
! total incoherent cross-section [1/cm]
        CASE('SIGI');   ARG(1)=abs(OBJ%SIGI)*10
! temperature [K]
        CASE('TEMP');   ARG(1)=abs(OBJ%TEMP)
! magnetic field [gauss]
        CASE('MAG');   ARG(1:3)=OBJ%MAG(1:3);LR=3
! strain [], direction in local coord.
        CASE('STRAIN');   ARG(1:3)=OBJ%STRAIN(1:3);LR=3
! isotropic preasure [MPa]
        CASE('PRESS');   ARG(1)=OBJ%PRESS
! transmission allowed
        CASE('TRANS');   ARG(1)=LOG2INT(OBJ%TRANS)
! scattering allowed
        CASE('SCATT');   ARG(1)=LOG2INT(OBJ%SCATT)
! thickness
        CASE('THICK');   ARG(1)=OBJ%THICK
! try FRAME parameters
        CASE DEFAULT; CALL FRAME_OUT_R(OBJ%FRAME,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SAMPLE_OUT_R

      end MODULE SAMPLES
