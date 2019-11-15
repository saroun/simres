!//////////////////////////////////////////////////////////////////////
!////  $Id: sample_single.f,v 1.23 2019/08/15 15:02:07 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.23 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: SAMPLE - SINGLE CRYSTAL
!////
!////////////////////////////////////////////////////////////////////////
      module SAMPLE_SINGLE
      use CONSTANTS
      use SAMPLES
      use SAMPLES_TRACE
      implicit none

      type TSCRYST; SEQUENCE
  ! input data
        TYPE(TSAMPLE) :: SAM        ! ancestor class
        real(KIND(1.D0)) :: CELL_S(3)   ! unit cell size [A]
        real(KIND(1.D0)) :: CELL_A(3)   ! unit cell angles [rad]
        real(KIND(1.D0)) :: VECA(3)     ! vector along z in the scattering plane [rlu]
        real(KIND(1.D0)) :: VECB(3)     ! another vector in the scattering plane [rlu]
        real(KIND(1.D0)) :: ZERO(3)     ! zero positions for gonio (defines rotation of the AB system w.r.t. frame coord.)
        real(KIND(1.D0)) :: MOS         ! mosaicity [rad]
        real(KIND(1.D0)) :: TAU(3)      ! Bragg point / Brillouin zone center [rlu]
        real(KIND(1.D0)) :: QHKL(3)     ! phonon q (and the reference point of a dispersion plane) [rlu]
        real(KIND(1.D0)) :: EN          ! phonon energy at QHKL [meV]
        real(KIND(1.D0)) :: GHKL(3)     ! gradient direction of a disp. plane [rlu]
        real(KIND(1.D0)) :: GRAD        ! gradient magnitude of a disp. plane [meV/rlu]
        real(KIND(1.D0)) :: E(3)        ! eigenvector (assume real)
  ! auxilliary data
        real(KIND(1.D0)) :: U0(3,3)     ! U0 matrix, which relates UB with the frame basis coordinates
        real(KIND(1.D0)) :: UB(3,3)     ! UB matrix
        real(KIND(1.D0)) :: UBINV(3,3)  ! UB^-1
        real(KIND(1.D0)) :: UBFRAME(3,3)  ! conversion from hkl to frame basis coordinates = U0.UB
        real(KIND(1.D0)) :: SANG(3,3)   ! metric tensor for rec. lattice in [A]
        real(KIND(1.D0)) :: SRLU(3,3)   ! metric tensor for rec. lattice in [rlu]
        real(KIND(1.D0)) :: THETA       ! scattering angle
        real(KIND(1.D0)) :: PSI         ! angle between Ki and Q
        real(KIND(1.D0)) :: CELLVOL     ! unit cell volume
      end type TSCRYST

! global instance of TSCRYST
      TYPE(TSCRYST),target :: SCRYST

      contains

C------------------------------------
      SUBROUTINE SCRYST_DEFAULT(OBJ)
! default values for Si 133, EN=0
C------------------------------------
      TYPE (TSCRYST) :: OBJ
      integer :: K
        call SAMPLE_DEFAULT(OBJ%SAM)
        OBJ%SAM%FRAME%CLASS=SCLS_SINGLE
        OBJ%SAM%SIGA=4.75D-4
        OBJ%SAM%SIGI=1.85D-5
        OBJ%CELL_S=(/(5.4309D0,K=1,3)/)
        OBJ%CELL_A=(/(PI/2.D0,K=1,3)/)
        OBJ%VECA=(/2.D0,1.D0,1.D0/)
        OBJ%VECB=(/1.D0,1.D0,1.D0/)
        OBJ%ZERO=0.D0
        OBJ%MOS=0.D0
        OBJ%TAU=(/1.D0,3.D0,3.D0/)
        OBJ%QHKL=(/0.D0,0.D0,0.D0/)
        OBJ%EN=0.D0
        OBJ%GHKL=(/2.D0,1.D0,1.D0/)
        OBJ%GRAD=0.D0
        OBJ%E=(/1.D0,0.D0,0.D0/)
      END SUBROUTINE SCRYST_DEFAULT


!------------------------------------
      SUBROUTINE SCRYST_PREPARE(OBJ,IERR)
! prepare dependent fields after input
!------------------------------------
      type(TSCRYST) :: OBJ
      integer,intent(out) :: IERR
        IERR=0
        call FRAME_INIT_MAT(OBJ%SAM%FRAME)
        call SCRYST_UBMATRIX(OBJ,IERR)
      END SUBROUTINE SCRYST_PREPARE


C------------------------------------
      SUBROUTINE SCRYST_INIT(OBJ)
C------------------------------------
      TYPE (TSCRYST) :: OBJ
        call SAMPLE_INIT(OBJ%SAM)
        ! write(*,*) 'SCRYST_INIT '
      end SUBROUTINE SCRYST_INIT

!-----------------------------------------------------------------
      SUBROUTINE SCRYST_UBMATRIX(OBJ,IERR)
! Update UB matrix and dependences
!------------------------------------------- ----------------------
      TYPE (TSCRYST) :: OBJ
      integer,intent(out) :: ierr
      real(kind(1.D0)) :: Z,AUX(3),R1(3,3),R2(3,3)
      integer :: i
      ierr=1
    ! calculate UB matrix
      call CalUBmatrix(OBJ,ierr)
      if (ierr.ne.0) return
    ! calculate zero matrix
      call getRotYXY(OBJ%ZERO,R1)
      ! call PRN_MATRIX('ZERO',R1)

      call getRotYXY((/0.D0,0.5D0*PI,0.5D0*PI/),R2)
      call M3XM3(-1,R1,R2,OBJ%U0)
      ! call PRN_MATRIX('U0',OBJ%U0)

    ! get UBFRAME = conversion from hkl to local sample coord.
      call M3XM3(1,OBJ%U0,OBJ%UB,OBJ%UBFRAME)
    ! update Q0 and QSC
      call SCRYST_GETQSC(OBJ,OBJ%SAM%QSC)
      OBJ%SAM%Q0=SQRT(OBJ%SAM%QSC(1)**2+OBJ%SAM%QSC(2)**2+OBJ%SAM%QSC(3)**2)
    ! convert excitation surface gradient to [meV.A] in local coord.
      call M3xV3(1,OBJ%UBFRAME,OBJ%GHKL,AUX)
      Z=SQRT(AdotB(OBJ,OBJ%GHKL,OBJ%GHKL))/(AUX(1)**2+AUX(2)**2+AUX(3)**2)
      do i=1,3
        OBJ%SAM%EXC(i)=AUX(i)*OBJ%GRAD*Z
      ENDDO
      end SUBROUTINE SCRYST_UBMATRIX

!---------------------------------------------------------
      subroutine SCRYST_GETQE(OBJ,Q,E)
! return nominal Q,E values
!---------------------------------------------------------
      TYPE (TSCRYST) :: OBJ
      real(kind(1.D0)), intent(out) :: Q,E
      real(kind(1.D0)) :: QAB(3)
      call SCRYST_GETQAB(OBJ,QAB)
      Q=SQRT(QAB(1)**2+QAB(2)**2+QAB(3)**2)
      E=OBJ%EN
      end subroutine SCRYST_GETQE

!---------------------------------------------------------
      subroutine SCRYST_GETQSC(OBJ,QSC)
! return Q in local (FRAME) coordinates
! assumes valid UB matrix
!---------------------------------------------------------
      TYPE (TSCRYST) :: OBJ
      real(kind(1.D0)),intent(out) :: QSC(3)
      real(kind(1.D0)) :: AUX(3)
      AUX=OBJ%TAU+OBJ%QHKL
      call M3xV3(1,OBJ%UBFRAME,AUX,QSC)
      end subroutine SCRYST_GETQSC

!---------------------------------------------------------
      subroutine SCRYST_GETQAB(OBJ,QAB)
! return Q in AB coordinates
! assumes valid UB matrix
!---------------------------------------------------------
      TYPE (TSCRYST) :: OBJ
      real(kind(1.D0)),intent(out) :: QAB(3)
      real(kind(1.D0)) :: AUX(3)
1     format('SCRYST_GETQAB ',a,': ',6(1x,G12.5))
      AUX=OBJ%TAU+OBJ%QHKL
      call M3xV3(1,OBJ%UB,AUX,QAB)
      !write(*,1) 'QHKL',AUX
      !write(*,1) 'QAB',QAB
      end subroutine SCRYST_GETQAB

!-----------------------------------------------------------------
      SUBROUTINE SCRYST_DEPENDENCES(OBJ,IERR)
! Calxulate variables dependent on SAM settings
! Assumes valid UB matrix and SAM data
!------------------------------------------- ----------------------
      TYPE (TSCRYST) :: OBJ
      integer,intent(out) :: ierr
      real(kind(1.D0)) :: STHETA,CTHETA,SPSI,CPSI,KI0,KF0
      KI0=OBJ%SAM%KI0
      KF0=OBJ%SAM%KF0
      if ((KI0.gt.0).and.(KF0.gt.0)) then
        ierr=2
    ! check triangle
        CTHETA=(KI0**2+KF0**2-OBJ%SAM%Q0**2)/2.D0/KI0/KF0
        CPSI=(KI0**2+OBJ%SAM%Q0**2-KF0**2)/2.D0/KI0/OBJ%SAM%Q0
        if (abs(CTHETA).le.1.D0) then
          STHETA=OBJ%SAM%SGN*SQRT(1.D0-CTHETA**2)
          SPSI=OBJ%SAM%SGN*SQRT(1.D0-CPSI**2)
          OBJ%THETA=atan2(STHETA,CTHETA)
          OBJ%PSI=atan2(SPSI,CPSI)
          ierr=0
        endif
      endif
      select case(ierr)
      case(1)
        call MSG_ERROR('SCRYST_DEPENDENCES','Problem with UB matrix',0,1)
      case(2)
        call MSG_ERROR('SCRYST_DEPENDENCES','Can''t close scattering triangle',0,1)
      end select
      end SUBROUTINE SCRYST_DEPENDENCES


!---------------------------------------------------------
      SUBROUTINE SCRYST_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSCRYST) :: OBJ
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
          call SCRYST_INP_R(OBJ,trim(ARG%DEF%ID),NUM,LR)
        case default
          write(*,*) 'SCRYST_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SCRYST_INP

!---------------------------------------------------------
      SUBROUTINE SCRYST_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSCRYST) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call SCRYST_OUT_R(OBJ,trim(ARG%DEF%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'SCRYST_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SCRYST_OUT


C---------------------------------------------------------
      SUBROUTINE SCRYST_INP_R(OBJ,PNAME,ARG,NARG)
C input data from ARG to OBJ for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TSCRYST),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! unit cell size [A]
        CASE('CELLS');  OBJ%CELL_S(1:3)=ARG(1:3);LR=3
! unit cell angles deg -> rad
        CASE('CELLA');  OBJ%CELL_A(1:3)=ARG(1:3)*deg;LR=3
! vectors defining scattering plane
        CASE('VECA');   OBJ%VECA(1:3)=ARG(1:3);LR=3
        CASE('VECB');   OBJ%VECB(1:3)=ARG(1:3);LR=3
! zero positions for gonio
        CASE('ZERO');   OBJ%ZERO(1:3)=ARG(1:3)*deg;LR=3
! mosaicity min -> rad
        CASE('MOS');    OBJ%MOS=ARG(1)*minute/R8LN2
! Bragg point or Brillouin zone center [rlu]
        CASE('TAU');    OBJ%TAU(1:3)=ARG(1:3);LR=3
! Bragg point or the reference point of a dispersion plane [rlu]
        CASE('QHKL');   OBJ%QHKL(1:3)=ARG(1:3);LR=3
! energy reference point of a dispersion plane [meV]
        CASE('EN');     OBJ%EN=ARG(1)
! gradient direction of a disp. plane [rlu]
        CASE('GHKL');   OBJ%GHKL(1:3)=ARG(1:3);LR=3
! gradient magnitude of a disp. plane [meV/rlu]
        CASE('GRAD');   OBJ%EN=ARG(1)
! phonon eigenvector
        CASE('EIG');    OBJ%E(1:3)=ARG(1:3);LR=3
       CASE DEFAULT; CALL SAMPLE_INP_R(OBJ%SAM,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SCRYST_INP_R

C---------------------------------------------------------
      SUBROUTINE SCRYST_OUT_R(OBJ,PNAME,ARG,NARG)
C output data from OBJ to ARG for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TSCRYST),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! unit cell size [A]
        CASE('CELLS');  ARG(1:3)=OBJ%CELL_S(1:3);LR=3
! unit cell angles deg -> rad
        CASE('CELLA');  ARG(1:3)=OBJ%CELL_A(1:3)/deg;LR=3
! vectors defining scattering plane
        CASE('VECA');   ARG(1:3)=OBJ%VECA(1:3);LR=3
        CASE('VECB');   ARG(1:3)=OBJ%VECB(1:3);LR=3
! zero positions for gonio
        CASE('ZERO');   ARG(1:3)=OBJ%ZERO(1:3)/deg;LR=3
! mosaicity min -> rad
        CASE('MOS');    ARG(1)=OBJ%MOS/minute*R8LN2
! Bragg point or Brillouin zone center [rlu]
        CASE('TAU');    ARG(1:3)=OBJ%TAU(1:3);LR=3
! Bragg point or the reference point of a dispersion plane [rlu]
        CASE('QHKL');   ARG(1:3)=OBJ%QHKL(1:3);LR=3
! energy reference point of a dispersion plane [meV]
        CASE('EN');     ARG(1)=OBJ%EN
! gradient direction of a disp. plane [rlu]
        CASE('GHKL');   ARG(1:3)=OBJ%GHKL(1:3);LR=3
! gradient magnitude of a disp. plane [meV/rlu]
        CASE('GRAD');   ARG(1)=OBJ%EN
! eigenvector
        CASE('EIG'); ARG(1:3)=OBJ%E(1:3);LR=3
        CASE DEFAULT; CALL SAMPLE_OUT_R(OBJ%SAM,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE SCRYST_OUT_R

C--------------------------------------------------------------------
      SUBROUTINE CalUBmatrix(OBJ,IERR)
C  Calculate UB matrix and its inverse for given TSCRYST object
C--------------------------------------------------------------------
      TYPE (TSCRYST) :: OBJ
      integer,intent(out) :: ierr
      real(KIND(1.D0)), PARAMETER :: EPS=1.0D-10
      real(KIND(1.D0)) :: XCC,C1,C2,C3,A(3),B(3),XBB(3,3),V1(3),V2(3),V3(3)
      REAL*8 U(3,3),COSA(3),SINA(3),COSB(3),SINB(3)
      INTEGER :: I
      integer,parameter :: e(2,3)=reshape((/2,3,3,1,1,2/),(/2,3/))

      IERR=1
      XCC=0.D0
      DO I=1,3
        A(I)=OBJ%CELL_S(I)
        IF(ABS(A(I)).LT.EPS) goto 99   ! check lattice spacing
        COSA(I)=COS(OBJ%CELL_A(I))
        SINA(I)=SIN(OBJ%CELL_A(I))
        XCC=XCC+COSA(I)**2
      ENDDO
      XCC=1.+2.*COSA(1)*COSA(2)*COSA(3)-XCC
      IERR=2
      IF(XCC.LE.0.) goto 99    ! check lattice angles
      XCC=SQRT(XCC)
      OBJ%CELLVOL=XCC*A(1)*A(2)*A(3)   ! this is unit cell volume
      DO I=1,3
        B(I)=2.D0*PI*SINA(I)/(A(I)*XCC)  ! length of the r.l. axes in Ang-1
        COSB(I)=(COSA(e(1,I))*COSA(e(2,I))-COSA(I))/(SINA(e(1,I))*SINA(e(2,I)))
        SINB(I)=SQRT(1.-COSB(I)**2)  ! angles btw. r.l. axes
      ENDDO

C XBB(i,j) are projections of rec. lat. base vectors on the orthonormal base:
C Let (a,b,c) and (a*,b*,c*) are direct and reciprocal lattice base vectors
C assume a parallel to a*
C XBB(i,1) ... parallel to a*
C XBB(i,2) ... parallel to c x a*
C XBB(i,3) ... parallel to c
C i.e. the columns are a*,b*,c* in the new orthonormal base
      XBB(1,1)=B(1)
      XBB(2,1)=0.D0
      XBB(3,1)=0.D0
      XBB(1,2)=B(2)*COSB(3)
      XBB(2,2)=B(2)*SINB(3)
      XBB(3,2)=0.D0
      XBB(1,3)=B(3)*COSB(2)
      XBB(2,3)=-B(3)*SINB(2)*COSA(1)
      XBB(3,3)=2.D0*PI/A(3)

! convert A1,A2 to the new orthonormal system:
      call M3xV3(1,XBB,OBJ%VECA,V1)
      call M3xV3(1,XBB,OBJ%VECB,V2)
! get V3 perpendicular to V1,V2  : V3 = V1 x V2
      call V3cV3(V1,V2,V3)
! get V2 perpendicular to V3,V1
      call V3cV3(V3,V1,V2)
! get norms of V1,V2,V3
      C1=V1(1)**2+V1(2)**2+V1(3)**2
      C2=V2(1)**2+V2(2)**2+V2(3)**2
      C3=V3(1)**2+V3(2)**2+V3(3)**2
      C1=SQRT(C1)
      C2=SQRT(C2)
      C3=SQRT(C3)
      IERR=3
      IF (C1*C2*C3.LT.EPS) goto 99 ! check scattering plane
! U(i,j) is the orthonormal system made from V1,V2,V3 on rows (also called AB system in RESTRAX)
      DO I=1,3
          U(1,I)=V1(I)/C1
          U(2,I)=V2(I)/C2
          U(3,I)=V3(I)/C3
      ENDDO
! UB matrix: conversion from (hkl) to AB system in [1/A]
      call M3xM3(1,U,XBB,OBJ%UB)
      call roundMatrix(OBJ%UB)
      ! call PRN_MATRIX('UB',OBJ%UB)
! inverse UB
      I=INVERT(3,OBJ%UB,3,OBJ%UBINV,3)

!  SANG(i,j) defines scalar product in the orthonormal basis for (hkl) vectors
!  X.dot.Y= X^T.SANG.Y in [A^-2]
      DO I=1,3
        OBJ%SANG(I,I)=B(I)**2
      END DO
      OBJ%SANG(1,2)=B(1)*B(2)*COSB(3)
      OBJ%SANG(1,3)=B(1)*B(3)*COSB(2)
      OBJ%SANG(2,3)=B(2)*B(3)*COSB(1)
      OBJ%SANG(2,1)=OBJ%SANG(1,2)
      OBJ%SANG(3,1)=OBJ%SANG(1,3)
      OBJ%SANG(3,2)=OBJ%SANG(2,3)

!  SRLU(i,j) defines scalar product in the orthonormal basis for (hkl) vectors
!  X.dot.Y= X^T.SRLU.Y in [rlu]
      DO I=1,3
        OBJ%SRLU(I,I)=1.D0
      END DO
      OBJ%SRLU(1,2)=COSB(3)
      OBJ%SRLU(1,3)=COSB(2)
      OBJ%SRLU(2,3)=COSB(1)
      OBJ%SRLU(2,1)=OBJ%SRLU(1,2)
      OBJ%SRLU(3,1)=OBJ%SRLU(1,3)
      OBJ%SRLU(3,2)=OBJ%SRLU(2,3)

      IERR=0
      return
c// handle errors
99    select case (IERR)
      case(1)
        call MSG_ERROR('UB matrix','Check Lattice Spacings',0,1)
      case(2)
        call MSG_ERROR('UB matrix','Check Cell Angles',0,1)
      case(3)
        call MSG_ERROR('UB matrix','Check Scattering Plane',0,1)
      end select

      END SUBROUTINE CalUBmatrix


!----------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION AdotB(OBJ,A,B)
! dot product of hkl vectors (result in [A^-2])
!----------------------------------------------------
      TYPE (TSCRYST) :: OBJ
      REAL(KIND(1.D0)),intent(in) :: A(3),B(3)
      REAL(KIND(1.D0)) :: C(3)
      call M3xV3(-1,OBJ%SANG,A,C)
      AdotB=C(1)*B(1)+C(2)*B(2)+C(3)*B(3)
      end FUNCTION AdotB

      end module SAMPLE_SINGLE