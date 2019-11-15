!//////////////////////////////////////////////////////////////////////
!////  $Id: spectrometer.f,v 1.42 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.42 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Basic instrument interface - TSPEC
!////  Provides high-level control of beamline components
!////  Other specialized interfaces (e.g. TAS, PWD, ...) are descendants of TSPEC
!////
!//////////////////////////////////////////////////////////////////////
      module SPECTROMETER
      use CONSTANTS
      use CLASSES
      use FIELDDEF
      use FIELDDATA
      use XMLINFO
      USE DATASETS
      USE SIMATH
      use FRAMES
      USE TRACINGDATA
      use FRAMES_TRACE
      implicit none

! TSPEC envelops data needed for basic control of a spectrometer
! NOTE: scattering triangle is uniquely defined by KI,KF,Q,PSI,SS
      type TSPEC; sequence
        TYPE(TFRAME) :: FRAME
        !character(LEN_ID)   :: ID   ! ID string
        !character(LEN_NAME) :: NAME ! descriptive name
        integer :: CLASS            ! class ID number
        real(KIND(1.D0)) :: KI,KF,Q ! scattering triangle
        real(KIND(1.D0)) :: EN,LAMBDA, TTHETA  ! 2nd version of scattering triangle definition
        real(KIND(1.D0)) :: PSI     ! rotation of KF axis (definition depends on KFMODE)
        integer :: SS               ! sign of sample take-off angle
        integer :: FIX              ! fixed KI (1) or KF (2), when changing energy transfer
        integer :: KFMODE           ! mode of the KF axis (flat or TOF), defines the meaning of PSI
        INTEGER :: FLATCONE         ! normal(0) or flat-cone (1) analyzer mode
        LOGICAL :: OffPlane         ! if true, off-plane Q is adjusted by lifting Kf from horizontal plane
                                    ! otherwise adjusted by sample rotation only
        LOGICAL :: ADJSPEC          ! automaticaly adjust for sample reflection
        LOGICAL :: ORISAM           ! automaticaly orient sample when relevant (single-crystals)
        LOGICAL :: QCONST           ! keep Q,E const. when adjusting ThetaS
        integer :: INPSET           ! input set: (0) ki or kf,fix,Q,EN (1) lambda, theta, EN
        TYPE(TDATASET) :: NOM       ! defines nominal setting
      end type TSPEC

! global instance of TSPEC
      TYPE(TSPEC),target :: SPEC

      contains

!--------------------------------------------------------
! Initialization
!--------------------------------------------------------
      subroutine SPEC_INIT(SPEC)
      TYPE(TSPEC) :: SPEC
      call FRAME_INIT(SPEC%FRAME)
      end subroutine SPEC_INIT

!--------------------------------------------------------
! Set default values
!--------------------------------------------------------
      subroutine SPEC_DEFAULT(SPEC)
      TYPE(TSPEC) :: SPEC
        SPEC%CLASS=DCLS_SPEC
        SPEC%KI=2.5D0
        SPEC%KF=2.5D0
        SPEC%Q=1.D0
        SPEC%PSI=0.D0
        SPEC%SS=-1
        SPEC%FIX=1
        SPEC%LAMBDA=2.D0
        SPEC%TTHETA=PI/2
        SPEC%EN=0.D0
        SPEC%FLATCONE=0
        SPEC%INPSET=0
        SPEC%OffPlane=.false.
        SPEC%KFMODE=0
        SPEC%ADJSPEC=.true.
        SPEC%ORISAM=.true.
        SPEC%QCONST=.false.
        call FRAME_CLEAR(SPEC%FRAME)
        ! make sure distance is zero, but this should be done in FRAME_DEFAULT
        SPEC%FRAME%DIST=0.D0
      end subroutine SPEC_DEFAULT

!-----------------------------------------------------------------
      SUBROUTINE SPEC_ADJAX(SPEC,KI0,KF0,Q0,PHI,KFMODE,IERR)
! Adjust angles for given KI0,KF0,Q0,PHI,KFMODE
! 1) Rotates AX to match the scattering angle
! 2) Does not move gonio, only updates transformation matrices
!------------------------------------------- ----------------------
      TYPE (TSPEC) :: SPEC
      real(kind(1.D0)),intent(in) :: KI0,KF0,Q0,PHI
      integer,intent(in) :: KFMODE
      integer,intent(out) :: ierr
      real(kind(1.D0)) :: AX(3), AXIS(3)
      call GetKfAngles(KI0,KF0,Q0,SPEC%SS,PHI,KFMODE,AX,IERR)
      if (IERR==0) then
        select case(KFMODE)
        ! PHI is sagital angle (X-craddle)
        case(KFMODE_FLAT)
          AXIS=(/AX(1),AX(2),SPEC%FRAME%AX(3)/)
        ! PHI is azimuthal angle (rotation of Kf around Ki)
        case(KFMODE_TOF)
          AXIS=AX
        end select
        call SetAxisAngles(SPEC%FRAME,AXIS)
      endif
      end SUBROUTINE SPEC_ADJAX

!-------------------------------------------------
      SUBROUTINE SPEC_ADJUST(SPEC,IERR)
! Update position of the component using nominal trajectory represented by NEUT
! Call FRAME_INIT_MAT first.
! Then use SLIT_PRE & SLIT_POST procedures to propagate NEUT
!-------------------------------------------------
      TYPE(TSPEC) :: SPEC
      integer,intent(out) :: IERR
      REAL(kind(1.D0)) :: Z
      integer :: i
1     format(a,': ',7(G13.6,1x))
      IERR=0
      call FRAME_INIT_MAT(SPEC%FRAME)
      !write(*,1) 'SPEC_ADJUST to a:  ',NEUT%K0*NEUT%T, NEUT%R
      call SLIT_PRE(SPEC%FRAME)
      !write(*,1) 'SPEC_ADJUST to b:  ',NEUT%K0*NEUT%T, NEUT%R
      call SLIT_TURN(SPEC%FRAME)
      !write(*,1) 'SPEC_ADJUST to c:  ',NEUT%K0*NEUT%T, NEUT%R
      ! change energy
      if (NEUT%K0.gt.0.D0) then
        Z=SPEC%KF/NEUT%K0
        do i=1,3
          NEUT%K(i)=NEUT%K(i)*Z
        enddo
        NEUT%K0=SPEC%KF
      endif
      call SLIT_POST(SPEC%FRAME)
      end SUBROUTINE SPEC_ADJUST

!---------------------------------------------------------
      SUBROUTINE SPEC_SET_THETA(SPEC,TTHETA)
! setter for scattering angle, adjusts Q
! requires KI,KF
!---------------------------------------------------------
      TYPE(TSPEC),intent(out) :: SPEC
      REAL(KIND(1.D0)),intent(in) :: TTHETA
      if (TTHETA.EQ.0.D0) then
        SPEC%Q=abs(SPEC%KF-SPEC%KI)
      else if ((SPEC%KI.gt.0.D0).and.(SPEC%KF.gt.0.D0)) then
        SPEC%Q=SQRT(SPEC%KI**2+SPEC%KF**2-2.D0*SPEC%KI*SPEC%KF*cos(TTHETA))
      endif
      SPEC%TTHETA=TTHETA
      end SUBROUTINE SPEC_SET_THETA

!---------------------------------------------------------
      SUBROUTINE SPEC_SET_THETA_QCONST(SPEC,THETA,IERR)
! setter for scattering angle, keeping Q,E constant
! change wavelength instead
!---------------------------------------------------------
      TYPE(TSPEC),intent(inout) :: SPEC
      REAL(KIND(1.D0)),intent(in) :: THETA
      integer ,intent(out) :: IERR
      REAL(KIND(1.D0)) :: CC,SS,E,Q2,KI1,KI2,KI,Z,DT,Z1,Z2
      ierr=1
      CC=cos(theta)
      SS=sin(theta)
      E=SPEC%KI**2-SPEC%KF**2
      Q2=SPEC%Q**2
      if (Q2<1.D-12) then
        SPEC%KF=SPEC%KI
        SPEC%Q=0.D0
        return
      endif
      if (SS**2<1.D-12) then
        KI=(Q2+E)**2/Q2/4.D0
      else
        DT=Q2**2-(E*SS)**2
        Z=Q2+E*SS**2
        if (DT.GT.0.D0) then
          KI1=(Z+CC*SQRT(DT))/(2.D0*SS**2)
          KI2=(Z-CC*SQRT(DT))/(2.D0*SS**2)
          Z1=abs(2.D0*KI1-E-2*SQRT(abs(KI1)*abs(KI1-E))*CC-Q2)
          Z2=abs(2.D0*KI2-E-2*SQRT(abs(KI2)*abs(KI2-E))*CC-Q2)
          if (Z1<Z2) then
            KI=KI1
          else
            KI=KI2
          endif
!1       format(a,6(2x,G12.5))
          !write(*,*) '  Z, DT,KI=',Z,DT,KI,twopi/sqrt(abs(KI))
          !write(*,*) '  CC,SS=',CC,SS
          !write(*,*) '  KI1,KI2=',KI1,KI2
          !write(*,*) '  Z1,Z2,Q2=',Z1,Z2,Q2
        else if (DT==0.D0) then
          KI=Z/(2.D0*SS)
        else
          KI=-1.D0
        endif
        if (KI>max(0.D0,E)) then
          SPEC%KI=SQRT(KI)
          SPEC%KF=SQRT(KI-E)
          ierr=0
        endif
      endif
      end SUBROUTINE SPEC_SET_THETA_QCONST

!---------------------------------------------------------
      SUBROUTINE SPEC_SET_LAMBDA(SPEC,LAMBDA,EN)
! adjust KI,KF for wavelength and energy transfer
! sets also FIX=0, so that KI=2pi/lambda
!---------------------------------------------------------
      TYPE(TSPEC),intent(out) :: SPEC
      REAL(KIND(1.D0)),intent(in) :: LAMBDA,EN
      integer :: ierr
      ierr=1
      if (LAMBDA.gt.0) then
        SPEC%FIX=0
        SPEC%KI=2.D0*PI/LAMBDA
        call SPEC_SET_EN(SPEC,EN,ierr)
        !  in case of error, set EN=0, which always works
        if (ierr.ne.0) then
          call SPEC_SET_EN(SPEC,0.D0,ierr)
        endif
        ! if error, at least assign the lambda,EN values
        ! but KI,KF may not be consistet
        if (ierr.ne.0) then
          SPEC%LAMBDA=LAMBDA
          SPEC%EN=EN
        endif
      endif
      end SUBROUTINE SPEC_SET_LAMBDA

!---------------------------------------------------------
      SUBROUTINE SPEC_GET_THETA(SPEC,THETA,IERR)
! getter for scattering angle
!---------------------------------------------------------
      TYPE(TSPEC),intent(in) :: SPEC
      REAL(KIND(1.D0)),intent(out) :: THETA
      integer ,intent(out) :: IERR
      REAL(KIND(1.D0)) :: CTHETA,STHETA
      ierr=1
      if(SPEC%Q<=0.D0) then
        THETA=0.D0
        ierr=0
      else
        CTHETA=(SPEC%KI**2+SPEC%KF**2-SPEC%Q**2)/2.D0/SPEC%KI/SPEC%KF
        if (abs(CTHETA).lt.1.D0) then
          STHETA=SPEC%SS*SQRT(1.D0-CTHETA**2)
          THETA=atan2(STHETA,CTHETA)
          ierr=0
        endif
      endif
      end SUBROUTINE SPEC_GET_THETA


!---------------------------------------------------------
      SUBROUTINE SPEC_SETQ(SPEC,Q,IERR)
! setter for Q
!---------------------------------------------------------
      TYPE(TSPEC),intent(inout) :: SPEC
      REAL(KIND(1.D0)),intent(in) :: Q
      integer ,intent(out) :: IERR
      REAL(KIND(1.D0)) :: CTHETA,STHETA
1     format(a,': ',6(G12.5,1x))
      ierr=1
      if(Q<=0.D0) then
        SPEC%TTHETA=0.D0
        ierr=0
      else
        CTHETA=(SPEC%KI**2+SPEC%KF**2-Q**2)/2.D0/SPEC%KI/SPEC%KF
        if (abs(CTHETA).lt.1.D0) then
          STHETA=SPEC%SS*SQRT(1.D0-CTHETA**2)
          SPEC%TTHETA=atan2(STHETA,CTHETA)
          SPEC%Q=Q
          ierr=0
        endif
      endif
      end SUBROUTINE SPEC_SETQ

!---------------------------------------------------------
      SUBROUTINE SPEC_SET_EN(SPEC,EN,IERR)
! Setter for energy transfer. makes KI, KF,EN and LAMBDA
! consistent, or returns ierr!=0
!---------------------------------------------------------
      TYPE(TSPEC),intent(inout) :: SPEC
      REAL(KIND(1.D0)),intent(in) :: EN
      integer ,intent(out) :: IERR
      REAL(KIND(1.D0)) :: ZI,ZF
      ierr=1
      ZI=0.D0
      ZF=0.D0
      select case(SPEC%FIX)
      case(0)
        ZF=SPEC%KI**2-EN/HSQOV2M
        if (ZF.gt.0) then
          ZF=SQRT(ZF)
          ZI=SPEC%KI
          if (ZI.gt.0) ierr=0
        endif
      case(1)
        ZI=SPEC%KF**2+EN/HSQOV2M
        if (ZI.gt.0) then
          ZI=SQRT(ZI)
          ZF=SPEC%KF
          if (ZF.gt.0) ierr=0
        endif
      end select
      !write(*,*) 'SPEC_SET_EN ',SPEC%FIX,ZI,ZF,EN,IERR
      if (ierr==0) then
        SPEC%KI=ZI
        SPEC%KF=ZF
        SPEC%EN=EN
        SPEC%LAMBDA=twopi/ZI
      endif
      end SUBROUTINE SPEC_SET_EN

!---------------------------------------------------------
      real(KIND(1.D0)) function SPEC_GET_EN(SPEC)
! setter for energy transfer
!---------------------------------------------------------
      TYPE(TSPEC),intent(in) :: SPEC
      SPEC_GET_EN=HSQOV2M*(SPEC%KI**2-SPEC%KF**2)
      end function SPEC_GET_EN

!---------------------------------------------------------
      SUBROUTINE SPEC_VALIDATE_POS(SPEC,IERR)
! call after change of input data to make all
! position parameters consistent
!---------------------------------------------------------
      TYPE(TSPEC),intent(inout) :: SPEC
      integer, intent(out) :: IERR
      IERR=0
      select case(SPEC%INPSET)
      ! diffraction input mode, set lambda, theta, EN=0
      case(0)
        SPEC%FIX=0
        SPEC%KI=TWOPI/SPEC%LAMBDA
        SPEC%KF=SPEC%KI
        SPEC%EN=0.D0
        call SPEC_SET_THETA(SPEC,SPEC%TTHETA)
      ! TAS input mode, adjust KI,KF for given EN, calc. TTHETA for given Q
      case(1)
        call SPEC_SET_EN(SPEC,SPEC%EN,IERR)
        if (ierr.ne.0) return
        call SPEC_GET_THETA(SPEC,SPEC%TTHETA,IERR)
        if (ierr.ne.0) return
      end select
      end SUBROUTINE SPEC_VALIDATE_POS



!---------------------------------------------------------
      SUBROUTINE SPEC_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSPEC) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      character(FIELD_BUFFER_STR) :: SARG
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2STR(ARG,0,SARG,LR)
        !  write(*,*) 'SPEC_INP ',trim(PNAME),' [',trim(ARG%DEF%TID),']','=',trim(SARG)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call SPEC_INP_R(OBJ,trim(ARG%ID),NUM,LR)
        case default
          write(*,*) 'SPEC_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SPEC_INP


!---------------------------------------------------------
      SUBROUTINE SPEC_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TSPEC) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call SPEC_OUT_R(OBJ,trim(ARG%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case default
          write(*,*) 'SPEC_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE SPEC_OUT

C---------------------------------------------------------
      SUBROUTINE SPEC_INP_R(SPEC,PNAME,ARG,NARG)
C parameter IO for instrument INTERFACE
C input data from ARG to OBJ for parameter namespace PNAME
! Setting of various instrument parameters
! Beamline components are not changed (this is done through INST_ADJUST)
! Only basic value checking is done here
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... TSPEC structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TSPEC) :: SPEC
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      integer :: IERR
      character(32) :: CNUM
      LR=1
      ierr=0
  !    write(*,*) 'SPEC_INP ',trim(PNAME),ARG(1)
      SELECT CASE (trim(PNAME))
! intrinsic variables
      case('INPSET'); SPEC%INPSET=NINT(ARG(1))
      case('KI')
        if (ARG(1).gt.0) then
        !  call SPEC_GET_EN(SPEC,EN)
          SPEC%KI=ARG(1)
        !  call SPEC_SET_EN(SPEC,EN,IERR)
        endif
      case('KF')
        if (ARG(1).gt.0) then
        !  call SPEC_GET_EN(SPEC,EN)
          SPEC%KF=ARG(1);
        !  call SPEC_SET_EN(SPEC,EN,IERR)
        endif
      case('Q0')
        if (ARG(1).gt.0) SPEC%Q=ARG(1)
      case('PSI')
        if ((ARG(1).ge.-90.D0).and.(ARG(1).le.90.D0)) then
          SPEC%PSI=ARG(1)*deg
        else
          ierr=1
        endif
      case('FIX')
       ! IAUX=SPEC%FIX
        SPEC%FIX=NINT(ARG(1))
        if (SPEC%FIX>1) SPEC%FIX=1
        if (SPEC%FIX<0) SPEC%FIX=0
      !  call SPEC_GET_EN(SPEC,EN)
      !  call SPEC_SET_EN(SPEC,EN,IERR)
      case('SS')
            if (ARG(1).le.0.D0) then
              SPEC%SS=-1
            else
              SPEC%SS=1
            endif
      case('FLATCONE');SPEC%FLATCONE=NINT(ARG(1))
      case('KFMODE');SPEC%KFMODE=NINT(ARG(1))
      case('OFFPL'); SPEC%OffPlane=(NINT(ARG(1)).eq.1)
      case('ADJ'); SPEC%ADJSPEC=(NINT(ARG(1)).eq.1)
      case('ORISAM'); SPEC%ORISAM=(NINT(ARG(1)).eq.1)
      case('QCONST'); SPEC%QCONST=(NINT(ARG(1)).eq.1)
! derived variables
      case('LAMBDA'); SPEC%LAMBDA=ARG(1)
        ! call SPEC_SET_LAMBDA(SPEC,ARG(1),IERR)
        ! write(*,*) 'SPEC_INP: LAMBDA=',ARG(1)
      case('THETA'); SPEC%TTHETA=ARG(1)*deg
        ! call SPEC_SET_THETA(SPEC,ARG(1)*deg,IERR)
    !    write(*,*) 'SPEC_INP: THETA=',ARG(1)
      case('EN'); SPEC%EN=ARG(1)
        ! call SPEC_SET_EN(SPEC,ARG(1),IERR)
      end select
      if (ierr.ne.0) then
        call FLOAT2STR(ARG(1),CNUM)
        call MSG_ERROR('SPEC_INP','unable to set '//trim(PNAME)//' = '//trim(CNUM),0,1)
        NARG=0
      else
        NARG=LR
      endif
      END SUBROUTINE SPEC_INP_R


C---------------------------------------------------------
      SUBROUTINE SPEC_OUT_R(SPEC,PNAME,ARG,NARG)
C input data from ARG to OBJ for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TSPEC),intent(in) :: SPEC
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
      case('INPSET'); ARG(1)=SPEC%INPSET
      case('KI'); ARG(1)=SPEC%KI
      case('KF');ARG(1)=SPEC%KF
      case('Q0');ARG(1)=SPEC%Q
      case('PSI');ARG(1)=SPEC%PSI/deg
      case('FIX');ARG(1)=SPEC%FIX
      case('SS')
            if (SPEC%SS<0) then
              ARG(1)=0.D0
            else
              ARG(1)=1.D0
            endif
      case('FLATCONE');ARG(1)=SPEC%FLATCONE
      case('KFMODE');ARG(1)=SPEC%KFMODE
! reference x-axis remains horizontal
      case('OFFPL'); ARG(1)=LOG2INT(SPEC%OffPlane)
! auto adjust
      case('ADJ'); ARG(1)=LOG2INT(SPEC%ADJSPEC)
      case('ORISAM'); ARG(1)=LOG2INT(SPEC%ORISAM)
      case('QCONST'); ARG(1)=LOG2INT(SPEC%QCONST)
      case('THETA'); ARG(1)=SPEC%TTHETA/deg
  !      write(*,*) 'SPEC_OUT: THETA=',Z/deg
      case('EN');ARG(1)=SPEC%EN ! HSQOV2M*(SPEC%KI**2-SPEC%KF**2)
      case('LAMBDA'); ARG(1)=SPEC%LAMBDA
      end select
      NARG=LR
      END SUBROUTINE SPEC_OUT_R


      end module SPECTROMETER

