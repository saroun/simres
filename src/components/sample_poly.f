!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.44 $
!////     $Date: 2019/08/15 15:02:07 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Describes component class: SAMPLE - POLYCRYSTAL
!////
!////////////////////////////////////////////////////////////////////////
      module SAMPLE_POLY
      use CONSTANTS
      use SAMPLES
      use TRACINGDATA
      USE FRAMES_TRACE
      use SAMPLE_POLY_TABLE
      use STRAIN_TABLE
      use RNDGEN
      implicit none

      save

      type TPCRYST; SEQUENCE
        TYPE(TSAMPLE) :: SAM            ! ancestor class
        real(KIND(1.D0)) :: FHKL        ! structure factor per unit volume [fm.A^-3]
        real(KIND(1.D0)) :: DHKL        ! d-spacing [A]
        real(KIND(1.D0)) :: DD          ! rel. spread of d-spacing distribution
        real(KIND(1.D0)) :: GRAIN       ! grain size [um]
        INTEGER :: JHKL                 ! plane multiplicity
! auxilliary parameters, not accessible by user input
        real(KIND(1.D0)) :: SSCAT       ! scattering cross-section [1/cm] for the diffraction cone
        real(KIND(1.D0)) :: THB         ! Bragg angle
        CHARACTER(256) :: REFTAB        ! filename with the table of reflections
        INTEGER :: ID_STRAIN_TABLE  ! >0 if there is a macroscopic strain defined
        integer :: IREF                 ! selected reflection from the table. If 0, use all
      end type TPCRYST

! global instance of TPCRYST
      TYPE(TPCRYST),target :: PCRYST

      logical,private :: dbg=.false.
      integer,private :: dbgn=0

      contains

C------------------------------------
      SUBROUTINE PCRYST_SETDBG(b)
C------------------------------------
      logical :: b
        dbg=b
        if (.not.b) dbgn=0
      END SUBROUTINE PCRYST_SETDBG


C------------------------------------
      SUBROUTINE PCRYST_DEFAULT(OBJ)
! default values for Fe 211
C------------------------------------
      TYPE (TPCRYST) :: OBJ
        call SAMPLE_DEFAULT(OBJ%SAM)
        OBJ%SAM%FRAME%CLASS=SCLS_POLY
        OBJ%SAM%SIGA=0.0121D0
        OBJ%SAM%SIGI=0.0034D0
        OBJ%FHKL=0.804D0
        OBJ%DHKL=1.17D0
        OBJ%DD=0.D0
        OBJ%GRAIN=0.D0 ! no size effect
        OBJ%JHKL=24    ! multiplicity for 211
        call CLEAR_REF_TABLE
        OBJ%REFTAB='none'
        OBJ%ID_STRAIN_TABLE=0
        OBJ%IREF=0
      END SUBROUTINE PCRYST_DEFAULT


!------------------------------------
      SUBROUTINE PCRYST_PREPARE(OBJ,IERR)
! prepare dependent fields after input
!------------------------------------
      type(TPCRYST) :: OBJ
      integer,intent(out) :: IERR
        IERR=0
        call FRAME_INIT_MAT(OBJ%SAM%FRAME)
      END SUBROUTINE PCRYST_PREPARE


C------------------------------------
      SUBROUTINE PCRYST_INIT(OBJ)
C------------------------------------
      TYPE (TPCRYST) :: OBJ
        call SAMPLE_INIT(OBJ%SAM)
        dbgn=0
      end SUBROUTINE PCRYST_INIT


!---------------------------------------------------------
      subroutine PCRYST_GETQE(OBJ,Q,E)
! return nominal Q,E values
!---------------------------------------------------------
      TYPE (TPCRYST) :: OBJ
      real(kind(1.D0)), intent(out) :: Q,E
      Q=twopi/OBJ%DHKL
      E=0.D0
      end subroutine PCRYST_GETQE


!-----------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION PCRYST_SSCAT(OBJ, LAMBDA)
! return total scattering cross-section for one reflection of an ideal polycrystal
!-----------------------------------------------------------------------
      TYPE (TPCRYST) :: OBJ
      real(kind(1.D0)),intent(in) :: LAMBDA
      real(kind(1.D0)) :: S
  ! slightly relax the Bragg edge condition in order to allow for strain effects
      if (LAMBDA<2.D0*OBJ%DHKL) then
        S=0.5D-3*OBJ%JHKL*OBJ%DHKL*(LAMBDA*OBJ%FHKL)**2
      else
        S=0.D0
      endif
      PCRYST_SSCAT=S
      end FUNCTION PCRYST_SSCAT

!-----------------------------------------------------------------
      SUBROUTINE PCRYST_SET_THETA(OBJ,THETA)
! Set scattering angle => adjust DHKL
!-----------------------------------------------------------------------
      TYPE (TPCRYST) :: OBJ
      real(kind(1.D0)),intent(in) :: THETA
      real(kind(1.D0)) :: LAMBDA,STHETA
      STHETA=OBJ%SAM%SGN*ABS(SIN(THETA/2.D0))
    !  OBJ%SAM%SGN=SIGN(1.D0,STHETA)
      if ((OBJ%SAM%KI0.gt.0).and.(ABS(STHETA).gt.0)) then
        LAMBDA=2.D0*PI/OBJ%SAM%KI0
        OBJ%DHKL=LAMBDA/2.D0/ABS(STHETA)
        OBJ%SAM%Q0=twopi/OBJ%DHKL
        OBJ%THB=asin(STHETA)
        OBJ%SSCAT=PCRYST_SSCAT(OBJ,LAMBDA)
      ENDIF
      end SUBROUTINE PCRYST_SET_THETA

!-----------------------------------------------------------------
      SUBROUTINE PCRYST_SET_DHKL(OBJ,DHKL)
! Set scattering angle => adjust DHKL
!-----------------------------------------------------------------------
      TYPE (TPCRYST) :: OBJ
      real(kind(1.D0)),intent(in) :: DHKL
      real(kind(1.D0)) :: LAMBDA,STHETA
      if (DHKL.GT.0) then
        OBJ%DHKL=DHKL
        OBJ%SAM%Q0=twopi/OBJ%DHKL
        if (OBJ%SAM%KI0.gt.0) then
          LAMBDA=2.D0*PI/OBJ%SAM%KI0
          STHETA=OBJ%SAM%SGN*LAMBDA/2.D0/DHKL
          if (ABS(STHETA).LT.1.D0) then
            OBJ%THB=asin(STHETA)
            OBJ%SSCAT=PCRYST_SSCAT(OBJ,LAMBDA)
          endif
        ENDIF
      endif
      end SUBROUTINE PCRYST_SET_DHKL


!-----------------------------------------------------------------
      SUBROUTINE PCRYST_DEPENDENCES(OBJ,IERR)
! Calculate variables dependent on SAM settings
!------------------------------------------- ----------------------
      TYPE (TPCRYST) :: OBJ
      integer,intent(out) :: ierr
      real(kind(1.D0)) :: STHETA,LAMBDA
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      character(32) :: CNUM
      ierr=1
      LAMBDA=twopi/OBJ%SAM%KI0
      if (LAMBDA.gt.0) then
        ierr=2
    ! scatt. cross-section in mm^-1
        OBJ%SSCAT=PCRYST_SSCAT(OBJ,LAMBDA)
    ! set exit axis angles
        STHETA=OBJ%SAM%SGN*LAMBDA/2.D0/OBJ%DHKL
        if (abs(STHETA).le.1.D0) then
          OBJ%THB=asin(STHETA)
          ierr=0
        endif
      endif
      select case(ierr)
      case(1)
        call MSG_ERROR('PCRYST_DEPENDENCES','Lambda must be > 0',0,1)
      case(2)
        call FLOAT2STR(LAMBDA,CNUM)
        call MSG_ERROR('PCRYST_DEPENDENCES','Can''t set Bragg angle for given wavelength '//trim(CNUM),0,1)
      end select
      end SUBROUTINE PCRYST_DEPENDENCES


!---------------------------------------------------------
      SUBROUTINE PCRYST_INP(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TPCRYST) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call FIELD2ARRAY(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
          call PCRYST_INP_R(OBJ,trim(ARG%DEF%ID),NUM,LR)
        case(ID_FIELD_STR)
          call FIELD2STR(ARG,0,SARG,LR)
          call PCRYST_INP_S(OBJ,trim(ARG%ID),SARG,LR)
        case default
          write(*,*) 'PCRYST_INP: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE PCRYST_INP

!---------------------------------------------------------
      SUBROUTINE PCRYST_OUT(OBJ,PNAME,ARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... TFIELD parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      type(TPCRYST) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: NARG
      integer :: LR,NA
      real(KIND(1.D0)) :: NUM(FIELD_BUFFER_FLOAT)
      CHARACTER(FIELD_BUFFER_STR) :: SARG
      LR=0
      NARG=0
      select case(ARG%DEF%TID)
        case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
          call PCRYST_OUT_R(OBJ,trim(ARG%DEF%ID),NUM,LR)
          call ARRAY2FIELD(ARG,0,NUM,FIELD_BUFFER_FLOAT,LR)
        case(ID_FIELD_STR)
          call PCRYST_OUT_S(OBJ,trim(ARG%ID),SARG,NA)
          call STR2FIELD(ARG,0,SARG,LR)
        case default
          write(*,*) 'PCRYST_OUT: undefined field type '//trim(ARG%DEF%TID)
      end select
      NARG=LR
      END SUBROUTINE PCRYST_OUT


!---------------------------------------------------------
      SUBROUTINE PCRYST_INP_S(OBJ,PNAME,SARG,NARG)
! input data from ARG to OBJ for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TPCRYST),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*),intent(in) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
        CASE('REFTAB')
          call PCRYST_READ_TAB(OBJ,trim(SARG))
      END SELECT
      NARG=LR
      END SUBROUTINE PCRYST_INP_S


C---------------------------------------------------------
      SUBROUTINE PCRYST_INP_R(OBJ,PNAME,ARG,NARG)
C input data from ARG to OBJ for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TPCRYST),intent(out) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(in) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! structure factor per unit volume [fm.A^-3]
        CASE('FHKL');   OBJ%FHKL=ARG(1)
! d-spacing [A]
        CASE('DHKL');   call PCRYST_SET_DHKL(OBJ,ARG(1))
! rel. dispersion of d-spacing
        CASE('DD');     OBJ%DD=ARG(1)
! grain size [um]
        CASE('GRAIN');  OBJ%GRAIN=ARG(1)*1.D4
! plane multiplicity
        CASE('JHKL');   OBJ%JHKL=NINT(ARG(1))
! scattering angle
        CASE('THETA');  call PCRYST_SET_THETA(OBJ,ARG(1)*deg)

! strain table index
        CASE('STRTAB'); OBJ%ID_STRAIN_TABLE = MAX(0,NINT(ARG(1)))
! reflection table index
        CASE('IREF'); OBJ%IREF = MAX(0,NINT(ARG(1)))
        CASE DEFAULT; CALL SAMPLE_INP_R(OBJ%SAM,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE PCRYST_INP_R


!---------------------------------------------------------
      SUBROUTINE PCRYST_OUT_S(OBJ,PNAME,SARG,NARG)
! output data from OBJ to ARG for parameter namespace PNAME
! INPUT
!    PNAME   .... parameter name
!    ARG     .... REAL*8 parameter values
! INPUT
!    OBJ     .... SLIT structure
!    NARG    .... number of items read from ARG
!---------------------------------------------------------
      TYPE(TPCRYST),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      CHARACTER(*), intent(out) :: SARG
      integer,intent(out) :: NARG
      integer :: LR
      LR=0
      SELECT CASE (trim(PNAME))
      ! backward compatibility requires TABNAME ...
        CASE('REFTAB')
          if (len_trim(OBJ%REFTAB).gt.0) SARG=trim(OBJ%REFTAB)
          LR=1
      END SELECT
      NARG=LR
      END SUBROUTINE PCRYST_OUT_S


C---------------------------------------------------------
      SUBROUTINE PCRYST_OUT_R(OBJ,PNAME,ARG,NARG)
C output data from OBJ to ARG for parameter namespace PNAME
C INPUT
C    PNAME   .... parameter name
C    ARG     .... REAL*8 parameter values
C INPUT
C    OBJ     .... SLIT structure
C    NARG    .... number of items read from ARG
C---------------------------------------------------------
      TYPE(TPCRYST),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      real(KIND(1.D0)), intent(out) :: ARG(*)
      integer,intent(out) :: NARG
      integer :: LR
      LR=1
      SELECT CASE (trim(PNAME))
! structure factor per unit volume [fm.A^-3]
        CASE('FHKL');   ARG(1)=OBJ%FHKL
! d-spacing [A]
        CASE('DHKL');   ARG(1)=OBJ%DHKL
! rel. dispersion of d-spacing
        CASE('DD');     ARG(1)=OBJ%DD
! grain size [um]
        CASE('GRAIN');  ARG(1)=OBJ%GRAIN*1.D-4
! plane multiplicity
        CASE('JHKL');   ARG(1)=OBJ%JHKL
! scattering angle
        CASE('THETA');   ARG(1)=2.D0*OBJ%THB/deg
! strain table index
        CASE('STRTAB');   ARG(1)=OBJ%ID_STRAIN_TABLE
! reflection table index
        CASE('IREF'); ARG(1)=OBJ%IREF
        CASE DEFAULT; CALL SAMPLE_OUT_R(OBJ%SAM,trim(PNAME),ARG(1),LR)
      END SELECT
      NARG=LR
      END SUBROUTINE PCRYST_OUT_R


C------------------------------------------------------
      LOGICAL*4 FUNCTION PCRYST_GO_OBSOLETE(OBJ,NEU)
C diffraction by polycrystal = scatter at |Q|=const,E=const
C Works in Q-coordinates
C------------------------------------------------------
      TYPE(TPCRYST),intent(in) :: OBJ
      TYPE(NEUTRON) :: NEU
      REAL(kind(1.D0)) :: stb,ctb,sinphi,cosphi,s2tb,c2tb,kf0,Q
      REAL(kind(1.D0)) :: TM(3,3),X1(3),PHI
1     format(a,': ',6(G10.4,1x))
      kf0=NEU%K0
      call getRotZX(NEU%K,(/1.D0,0.D0,0.D0/),TM)

! get Q-vector
  ! include strain broadening (Gaussian)
      Q=2*pi/OBJ%DHKL
      if (OBJ%DD.GT.0.D0) Q=Q*(1.D0-OBJ%DD*GASDEV()/R8LN2)
  ! include size effect (Lorentzian)
      if (OBJ%GRAIN.GT.0.D0) Q=Q+2.D0*pi/OBJ%GRAIN*tan(pi/2*(2*RAN1()-1))
! get Bragg angle
      stb=Q/2/kf0
      if (abs(stb).ge.1.D0) GOTO 99 ! Bragg edge
      ctb=SQRT(1.D0-stb**2)
      s2tb=2.D0*stb*ctb
      c2tb=ctb**2-stb**2
      !if (dbg) write(*,1) 'PCRYST_GO KI',NEU%K
      if (dbg) write(*,1) 'PCRYST_GO Q',Q
      !if (dbg) write(*,1) 'PCRYST_GO theta',atan2(stb,ctb)/deg
      !if (dbg) write(*,1) 'PCRYST_GO 2.theta',atan2(s2tb,c2tb)/deg
! get azimuthal angle
      PHI=XRND(IXRND(OBJ%SAM%IRND(2)))
      if (ABS(PHI).ge.PI) goto 99
      sinphi=sin(PHI)
      cosphi=cos(PHI)
      X1(1)=kf0*s2tb*cosphi
      X1(2)=kf0*s2tb*sinphi
      X1(3)=kf0*c2tb
      !if (dbg) write(*,1) 'PCRYST_GO phi',atan2(sinphi,cosphi)/deg
      !if (dbg) write(*,1) 'PCRYST_GO KF',X1
! transform it to the current coordinate system
      call M3xV3(-1,TM,X1,NEU%K)
      ! if (dbg) write(*,1) 'PCRYST_GO SSCAT',OBJ%SSCAT
      ! dbg=.false.
      NEU%P=NEU%P*OBJ%SSCAT/2.D0/PI
      PCRYST_GO_OBSOLETE=.true.
      RETURN
99    PCRYST_GO_OBSOLETE=.false.
      NEU%P=0.D0
      END FUNCTION PCRYST_GO_OBSOLETE

!----------------------------------------------------------------
      LOGICAL*4 FUNCTION PCRYST_GO_INST(OBJ,PATH,SIGMA)
! diffraction by polycrystal = scatter at |Q|=const,E=const
! Works in Q-coordinates
! input:
! PATH = path through the sample along incident beam
! Return:
! true if not absorbed
! SIGMA = total diffraction + absorption + incoherent cross-section [mm^-1]
! NOTE:
! simres integrates over thickness and azimuthal angle
! factor exp(-SIGMA*PATH) is added by the calling routine
!----------------------------------------------------------------
      TYPE(TPCRYST),intent(in) :: OBJ
      REAL(kind(1.D0)),intent(in) :: PATH
      REAL(kind(1.D0)),intent(out) :: SIGMA
      REAL(kind(1.D0)) :: Q
      REAL(kind(1.D0)) :: SC(0:PCRYST_MAXREF)
      integer :: ISC(1:PCRYST_MAXREF)
      REAL(kind(1.D0)) :: lambda,XR,PHI,P0, Z
      integer :: NS,i
1     format(a,': ',6(G13.6,1x))
      if (.not.(OBJ%SAM%TRANS.or.OBJ%SAM%SCATT)) goto 99
      lambda=TWOPI/NEUT%K0
  ! get azimuthal angle
      dbg=(dbg.and.(dbgn<20))
      dbgn=dbgn+1
      if (dbg) write(*,1) 'PCRYST_GO_INST lambda',lambda

  ! get scattering cross-sections for all reflections, cumulative
      SC(0)=0.D0
      ISC = 0
      NS=0
      ! from table
      if (PCRYST_NREF>0) then
        do i=1,PCRYST_NREF
          Z = REF_TABLE_SSCATT(i,lambda)
          if (Z>0.D0) then
            NS=NS+1
            SC(NS)=SC(NS-1)+Z
            ISC(NS) = i
            if (dbg) write(*,1) '    SIGMA',NS,REFTAB(1,NS), SC(NS)-SC(NS-1)
          endif
        enddo
        ! total removal cross-section (absorption + scattering)
        SIGMA=SC(NS)+REF_TABLE_SABS(lambda)
        if (dbg) write(*,1) '    SIGMA_TOT',SIGMA
      ! from input parameters
      else
        Z = PCRYST_SSCAT(OBJ,lambda)
        if (Z>0.D0) then
          NS=1
          SC(NS)=SC(NS-1)+Z
        endif
        ! total removal cross-section (absorption + scattering)
        SIGMA=SC(NS)+OBJ%SAM%SIGA*lambda+OBJ%SAM%SIGI
      endif

      if (dbg) write(*,1) 'PCRYST_GO_INST lambda,SIGMA',lambda,SC(NS)

  ! if only transmission is allowed:
  ! Simres does not integrate over sample thickness and phi if SCATT=false
  ! therefore no weight is necessary
      if (.not.OBJ%SAM%SCATT) goto 110

  ! clip PHI here, so that sampling limits can't exceed +-PI
      PHI=XRND(IXRND(OBJ%SAM%IRND(2)))
      if (dbg) write(*,1) '    PHI/PI',PHI/PI
      if (ABS(PHI).ge.PI) goto 99
      P0=1.D0
  ! both scattering and transmission:
  ! choose one of them with equal probability
      if (OBJ%SAM%TRANS.and.OBJ%SAM%SCATT) then
        P0=2.D0
        XR=RAN1()
        ! transmission
        if (XR<0.5) then
          NEUT%P=NEUT%P*P0/TWOPI/PATH
          goto 110
        endif
        ! no reflection, exit
        if (NS<=0) goto 99
        ! else continue with scattering
      endif
  ! random event selection
      XR=RAN1()*SC(NS)
  ! get Q
      if (PCRYST_NREF>0) then
        ! pre-selected single reflection from the table
        if ((OBJ%IREF>0).and.(OBJ%IREF<=PCRYST_NREF)) then
          i=1
          do while ((i<NS).and.(OBJ%IREF.ne.ISC(i)))
            i=i+1
          enddo
          if (OBJ%IREF.ne.ISC(i)) then
            goto 99
          endif
        else
        ! Russian rulette - choose reflection
          i=1
          do while ((XR>SC(i)).and.(i<NS))
            i=i+1
          enddo
        endif
        Q=TWOPI/REFTAB(1,i)
      else
        Q=TWOPI/OBJ%DHKL
      endif
  ! do scattering
      if (dbg) write(*,1) '    d',i,TWOPI/Q
      call PCRYST_DIFF(OBJ,Q,PHI)
! weight event
      NEUT%P=NEUT%P*P0*SC(NS)/TWOPI
      PCRYST_GO_INST=.true.
      !dbg=.false.
      RETURN
! transmission:
110   PCRYST_GO_INST=.true.
      !dbg=.false.
      RETURN
! absorption
99    PCRYST_GO_INST=.false.
      NEUT%P=0.D0
      !dbg=.false.
      END FUNCTION PCRYST_GO_INST


!----------------------------------------------------------------
      subroutine PCRYST_DIFF(OBJ,Q,PHI)
! diffraction for a single plane given by Q and PHI (azimuthal anngle)
! Works in Q-coordinates:    z//KI, y// (KI x Q0)
!----------------------------------------------------------------
      TYPE(TPCRYST),intent(in) :: OBJ
      REAL(kind(1.D0)),intent(in) :: Q,PHI
      REAL(kind(1.D0)) :: KI,stb,ctb,sinphi,cosphi,s2tb,c2tb
      REAL(kind(1.D0)) :: TM(3,3),TD(3,3),C(3),VQ(3),R(3),AUX(3),E(6),SP(3)
      REAL(kind(1.D0)) :: depth,dQ,Q0,dP
      integer :: i
1     format(a,': ',6(G13.6,1x))
      sinphi=sin(PHI)
      cosphi=cos(PHI)
      dQ=0.D0
      KI=NEUT%K0
  ! get conversion matrix from local to Z//NEU%K, Y//(NEU%K x (1,0,0)_loc)
      call getRotZX(NEUT%K,(/1.D0,0.D0,0.D0/),TM)
  ! include strain
      if (STRAIN_ISDEF(OBJ%ID_STRAIN_TABLE)) then
        stb=Q/2/KI
      ! get direction cosines of Q
        if (abs(stb).ge.1.D0) then  ! backscattering
          VQ(1)=0.D0
          VQ(2)=0.D0
          VQ(3)=-Q
        ELSE
          ctb=SQRT(1.D0-stb**2)
          s2tb=2.D0*stb*ctb
          c2tb=ctb**2-stb**2
        ! get Q vector
          VQ(1)=KI*s2tb*cosphi
          VQ(2)=KI*s2tb*sinphi
          VQ(3)=KI*(c2tb-1.D0)
        endif
        ! transform to local coordinates
        call M3xV3(-1,TM,VQ,AUX) ! to y // Ki x Q
        call M3xV3(-1,OBJ%SAM%QMAT,AUX,VQ) ! from y // Ki x Q to local       
        ! get depth under surface
        call M3xV3(-1,OBJ%SAM%QMAT,NEUT%R,R)
        ! depth = CONTAINER_DEPTH(OBJ%SAM%FRAME, R)
        call CONTAINER_COORD(OBJ%SAM%FRAME, R, depth, TD)
        !call PCRYST_SETDBG((OBJ%SAM%FRAME%COUNT<10))
        ! calculate direction cosines
        do i=1,3
          !C(i)=VQ(i)/Q
          C(i) = V3xV3(VQ, TD(:,i))/Q
        enddo
        ! get strain and probability
        call STRAIN_GET(dbg, OBJ%ID_STRAIN_TABLE,depth, E, SP )
        ! transform strain tensor elements to delta_Q 
        dQ=E(1)*C(1)**2+E(2)*C(2)**2+E(3)*C(3)**2
        dQ=dQ+E(4)*C(2)*C(3)+E(5)*C(1)*C(3)+E(6)*C(1)*C(2)
        dQ=Q*dQ
        dP = 1.D0 + SP(1)*abs(C(1))+SP(2)*abs(C(2))+SP(3)*abs(C(3))
        NEUT%P = NEUT%P*dP
        if (dbg) write(*,1) '    R',R
        if (dbg) write(*,1) '    depth',depth
		if (dbg) write(*,1) '    TD(1)',TD(:,1)
		if (dbg) write(*,1) '    TD(2)',TD(:,2)
		if (dbg) write(*,1) '    TD(3)',TD(:,3)
		if (dbg) write(*,1) '    Q',VQ/Q
        if (dbg) write(*,1) '    C',C
        if (dbg) write(*,1) '    E',E
        if (dbg) write(*,1) '    SP',SP
        if (dbg) write(*,1) '    dQ/Q, dP ',dQ/Q, dP
      endif
  ! include strain broadening (Gaussian)
      if (OBJ%DD.GT.0.D0) dQ=dQ+Q*OBJ%DD*GASDEV()/R8LN2
  ! include size effect (Lorentzian)
      if (OBJ%GRAIN.GT.0.D0) dQ=dQ+TWOPI/OBJ%GRAIN*tan(pi/2*(2*RAN1()-1.0))
      if (dbg) write(*,1) '    Q,dQ,KI',Q,dQ,KI
      Q0=Q-dQ
  ! get Bragg angle
      stb=Q0/2/KI
      if (dbg) write(*,1) '    sin(thb)',stb,Q0
      if (abs(stb).ge.1.D0) GOTO 99 ! Bragg edge
      ctb=SQRT(1.D0-stb**2)
      s2tb=2.D0*stb*ctb
      c2tb=ctb**2-stb**2
  ! get KF vector
      AUX(1)=KI*s2tb*cosphi
      AUX(2)=KI*s2tb*sinphi
      AUX(3)=KI*c2tb
      !if (dbg) write(*,1) '    K1',X1
  ! transform it to the current coordinate system
      call M3xV3(-1,TM,AUX,NEUT%K)
      if (dbg) write(*,1) '    K2',NEUT%K
      return
99    if (.not.OBJ%SAM%TRANS) then
        NEUT%P=0.D0
      endif
      END subroutine PCRYST_DIFF

!---------------------------------------------------------------
      SUBROUTINE PCRYST_READ_TAB(OBJ,FNAME)
! read source lookup table
!---------------------------------------------------------------
      TYPE(TPCRYST),intent(out) :: OBJ
      CHARACTER*(*),intent(in) :: FNAME
      integer :: ilin
      call READ_REF_TABLE(trim(FNAME),ilin)
      ! if (ilin.gt.0) call REPORT_FLUX_TABLE(ilin)
      OBJ%REFTAB=trim(REF_TABLE_NAME)
      end SUBROUTINE PCRYST_READ_TAB



      end module SAMPLE_POLY
