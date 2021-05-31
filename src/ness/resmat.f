!//////////////////////////////////////////////////////////////////////
!////  $Id: resmat.f,v 1.124 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.124 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for matrix description of neutron transport
!////  Theoretical background:
!////  see I. Ionita, A.D. Stoica, J. Appl. Cryst. (2000). 33, 1067-1074
!////
!//////////////////////////////////////////////////////////////////////
      MODULE RESMAT
  !    USE CONSTANTS
      use SIMATH
      USE TRACINGDATA
      USE COMPONENTS
      USE GENERATOR
      use FILETOOLS
      use XMLINFO
      use XTALS
      USE XTALS_REF
      USE SOURCES_TABLE
      implicit none

      real(KIND(1.D0)), parameter :: s_step=0.288675134595
      real(KIND(1.D0)), parameter :: s_tria=0.408248290464
      real(KIND(1.D0)), parameter :: s_circ=0.25D0
      real(KIND(1.D0)), parameter :: s_gaus=0.424660900144

      integer, parameter,private :: MD=128
      integer,private :: NTRM ! number of random variables
      integer,private :: NTRM_PRIM ! number of random variables for primary beamline


!C  TRM  ... transmission matrix
!C  CMAT ... coordinate conversion matrix
!C  CON  ... constraints
      real(KIND(1.D0)),private :: TRM(MD,MD),CMAT(MD,MD),CON(MD,MD)
  !    real(KIND(1.D0)) :: TRMUP(MD,MD),CMATUP(MD,MD),CONUP(MD,MD)

    ! partial matrices for primary spectrometer
    !  real(KIND(1.D0)),private :: TRM1(MD,MD),CON1(MD,MD)

!C Matrices created in RESMAT_EVAL for further use :
!C   CMATUP=CMAT^-1
!C   TRMUP .. transmission matrix for up-stream transport = CMATUP^T.TRM.CMATUP
!C   COVUP,COV  ... covariance matrices for down and up-stream (=TRM^-1 and TRMUP^-1)
!C   DIAG,DUP ... vectors with diagonalized COVUP and COV (eigenvalues)
!C   UMAT,UUP ... matrices diagonalizing COVUP and  COV (eigenvectors)
!C   DCOVUP,DETCOV ... determinants of COVUP, COV
      real(KIND(1.D0)),private :: COV(MD,MD),UMAT(MD,MD),DIAG(MD),DETCOV
  !    real(KIND(1.D0)) :: COVUP(MD,MD),UUP(MD,MD),DUP(MD),DCOVUP
  ! legend for matrix elements
      character*2048 MATLEG

      real(KIND(1.D0)),private :: KNOMINAL

    ! store shift of the component center (stage) for the 1st and last pairs of components
    ! it is used to calculate initial estimation of means for random event generator (TMEAN)
      real(KIND(1.D0)) :: XY_MEAN(2),KXY_MEAN(2),KZ_MEAN

  ! debug options
      logical, private :: DBG=.false.     ! flag
      integer,parameter,private :: IUDBG0=0   ! Log-file unit
      integer, private :: IUDBG=6   ! Log-file unit

      CONTAINS

!-----------------------------------------------------
!C Initialize arrays
!C Basic variable space is (x,y,z,kx,ky,kz)
!C All other variables (i>6) are those used by components
!-----------------------------------------------------
      subroutine RESMAT_INIT
      integer :: i
        TRM=0.D0
        CON=0.D0
        call UNIMAT(CMAT,MD,MD)
        NTRM=0
        NTRM_PRIM=0
        do i=1,CRND
          IXRND(i)=0
        enddo
        MATLEG=' '
        KNOMINAL=0.D0
        XY_MEAN=0.D0
        KXY_MEAN=0.D0
        DBG=(IUDBG0.gt.0)
        if (DBG) then
          IUDBG=IUDBG0
          if (IUDBG0.NE.6) IUDBG=OPENFILEUNIT(trim(OUTPATH)//'/RESMAT.log',.false.)
          DBG=(IUDBG.gt.0)
          if (DBG) write(*,*) 'RESMAT_INIT OK'
        endif
      end subroutine RESMAT_INIT


!-----------------------------------------------------------------------
! Evaluate matrices:
! Integrate over variables other than x,y,z,kx,ky,kz
! Calculate matrices required for further use (see file header)
!-----------------------------------------------------------------------
      subroutine RESMAT_EVAL
      INTEGER :: i
      REAL(KIND(1.D0)) :: AUX(MD,MD),AUX1(MD,MD),DUM
      INTEGER :: j,jmin
1     format(a,6(1x,G12.5))
      ! DBG=.false.
      if (DBG) call XML_RSXDUMP(IUDBG,' ',1)
! Constraints on z=0 at starting surface
! (either source or detector, depending on tracing direction).
      CON(3,3)=CON(3,3)+1.D0

! Set TOF mode on if there is a time-limited component
      write(*,1) 'Evaluating transmission matrix, dim=',NTRM
      if (TRM(7,7).gt.1.D-11) then
        if (TROPT%IMODE.ne.tr_secondary) IXRND(7)=7
        if (.not.TOF_MODE) then
          TOF_MODE=.true.
          call MSG_INFO('TOF mode is ON ',1)
        endif
! otherwise set constraint t=0
      else
        CON(7,7)=CON(7,7)+1.D0
        if (TOF_MODE) then
          TOF_MODE=.false.
          call MSG_INFO('TOF mode is OFF ',1)
        endif
      endif

! integrate over z
      call IntegrateVariable(3)

! integrate all other constrained variables

      jmin=5
      if (TROPT%IMODE.eq.tr_secondary) jmin=0
      j=NTRM
      do while (j.gt.jmin)
        call IntegrateVariable(j)
        j=j-1
      enddo

  ! check it !
      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT integrated ',MATLEG,MATLEG,CMAT(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM integrated ',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CON integrated ',MATLEG,MATLEG,CON(1,1),MD)
      endif

! COV
      if (DBG) write(*,*) 'RESMAT_EVAL  INVERT start'
      !write(*,*) 'RESMAT_EVAL  INVERT start ',NTRM
      J=INVERT(NTRM,TRM,MD,COV,MD)

      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV ',MATLEG,MATLEG,COV(1,1),MD)
      endif
      if (DBG) write(*,*) 'RESMAT_EVAL  INVERT OK'

! determinants
    !  DETCOV=DETNAP(COV,MD,NTRM)
      DETCOV=DETJAC(COV,MD,NTRM)
      if (DBG) write(*,*) 'RESMAT_EVAL  DET OK'


! DIAG,UMAT
      CALL JACOBI(COV,NTRM,MD,DIAG,UMAT,I)
      if (DBG) write(*,*) 'RESMAT_EVAL  JACOBI OK ',I

  !    if (TROPT%DIR.eq.tr_upstream) then
  !      call MxM(1,NTRM,MD,UMAT,AUX)
  !      UMAT=AUX
  !    endif

  ! check it !
      if (DBG) then
        call MxM(1,NTRM,MD,COV,UMAT,AUX)
        call MxM(-1,NTRM,MD,UMAT,AUX,AUX1)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV ',MATLEG,MATLEG,COV(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'UT.COV.U ',MATLEG,MATLEG,AUX1(1,1),MD)
      endif

      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV',MATLEG,MATLEG,COV(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'UMAT',MATLEG,MATLEG,UMAT(1,1),MD)
        call XML_ARRAY(IUDBG,NTRM,'DIAG',MATLEG,DIAG(1))
        DUM=1.D0
        do j=1,NTRM
          DUM=DUM*DIAG(j)
        enddo
        if (DETCOV.gt.0) write(IUDBG,*) 'RESMAT_EVAL  DET=',DETCOV,' SQRT(DET)=',SQRT(DETCOV)
        if (DUM.gt.0) write(IUDBG,*)    'RESMAT_EVAL  DET=',DUM,   ' SQRT(DET)=',SQRT(DUM)
      endif

        ! close log flie
      if (DBG) then
        call XML_RSXDUMP(IUDBG,' ',0)
        if (IUDBG.ne.6) close(IUDBG)
        IUDBG=6
        write(*,*) 'RESMAT_EVAL OK'
      endif
      write(*,1) 'Transmission matrix reduced to dim=',NTRM


      TRACING_VARID=MATLEG
      end subroutine RESMAT_EVAL


!-----------------------------------------------------------------------
!C Reverse tracing direction to up-stream, i.e. invert transport matrix etc.
!-----------------------------------------------------------------------
      subroutine RESMAT_REVERSE
      REAL(KIND(1.D0)) :: AUX(MD,MD)
      integer :: I
    ! INVERT works with a local copy of CMAT,
    ! so we can use the same variable for output and input
      I=INVERT2(NTRM,CMAT,MD)
    ! transform transmission matrix
      AUX=0.D0
      call BTAB(TRM,CMAT,AUX,NTRM)
      TRM=AUX
    ! transform constraint matrix
      AUX=0.D0
      call BTAB(CON,CMAT,AUX,NTRM)
      CON=AUX
      call UNIMAT(CMAT,MD,MD)
      end subroutine RESMAT_REVERSE

!--------------------------------------------------------------
!C Return eigenvalues and eigenvectors of the covariance matrix
!C does not evaluate! RESMAT_EVAL must be called before
!C   EIVALS  ... vector diagonalized covariance matrix (eigenvalues)
!C   EIVECT  ... matrix diagonalizing the covariance matrix (eigenvectors)
!C   IDIR     ... transport direction
!C   NL,N ... leading dimension and rank of the matrices
!C NOTE: the matrix rank, N is for output !
!---------------------------------------------------------------
      subroutine RESMAT_EIGEN(EIVECT,EIVALS,NL,N)
      integer, intent(in) :: NL
      REAL(KIND(1.D0)),intent(out) :: EIVALS(NL),EIVECT(NL,NL)
      integer, intent(out) :: N
      integer :: i,j
      do j=1,NTRM
        do i=1,NTRM
          EIVECT(i,j)=UMAT(i,j)
        enddo
        EIVALS(j)=DIAG(j)
      enddo
      N=NTRM
      end subroutine RESMAT_EIGEN

!--------------------------------------------------------
! Calculate initial mean position and momentum for the event generator
!---------------------------------------------------------
      subroutine RESMAT_EVAL_SHIFT
      TYPE(PCRYSTAL) :: CR
      TYPE(PXTAL) :: XT
      REAL(KIND(1.D0)) :: X1,Y1,DIST1,X2,Y2,DIST2
      type(TCOMOBJ) :: COMOBJ
1     format(a,': ',6(1x,G12.5))
      if (TROPT%DIR.eq.tr_upstream) then
        if (TROPT%IMODE.eq.tr_primary) then
          if (SAMOBJ%ICLS.gt.0) then
            ! X2=SAMOBJ%P_FRAME%STA(1)
            ! Y2=SAMOBJ%P_FRAME%STA(2)
            ! DIST2=SAMOBJ%P_FRAME%DIST
            ! call FRAME_GET_CTR(SAMOBJ%P_FRAME,X2,Y2,DIST2)
! The misfit correctioon does not work correctly e.g. for strain scanning.
! It's better to skip it in the case of sample misalignment.
! The matrix model needs improvement to include misfits in a genreal case
            X2=0.D0
            Y2=0.D0
          else
            COMOBJ=getComByIndex(BEAMLINE_NC(1),1)
            call GET_CTR_SHIFT(COMOBJ,X2,Y2,DIST2)
          endif
          !write(*,1) 'up, XY_MEAN',X2,Y2,XY_MEAN
          if (XY_MEAN(1).eq.0.D0) XY_MEAN(1)=X2
          if (XY_MEAN(2).eq.0.D0) XY_MEAN(2)=Y2
          if (BEAMLINE_NC(1)>1) then
            COMOBJ=getComByIndex(BEAMLINE_NC(1)-1,1)
            call GET_CTR_SHIFT(COMOBJ,X1,Y1,DIST1)
            if (DIST2>0.D0) then
              !write(*,1) 'up, KXY_MEAN',(X2-X1)/DIST2,(Y2-Y1)/DIST2,KXY_MEAN
              !write(*,1) '  is zero',(KXY_MEAN(1).ne.0.D0),(KXY_MEAN(2).ne.0.D0)
              if (KXY_MEAN(1).eq.0.D0) KXY_MEAN(1)=(X2-X1)/DIST2
              if (KXY_MEAN(2).eq.0.D0) KXY_MEAN(2)=(Y2-Y1)/DIST2
            endif
          endif
        endif
      else
        COMOBJ=getComByIndex(1,1)
        call GET_CTR_SHIFT(COMOBJ,X1,Y1,DIST1)
        !write(*,1) 'down, XY_MEAN',X2,Y2,XY_MEAN
        XY_MEAN=(/X1,Y1/)
        if (BEAMLINE_NC(1)>1) then
          COMOBJ=getComByIndex(2,1)
          call GET_CTR_SHIFT(COMOBJ,X2,Y2,DIST2)
          if (DIST2>0.D0) then
            !write(*,1) 'down,  KXY_MEAN',(X2-X1)/DIST2,(Y2-Y1)/DIST2,KXY_MEAN
            KXY_MEAN=(/X2-X1,Y2-Y1/)/DIST2
          endif
        endif
      endif
      KZ_MEAN=0.D0
      if (MONO_N.gt.0) then
        select case(MONOCHROMATORS(1)%CCLS)
        case(CCLS_CRYSTAL)
    ! CRYSTAL
          call CRYSTAL_GET(MONOCHROMATORS(1)%INST,CR)
          KZ_MEAN=(TWOPI/CR%X%LAMBDA-KNOMINAL)/KNOMINAL
        case(CCLS_XTAL)
    ! XTAL
          call XTAL_GET(MONOCHROMATORS(1)%INST,XT)
          KZ_MEAN=(TWOPI/XT%X%LAMBDA-KNOMINAL)/KNOMINAL
        end select
      endif

    ! KZ_MEAN does not work
    ! => don't call XTAL_ADJUST_LAMBDA in XTAL_ADJUST

      KZ_MEAN=0.D0


      !write(*,1) 'RESMAT_EVAL_SHIFT  XY',XY_MEAN(1:2)
      !write(*,1) 'RESMAT_EVAL_SHIFT KXYZ',KXY_MEAN(1:2),KZ_MEAN
      !write(*,1) 'RESMAT_EVAL_SHIFT LAMBDA',TWOPI/(KZ_MEAN*KNOMINAL+KNOMINAL)

      end subroutine RESMAT_EVAL_SHIFT

!--------------------------------------------------------
! Construct matrices for given beamline (primaroy or secondary)
! Primary beamline may include sample container
! Secondary beamline starts with sample scattering (without container constraints)
!---------------------------------------------------------
      subroutine RESMAT_ADDCOMPONENTS
      integer :: i,NT
      logical :: dbg1=.false.
      call RESMAT_INIT
      if (DBG) call XML_RSXDUMP(IUDBG,' ',1)
  ! PRIMARY
      call RESMAT_ADDGENERATOR ! does nothing unless limits are applied
      if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS RESMAT_ADDGENERATOR  OK '
      do i=1,BEAMLINE_NC(1)
        call RESMAT_ADDCOMPONENT(i,1)
        if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS 1 ',i
      enddo
  ! shift to sample incident plane
      if (SAMOBJ%ICLS.gt.0) then
        call ADD_SHIFT(SAMOBJ%P_FRAME%DIST)
      endif

    ! invert tracing direction for up-stream mode
      if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS RESMAT_REVERSE START'
      if (TROPT%DIR.eq.tr_upstream) then
        call RESMAT_REVERSE
      endif
      if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS RESMAT_REVERSE OK '

  ! add constraint given by container size/shape
      if (SAMOBJ%ICLS.gt.0) then
        if (SAMOBJ%P_SAMPLE%SCATT) then
          call RESMAT_CONTAINER(SAMOBJ%P_FRAME)
          if (dbg1) write(*,*) 'RESMAT_CONTAINER '
        endif
      endif

      NT=NTRM

  ! save matrices for primary beamline, incl. sample container
    !  TRM1=TRM
    !  CON1=CON
    !  NTRM1=NTRM
  ! SECONDARY (only when MODE<>primary)
    ! start with sample without container
      if (TROPT%IMODE.ne.tr_primary) then
        if (SAMOBJ%ICLS.gt.0) then
          if (SAMOBJ%P_SAMPLE%SCATT) then
            call RESMAT_ADDSAMPLE(SAMOBJ%P_SAMPLE)
            if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS sample '
          endif
        endif
        do i=1,BEAMLINE_NC(2)
          call RESMAT_ADDCOMPONENT(i,2)
          if (dbg1) write(*,*) 'RESMAT_ADDCOMPONENTS 2 ',i
        enddo
      endif

      if (TROPT%IMODE.eq.tr_secondary) then
        do i=1,NT
          CON(i,i)=CON(i,i)+1.D0
        enddo
        do i=1,CRND
          if (IXRND(i).le.NT) IXRND(i)=0
        enddo
      endif
      if (DBG) call XML_RSXDUMP(IUDBG,' ',0)

      call RESMAT_EVAL_SHIFT
      call RESMAT_EVAL
      if (dbg1)  write(*,*) 'RESMAT_ADDCOMPONENTS RESMAT_EVAL OK '
      end subroutine RESMAT_ADDCOMPONENTS

!--------------------------------------------------------
! Add a component with given component's and beamline indices
!---------------------------------------------------------
      subroutine RESMAT_ADDCOMPONENT(INDX,IBEAM)
      integer,intent(in) :: INDX,IBEAM
      type(TCOMOBJ) :: OBJ
      OBJ=getComByIndex(INDX,IBEAM)
      if (NTRM>MD) then
        call MSG_ERROR('RESMAT','Maximum dimension exceeded.',1,1)
      endif
      !write(*,*) 'RESMAT_ADDCOMPONENT ',trim(OBJ%ID)
      select case (OBJ%CCLS)
      case(CCLS_FRAME)
        call RESMAT_ADDSLIT(OBJ)
      case(CCLS_SOURCE)
        call RESMAT_ADDSOURCE(OBJ)
      case(CCLS_DETECTOR)
        call RESMAT_ADDSLIT(OBJ)
      case(CCLS_GUIDE)
        call RESMAT_ADDGUIDE(OBJ)
      case(CCLS_SGUIDE)
        call RESMAT_ADDSGUIDE(OBJ)
      case(CCLS_CRYSTAL)
        call RESMAT_ADDCRYSTAL(OBJ)
      case(CCLS_XTAL)
        call RESMAT_ADDXTAL(OBJ)
      case(CCLS_DCHOPPER)
        call RESMAT_ADDCHOPPER(OBJ)
      case(CCLS_MONITOR)
        call RESMAT_ADDMONITOR(OBJ,INDX,IBEAM)
      end select
      end subroutine RESMAT_ADDCOMPONENT


!-----------------------------------------------------
      subroutine ADD_SHIFT(DIST)
! Add shift along beam axis by DIST
!-----------------------------------------------------
      real(KIND(1.D0)),intent(in) :: DIST
      real(KIND(1.D0)) :: C(7,7)
      INTEGER :: INDX(7)
  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)
      call UNIMAT(C,7,7)
      C(1,4)=DIST
      C(2,5)=DIST
      C(7,6)=-DIST/KNOMINAL
 ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)
      end subroutine ADD_SHIFT


!-----------------------------------------------------
! Get the component center coordinates
! in incident axis frame.
!-----------------------------------------------------
      subroutine GET_CTR_SHIFT(COMOBJ,X,Y,DIST)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      real(KIND(1.D0)),intent(out) :: X,Y,DIST
      TYPE(TFRAME),pointer :: OBJ
      TYPE(PFRAME) :: POBJ
      call PCOMOBJ_GET(COMOBJ,POBJ)
      OBJ => POBJ%X
      call FRAME_GET_CTR(OBJ,X,Y,DIST)
      end subroutine GET_CTR_SHIFT

!-----------------------------------------------------
! Add an event generator.
! All beamline defition should start with this
!-----------------------------------------------------
      subroutine RESMAT_ADDGENERATOR
 !     real(KIND(1.D0)) :: T(7,7),S,X,Y,DK
 !     INTEGER :: INDX(7)
  ! local indices correspond to (x,y,z,kx,ky,kz,t)
 !     INDX=(/1,2,3,4,5,6,7/)
! transmission matrix in local coordinates

! define 1st 7 random variables (phase-space + time)
      MATLEG='x|y|z|kx|ky|kz|t'
      NTRM=7
      KNOMINAL=EVENTGEN%K0
 !     T=0.D0

  ! horizontal divergence
!      if (EVENTGEN%HPROF.NE.0) then
!        T(4,4)=1.D0/(s_step*EVENTGEN%DIVH)**2
!      endif
  ! vertical divergence
!      if (EVENTGEN%VPROF.NE.0) then
!        T(5,5)=1.D0/(s_step*EVENTGEN%DIVV)**2
!      endif
  ! wavelength
!      if (EVENTGEN%KPROF.NE.0) then
!        DK=2.D0*MAX(abs(EVENTGEN%KMAX-EVENTGEN%K0),abs(EVENTGEN%K0-EVENTGEN%KMIN))
!        T(6,6)=1.D0/(s_step*DK)**2
!      endif

! Add CMAT^T.T.CMAT
!      call AddBTAB(TRM,T,7,7,CMAT,MD,NTRM,INDX(1))

!      if (DBG) then
!        write(IUDBG,*) 'ADD GENERATOR: ',KNOMINAL
!        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        ! READ(*,*)
!      endif

      end subroutine RESMAT_ADDGENERATOR

!-----------------------------------------------------------------------
! Add a source component to the beamline
! similar to SLIT, but considers also limits of wavelength and divergence
!------------------------------------------------------------------------
      subroutine RESMAT_ADDSOURCE(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(TFRAME),pointer :: OBJ
      TYPE(PSOURCE) :: POBJ
      real(KIND(1.D0)) :: T(7,7),C(7,7),S,tmin,tmax,kmin,kmax,divx,divy,DK
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD)
      INTEGER :: INDX(7)
      call SOURCE_GET(COMOBJ%INST,POBJ)
      if (.not.associated(POBJ%X)) return
      OBJ => POBJ%X%FRAME
  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)

! pre-process
      call GetPreTrans(OBJ,C)
! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

! transmission matrix in local coordinates
      T=0.D0
    ! get range given by lookup table or Maxwell distr.
      call FLUX_RANGE(POBJ%X%TEMP,kmin,kmax)
      call FLUX_RANGE_DIV(divx,divy)
      select case(OBJ%SHAPE)
      case(FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER)
            S=s_step
      case(FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID)
            S=s_circ
      case default
            S=s_step
      end select
  ! width x height
      T(1,1)=1.D0/(S*OBJ%SIZE(1))**2
      T(2,2)=1.D0/(S*OBJ%SIZE(2))**2
  ! horizontal divergence
      if (divx.gt.0.D0) then
        T(4,4)=1.D0/(s_step*divx)**2
      endif
  ! vertical divergence
      if (divy.gt.0.D0) then
        T(5,5)=1.D0/(s_step*divy)**2
      endif
  ! wavelength band
      if (kmax.gt.kmin) then
        DK=2.D0*MAX(abs(kmax-EVENTGEN%K0),abs(EVENTGEN%K0-kmin))
        T(6,6)=1.D0/(s_step*DK)**2
      endif
  ! pulse width
      if (POBJ%X%TYP.eq.src_pulsed) then
        call SOURCE_PULSE_RANGE(POBJ%X,tmin,tmax)
        !write(*,*) 'RESMAT source tmin,tmax=',tmin,tmax
        if (tmax.gt.tmin) then
          T(7,7)=1.D0/(s_step*(tmax-tmin))**2
        endif
      endif
    ! k-vector spread
  !    call SOURCE_RANGE(k0,dkx,dky,dkz)
  !    if (dkz.gt.0) T(6,6)=1.D0/(s_step*dkz)**2
  !    if (dkx.gt.0) T(4,4)=1.D0/(s_step*dkx)**2
  !    if (dky.gt.0) T(5,5)=1.D0/(s_step*dky)**2
    ! include transition to z=0, where T referes to
      call GetTransZ0(OBJ,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,7,7,CT,MD,NTRM,INDX)

! post-process
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ,C)
! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

      if (DBG) then
      !  write(IUDBG,*) 'ADDSOURCE: ','class=',OBJ%CLASS,' ',trim(OBJ%NAME)
        call XML_VALUE(IUDBG,'ADD SOURCE','string',' ',NTRM,0.D0,trim(OBJ%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,7,7,'T',MATLEG,MATLEG,T(1,1),7)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif

      end subroutine RESMAT_ADDSOURCE

!-----------------------------------------------------
! Add an empty (transparent) component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDNONE(OBJ)
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)) :: C(7,7),A(3,3)
      INTEGER :: INDX(7)
  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)

  ! transition to by OBJ%DIST
      call UNIMAT(C,7,7)

      C(1,4)=OBJ%DIST
      C(2,5)=OBJ%DIST
      C(7,6)=-OBJ%DIST/KNOMINAL
      call UNIMAT(A,3,3)
      A(3,3)=0.0
      call MltpBySubmatrix(0,C,7,3,A,3,3,(/1,2,3/))
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

      if (DBG) then
        ! write(IUDBG,*) 'ADDNONE: ','class=',OBJ%CLASS,' ',trim(OBJ%NAME)
        call XML_VALUE(IUDBG,'ADD NONE','string',' ',NTRM,0.D0,trim(OBJ%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif
      end subroutine RESMAT_ADDNONE

!-----------------------------------------------------
! Add a slit component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDSLIT(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(TFRAME),pointer :: OBJ
      TYPE(PFRAME) :: POBJ
      real(KIND(1.D0)) :: T(2,2),C(7,7),S
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD)
      INTEGER :: INDX(7)

      call PCOMOBJ_GET(COMOBJ,POBJ)
      OBJ => POBJ%X

  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)

      if (OBJ%TRANSPARENT.EQ.1) then
        call RESMAT_ADDNONE(OBJ)
        return
      endif

! pre-process
      call GetPreTrans(OBJ,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

! transmission matrix in local coordinates
      select case(OBJ%SHAPE)
      case(FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER)
            S=s_step
      case(FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID)
            S=s_circ
      case default
            S=s_step
      end select
      T=0.D0
      T(1,1)=1.D0/(S*OBJ%SIZE(1))**2
      T(2,2)=1.D0/(S*OBJ%SIZE(2))**2
    ! include transition to z=0, where the slit limits apply
      call GetTransZ0(OBJ,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,2,2,CT,MD,NTRM,(/1,2/))

! post-process
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

      if (DBG) then
       ! write(IUDBG,*) 'ADDSLIT: ','class=',OBJ%CLASS,' ',trim(OBJ%NAME)
        call XML_VALUE(IUDBG,'ADD SLIT','string',' ',NTRM,0.D0,trim(OBJ%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,2,2,'T',MATLEG,MATLEG,T(1,1),2)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CT',MATLEG,MATLEG,CT(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif

      end subroutine RESMAT_ADDSLIT


!-----------------------------------------------------
! Add a slit component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDMONITOR(COMOBJ,COMIDX,IBEAM)
      integer,intent(in) :: COMIDX,IBEAM
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(MONITOR),pointer :: OBJ
      TYPE(PMONITOR) :: POBJ
      real(KIND(1.D0)) :: T(2,2),C(7,7),S
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD)
      INTEGER :: INDX(7)
      logical :: LOG1,isLast
      call MONITOR_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
1     format(a,': ',6(1x,G12.5))

  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)

! pre-process
      call GetPreTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

! transmission matrix in local coordinates
      select case(OBJ%FRAME%SHAPE)
      case(FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER)
            S=s_step
      case(FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID)
            S=s_circ
      CASE DEFAULT
            S=s_step
      end select
      T=0.D0
      T(1,1)=1.D0/(S*OBJ%FRAME%SIZE(1))**2
      T(2,2)=1.D0/(S*OBJ%FRAME%SIZE(2))**2

! include clip area
      ! set true if the monitor is the last component on the primary section
      ! and tracing direction is up-stream
      isLast=(TROPT%DIR.eq.tr_upstream)
      isLast=((IBEAM==1).and.(COMIDX==BEAMLINE_NC(1)))
      ! only for primary section
      isLast=(isLast.and.(TROPT%IMODE.eq.tr_primary))
      LOG1=(OBJ%BLOCKAREA.and.(OBJ%CNTMODE==monit_mode_inner))
      !write(*,1) 'RESMAT_ADDMONITOR',COMIDX,IBEAM,isLast,LOG1
      !if (LOG1) then

      ! disable this, it does not work well
      if (.FALSE.) then
        select case(OBJ%IX)
        case(1,2,4,5)
          T(OBJ%IX,OBJ%IX)=1.D0/(s_step*OBJ%DX)**2
          if (isLast) then
            select case(OBJ%IX)
            case(1,2)
              XY_MEAN(OBJ%IX)=OBJ%X0
            case(4,5)
              KXY_MEAN(OBJ%IX-3)=OBJ%X0
            end select
          endif
          !write(*,1) '   IX',OBJ%IX,OBJ%X0,XY_MEAN,KXY_MEAN
        end select
        select case(OBJ%IY)
        case(1,2,4,5)
          T(OBJ%IY,OBJ%IY)=1.D0/(s_step*OBJ%DY)**2
          if (isLast) then
            select case(OBJ%IY)
            case(1,2)
              XY_MEAN(OBJ%IY)=OBJ%Y0
            case(4,5)
              KXY_MEAN(OBJ%IY-3)=OBJ%Y0
            end select
          endif
          !write(*,1) '   IY',OBJ%IY,OBJ%Y0,XY_MEAN,KXY_MEAN
        end select
      endif
    ! include transition to z=0, where the slit limits apply
      call GetTransZ0(OBJ%FRAME,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,2,2,CT,MD,NTRM,(/1,2/))

! post-process
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

      if (DBG) then
       ! write(IUDBG,*) 'ADDSLIT: ','class=',OBJ%CLASS,' ',trim(OBJ%NAME)
        call XML_VALUE(IUDBG,'ADD MONITOR','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,2,2,'T',MATLEG,MATLEG,T(1,1),2)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CT',MATLEG,MATLEG,CT(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif

      end subroutine RESMAT_ADDMONITOR


!-----------------------------------------------------
! Add a chopper component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDCHOPPER(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(DCHOPPER),pointer :: OBJ
      TYPE(PDCHOPPER) :: POBJ
      real(KIND(1.D0)) :: T(3,3),C(7,7),S
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD),DT,DW
      INTEGER :: INDX(7),i
      call DCHOPPER_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)

! pre-process
      call GetPreTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)
! transmission matrix in local coordinates
      select case(OBJ%FRAME%SHAPE)
      case(FRAME_SHAPE_BOX,FRAME_SHAPE_CYLINDER)
            S=s_step
      case(FRAME_SHAPE_DISC,FRAME_SHAPE_ELLIPSOID)
            S=s_circ
      CASE DEFAULT
            S=s_step
      end select
      T=0.D0
      T(1,1)=1.D0/(S*OBJ%FRAME%SIZE(1))**2
      T(2,2)=1.D0/(S*OBJ%FRAME%SIZE(2))**2
      ! no time constraint for frame overlappin chopper
      if (.not.OBJ%OVERLAP) then
        DT=0.0
        do i=1,OBJ%NW
          DW=abs(OBJ%PHASES(i))+abs(OBJ%WIN*OBJ%WIDTHS(i)/2.D0)
          DT=MAX(DT,2.D0*DW/abs(OBJ%FRQ)*HOVM)
        enddo
        T(3,3)=1.D0/(s_step*DT)**2
      endif
! include transition to z=0, where the slit limits apply
      call GetTransZ0(OBJ%FRAME,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,3,3,CT,MD,NTRM,(/1,2,7/))
! post-process
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)
      if (DBG) then
       ! write(IUDBG,*) 'ADDCHOPPER: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD DCHOPPER','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,3,3,'T',MATLEG,MATLEG,T(1,1),3)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CT',MATLEG,MATLEG,CT(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif
      end subroutine RESMAT_ADDCHOPPER



!-----------------------------------------------------
! Add constraint given by container size/shape
! Adds penetration depth as a new random variable.
! Does not call post-process transformations,
! = ends in local coordinates !!
!-----------------------------------------------------
      subroutine RESMAT_CONTAINER(OBJ)
      TYPE(TFRAME),intent(out) :: OBJ
      real(KIND(1.D0)) :: T(3,3),C(7,7),SX,SY,SZ,K0(3)
      INTEGER :: INDX(7),i
  ! local indices correspond to (x,y,z,kx,ky,kz,t)
      INDX=(/1,2,3,4,5,6,7/)
! pre-process
      call GetPreTrans(OBJ,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX)

! Add penetration depth, PD
      NTRM=NTRM+1
      MATLEG=trim(MATLEG)//'|pd'

    ! OBJ must know where to find corresponding random variable in the XRND array:
    ! it will be obtained as XRND(IXRND(OBJ%FRAME%IRNDP))
      IXRND(NTRM)=NTRM
      OBJ%IRNDP=NTRM

! add transition matrix along k by PD
      CALL M3XV3(OBJ%MLOC,OBJ%RLOC,(/0.D0,0.D0,1.D0/),K0)
      do i=1,3
        CMAT(i,NTRM)=K0(i)
      enddo

! transmission matrix in local coordinates
      select case(OBJ%SHAPE)
      case(FRAME_SHAPE_BOX)
        SX=s_step
        SY=s_step
        SZ=s_step
      case(FRAME_SHAPE_CYLINDER)
        SX=s_circ
        SY=s_step
        SZ=s_circ
      case(FRAME_SHAPE_DISC)
        SX=s_circ
        SY=s_circ
        SZ=s_step
      case(FRAME_SHAPE_ELLIPSOID)
        SX=s_circ
        SY=s_circ
        SZ=s_circ
      case default
        SX=s_step
        SY=s_step
        SZ=s_step
      end select
      T=0.D0
      T(1,1)=1.D0/(SX*OBJ%SIZE(1))**2
      T(2,2)=1.D0/(SY*OBJ%SIZE(2))**2
      T(3,3)=1.D0/(SZ*OBJ%SIZE(3))**2

! Update transmission matrix
      call AddBTAB(TRM,T,3,3,CMAT,MD,NTRM,(/1,2,3/))

! post-process
  ! no post-transformation for this subroutine !

      if (DBG) then
       ! write(IUDBG,*) 'ADDFRAME: ','class=',OBJ%CLASS,' ',trim(OBJ%NAME)
        call XML_VALUE(IUDBG,'ADD FRAME','string',' ',NTRM,0.D0,trim(OBJ%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif

      end subroutine RESMAT_CONTAINER


!-----------------------------------------------------------------------
! Add sample to the beamline
! ATTENTION: Starts with local coordinates at the entry. Pre-transformation
! should be done by preceding call to ADDCONTAINER !!
! added random variables are X=(theta, phi, dE)
! pd    ... penetration depth (along incident beam)
! theta ... horizontal scattering angle
! phi   ... vertical scattering angle
! dE    ... energy transfer (positive for down scattering)
! theta and phi are measured w.r.t. EXIT coordinates
!------------------------------------------------------------------------
      subroutine RESMAT_ADDSAMPLE(OBJ)
      TYPE(TSAMPLE) :: OBJ
      real(KIND(1.D0)) :: KI0,KF0
      real(KIND(1.D0)) :: MB(6,6),SC(6,6),B(6),Cout(3,3),U(3,3),R(3,3),R1(3,3),T(2,2)
      real(KIND(1.D0)) :: KI(3),KF(3),G(3),AUX(3),AUX1(3),EXC(3),Ef
      INTEGER :: INDX(7),i,j,NX
      INDX=(/1,2,3,4,5,6,7/)


      KI0=OBJ%KI0
      KF0=OBJ%KF0
      KNOMINAL=KI0

  ! conversion to Q-coordinates: z//KI, y // (KI x Q)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,OBJ%QMAT,3,3,(/1,2,3/))
      call MltpBySubmatrix(0,CMAT,MD,NTRM,OBJ%QMAT,3,3,(/4,5,6/))
  ! get also nominal KI KF in Q-coordinates
    !  CALL M3XV3(1,OBJ%QMAT,(/0.D0,0.D0,KI0/),KI)
    !  CALL M3XV3(-OBJ%FRAME%MEXI,OBJ%FRAME%REXI,(/0.D0,0.D0,KF0/),AUX)
    !  CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,AUX,AUX1)
    !  CALL M3XV3(1,OBJ%QMAT,AUX1,KF)
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,KI0/),AUX)
      AUX1=AUX+OBJ%QSC
      CALL M3XV3(1,OBJ%QMAT,AUX1,KF)
      CALL M3XV3(1,OBJ%QMAT,AUX,KI)

  ! transform OBJ%EXC from local to Q-coordinates and normalize
      CALL M3XV3(1,OBJ%QMAT,OBJ%EXC,EXC)
      EXC=2.D0*EXC/HSQOV2M
      Ef=HSQOV2M*KF0**2

  ! add scattering angle horizontal
      INDX=(/4,5,6,0,0,0,0/)
      NX=4
      NTRM=NTRM+1
      INDX(NX)=NTRM
      MATLEG=trim(MATLEG)//'|theta'
  ! add scattering angle vertical
      NTRM=NTRM+1
      NX=NX+1
      INDX(NX)=NTRM
      MATLEG=trim(MATLEG)//'|phi'
  ! add energy transfer
      NTRM=NTRM+1
      NX=NX+1
      INDX(NX)=NTRM
      MATLEG=trim(MATLEG)//'|dE'
  ! make relevant random variables non-integrable (used by SAMPLE tracing procedure)
      select case(OBJ%TYP)
      case(isam_INELASTIC)
        do i=1,3
          IXRND(NTRM-3+i)=NTRM-3+i
          OBJ%IRND(i)=NTRM-3+i
        enddo
      case(isam_ELASTIC,isam_PHONON,isam_VANAD)
        do i=1,2
          IXRND(NTRM-3+i)=NTRM-3+i
          OBJ%IRND(i)=NTRM-3+i
        enddo
      case(isam_POWDER)
        IXRND(NTRM-1)=NTRM-1
        OBJ%IRND(2)=NTRM-1
      end select

! impose limits on scattering angles
      T=0.D0
      T(1,1)=1.D0/(s_step*PI)**2
      T(2,2)=1.D0/(s_step*PI)**2
! Update transmission matrix
      call AddBTAB(TRM,T,2,2,CMAT,MD,NTRM,(/NTRM-2,NTRM-1/))

! Add scattering law as a constraint
! given in Q-coordinates, i.e. KI(1)=0,KI(2)=0,KI(3)=KI0
      B=0.D0
      select case(OBJ%TYP)
      case(isam_ELASTIC,isam_VANAD)
        B(6)=1.D0
  ! powder diffraction
  ! G0.(dKf - dKi) = 0
      case(isam_POWDER)
        G=KF-KI
        ! B(1)=G(1)-G(3)*KF(1)
        ! B(3)=-G(1)*KI0
        ! B(4)=G(1)*KF(1)+G(3)*KF(3)-G(3)*KI0
        ! bug ? correction:
        B(1)=-G(1)*KI0
        B(3)=G(1)*KF(1)+G(3)*KF(3)-G(3)*KI0
        B(4)=G(1)*KF(3)-G(3)*KF(1)
        !------------
        B(6)=1.D0
  ! inelastic coherent (phonon dispersion)
  ! (KI-EXC).dKi - (KF-EXC).dKf = 0
      case(isam_PHONON)
        AUX=KF-EXC
        B(1)=-EXC(1)*KI0
        B(2)=-EXC(2)*KI0
        B(3)=AUX(3)*KI0  - (KI0/KF0)**2*(AUX(1)*KF(1)+AUX(3)*KF(3))
        B(4)=-(AUX(1)*KF(3)-AUX(3)*KF(1))
        B(5)=EXC(2)*KF(1)
        B(6)=(AUX(1)*KF(1)+AUX(3)*KF(3))/2.D0/Ef
      end select
    ! create constraint matrix
      do j=1,6
        DO i=j,6
          MB(i,j)=B(i)*B(j)
          MB(j,i)=MB(i,j)
        enddo
      enddo
! Update constraint matrix
      call AddBTAB(CON,MB,6,NX,CMAT,MD,NTRM,INDX(1))

  ! add additional kf component due to theta,phi,dE  variables
    !  call UNIMAT(SC,6,6)
      SC=0.D0
      SC(1,3)=KF(1)*(KI0/KF0)**2/KF0
      SC(1,4)=KF(3)/KF0
      SC(1,6)=-0.5D0/Ef*KF(1)/KF0
      SC(3,3)=KF(3)*(KI0/KF0)**2/KF0
      SC(3,4)=-KF(1)/KF0
      SC(3,6)=-0.5D0/Ef*KF(3)/KF0

      ! SC(2,5)=+KF(1)/KF0
      select case(OBJ%COORD)
      case(csam_PWD)
        SC(2,5)=+KF(1)/KF0
      case(csam_SANS)
        SC(2,5)=1.D0
      end select

      call MltpBySubmatrix(0,CMAT,MD,NTRM,SC,6,NX,INDX(1))
  ! CMAT now transforms ki->kf in Q-coordinates
  ! transform to exit coordinates by REXI.RLOC^T.QMAT^-1
      call UNIMAT(U,3,3)
      Call M3XM3(-1,OBJ%QMAT,U,R)
      Call M3XM3(-OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,R,R1)
      Call M3XM3(OBJ%FRAME%MEXI,OBJ%FRAME%REXI,R1,Cout)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cout,3,3,(/1,2,3/))
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cout,3,3,(/4,5,6/))
  ! change KNOMINAL to KF0
      KNOMINAL=OBJ%KF0
      if (DBG) then
      !  write(IUDBG,*) 'ADDSAMPLE: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD SAMPLE','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
      !  write(IUDBG,*) 'IRND =',OBJ%IRND(1:3)
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
      !  read(*,*)
      endif

      end subroutine RESMAT_ADDSAMPLE


!-----------------------------------------------------
! Add a guide component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDGUIDE(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(GUIDE),pointer :: OBJ
      TYPE(PGUIDE) :: POBJ
      integer :: INDX(9),NX,NH,NV
      real(KIND(1.D0)) :: T(9,9),C(7,7)
      real(KIND(1.D0)) :: beta,nc,z,GH,GV,LAMBDA
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD)
      call GUIDE_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
! empty component: transition only
      if (OBJ%FRAME%TRANSPARENT.EQ.1) then
        call RESMAT_ADDNONE(OBJ%FRAME)
        return
      endif

  ! local indices correspond to (x,y,z,kx,ky,kz,t,gh,gv)
  ! gh,gv are random reflection angles (hor. and vert.)
      INDX=(/1,2,3,4,5,6,7,0,0/)
      NX=7
      GH=0.D0
      GV=0.D0
      NH=0
      NV=0
! pre-process
      call GetPreTrans(OBJ%FRAME,C)
! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))
! transmission matrix in local coordinates
      T=0.D0
  ! horizontal
      beta=OBJ%W2/OBJ%FRAME%SIZE(1)
      nc=OBJ%NLH-1
      z=(nc+1)**2/(s_step*OBJ%W2)**2/(2*nc**2+1.0)
      T(1,1)=z*(((beta-1)*nc)**2+beta**2+1.0)
      T(1,4)=z*(1-(beta-1)*nc**2)*OBJ%FRAME%SIZE(3)
      T(4,1)=T(1,4)
      T(4,4)=z*(1+nc**2)*OBJ%FRAME%SIZE(3)**2
  ! vertical
      beta=OBJ%H2/OBJ%FRAME%SIZE(2)
      nc=OBJ%NLV-1
      z=(nc+1)**2/(s_step*OBJ%H2)**2/(2*nc**2+1)
      T(2,2)=z*(((beta-1)*nc)**2+beta**2+1.0)
      T(2,5)=z*(1-(beta-1)*nc**2)*OBJ%FRAME%SIZE(3)
      T(5,2)=T(2,5)
      T(5,5)=z*(1+nc**2)*OBJ%FRAME%SIZE(3)**2

  ! reflecting guide/bender
  ! a rough approximation, increasing bragg width for kx,ky
      LAMBDA=2.D0*PI/KNOMINAL
      if (OBJ%TYP.GT.0) then
        GH=LAMBDA*SQRT((OBJ%GHLU**2+OBJ%GHLD**2+OBJ%GHRU**2+OBJ%GHRD**2)/4.D0)
        GV=LAMBDA*SQRT((OBJ%GVT**2+OBJ%GVB**2)/2.D0)
      !  write(*,*) 'ADDGUIDE LAMBDA= ',LAMBDA
      !  write(*,*) 'ADDGUIDE 1: ',T(4,4),GH,GV
    ! horizontal reflections
        if (GH.GT.0) then
          if (OBJ%FRAME%SIZE(3).gt.0.D0) GH=GH+ABS(OBJ%FRAME%SIZE(1)-OBJ%W2)/OBJ%FRAME%SIZE(3)
          NTRM=NTRM+1
          NX = NX+1
          NH = NX
          INDX(NX)=NTRM
          T(NX,NX)=(1.D0/s_step/GH)**2
          T(4,4)=1.D0/(1.D0/T(NX,NX)+1.D0/T(4,4))
    !      write(*,*) 'ADDGUIDE 2: ',T(4,4),T(NX,NX)
          MATLEG=trim(MATLEG)//'|gh'
        endif
    ! vertical reflections
        if (GV.GT.0) then
          if (OBJ%FRAME%SIZE(3).gt.0.D0) GV=GV+ABS(OBJ%FRAME%SIZE(2)-OBJ%H2)/OBJ%FRAME%SIZE(3)
          NTRM=NTRM+1
          NX = NX+1
          NV = NX
          INDX(NX)=NTRM
          T(NX,NX)=(1.D0/s_step/GV)**2
          T(5,5)=1.D0/(1.D0/T(NX,NX)+1.D0/T(5,5))
          MATLEG=trim(MATLEG)//'|gv'
        endif
      endif

    ! include transition to z=0, where T referes to
      call GetTransZ0(OBJ%FRAME,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,9,NX,CT,MD,NTRM,INDX(1))

! post-process
    ! reflecting guide/bender: add random reflection at the guide end
      if (OBJ%TYP.GT.0) then
      ! horizontal reflections
        if (GH.GT.0) then
          CMAT(4,INDX(NH))=CMAT(4,INDX(NH))+1.D0
          CMAT(1,INDX(NH))=CMAT(1,INDX(NH))-OBJ%FRAME%SIZE(3)
        endif
      ! vertical reflections
        if (GV.GT.0) then
          CMAT(5,INDX(NV))=CMAT(5,INDX(NV))+1.D0
          CMAT(2,INDX(NV))=CMAT(2,INDX(NV))-OBJ%FRAME%SIZE(3)
        endif
      endif
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))

      if (DBG) then
      !  write(IUDBG,*) 'ADDGUIDE: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD GUIDE','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,7,7,'T',MATLEG,MATLEG,T(1,1),9)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
      !  READ(*,*)
      endif
      end subroutine RESMAT_ADDGUIDE


!-----------------------------------------------------
! Add a guide component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDSGUIDE(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(SGUIDE),pointer :: OBJ
      TYPE(PSGUIDE) :: POBJ
      integer :: INDX(9),NX,i,NH, NV
      real(KIND(1.D0)) :: T(9,9),C(7,7)
      real(KIND(1.D0)) :: beta,z,GH,GV,LAMBDA,GL,W,H
      real(KIND(1.D0)) :: C4(4,4),CT(MD,MD)
      call SGUIDE_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
! empty component: transition only
      if (OBJ%FRAME%TRANSPARENT.EQ.1) then
        call RESMAT_ADDNONE(OBJ%FRAME)
        return
      endif

  ! local indices correspond to (x,y,z,kx,ky,kz,t,gh,gv)
  ! gh,gv are random reflection angles (hor. and vert.)
      INDX=(/1,2,3,4,5,6,7,0,0/)
      NX=7
      GH=0.D0
      GV=0.D0
      NH=0
      NV=0
! pre-process
      call GetPreTrans(OBJ%FRAME,C)
! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))
! transmission matrix in local coordinates
      T=0.D0
  ! horizontal
      W=OBJ%W(OBJ%NSEG)      ! width of the exit window
      GL=OBJ%WIN(3,OBJ%NSEG) ! total length
      beta=W/OBJ%FRAME%SIZE(1)
      z=1.D0/(s_step*W)**2
      T(1,1)=z*(beta**2+1.0)
      T(1,4)=z*GL
      T(4,1)=T(1,4)
      T(4,4)=z*GL**2
  ! vertical
      H=OBJ%H(OBJ%NSEG)      ! height of the exit window
      beta=H/OBJ%FRAME%SIZE(2)
      z=1.D0/(s_step*H)**2
      T(2,2)=z*(beta**2+1.0)
      T(2,5)=z*GL
      T(5,2)=T(2,5)
      T(5,5)=z*GL**2

  ! reflecting guide/bender
  ! a rough approximation, increasing bragg width for kx,ky
      LAMBDA=2.D0*PI/KNOMINAL
      GH=0.0
      GV=0.0
      do i=1,OBJ%NSEG
        GH=GH+OBJ%MC(1,i)**2+OBJ%MC(2,i)**2
        GV=GV+OBJ%MC(3,i)**2+OBJ%MC(4,i)**2
      enddo
    ! horizontal reflections
      if (GH.GT.0) then
        GH=LAMBDA*SQRT(GH/OBJ%NSEG)*GammaNi
        if (GL.gt.0.D0) GH=GH+ABS(OBJ%FRAME%SIZE(1)-W)/GL
        NTRM=NTRM+1
        NX=NX+1
        NH=NX
        INDX(NX)=NTRM
        T(NX,NX)=(1.D0/s_step/GH)**2
        T(4,4)=1.D0/(1.D0/T(NX,NX)+1.D0/T(4,4))
    !      write(*,*) 'ADDGUIDE 2: ',T(4,4),T(8,8)
        MATLEG=trim(MATLEG)//'|gh'
      endif
    ! vertical reflections
      if (GV.GT.0) then
        GV=LAMBDA*SQRT(GV/OBJ%NSEG)*GammaNi
        if (GL.gt.0.D0) GV=GV+ABS(OBJ%FRAME%SIZE(2)-H)/GL
        NTRM=NTRM+1
        NX=NX+1
        NV=NX
        INDX(NX)=NTRM
        T(NX,NX)=(1.D0/s_step/GV)**2
        T(5,5)=1.D0/(1.D0/T(NX,NX)+1.D0/T(5,5))
        MATLEG=trim(MATLEG)//'|gv'
      endif
    ! include transition to z=0, where T referes to
      call GetTransZ0(OBJ%FRAME,C4)
      CT=CMAT
      call MltpBySubmatrix(0,CT,MD,NTRM,C4,4,4,(/1,2,3,7/))
! Update transmission matrix
      call AddBTAB(TRM,T,9,NX,CT,MD,NTRM,INDX(1))

! post-process
    ! reflecting guide/bender: add random reflection at the guide end
      ! horizontal reflections
      if (GH.GT.0) then
        CMAT(4,INDX(NH))=CMAT(4,INDX(NH))+1.D0
        CMAT(1,INDX(NH))=CMAT(1,INDX(NH))-GL
      endif
      ! vertical reflections
      if (GV.GT.0) then
        CMAT(5,INDX(NV))=CMAT(5,INDX(NV))+1.D0
        CMAT(2,INDX(NV))=CMAT(2,INDX(NV))-GL
      endif

  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
  ! update CMAT, CMAT=C.CMAT
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))

      if (DBG) then
      !  write(IUDBG,*) 'ADDGUIDE: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD GUIDE','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,7,7,'T',MATLEG,MATLEG,T(1,1),9)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
      !  READ(*,*)
      endif
      end subroutine RESMAT_ADDSGUIDE


!-----------------------------------------------------
! Add a crystal component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDCRYSTAL(COMOBJ)
      use CRYSTALS
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(PCRYSTAL) :: POBJ
      TYPE(CRYSTAL),pointer :: OBJ
      integer :: i,j,INDX(10),NX
      real(KIND(1.D0)) :: T(10,10),C(7,7),Cin(5,5),Cout(5,5),Cg(10,10),B(10), MB(10,10), DG(3,3)
      real(KIND(1.D0)) :: K0(3),KG(3),KABS,KSQR,cosG
1     format(a,': ',10(G12.5,1x))

      call CRYSTAL_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
! empty component: transition only
      if ((OBJ%TYP.eq.ctyp_filter).or.(OBJ%FRAME%TRANSPARENT.EQ.1)) then
        call RESMAT_ADDNONE(OBJ%FRAME)
        return
      endif

  ! local indices correspond to (x,y,z,kx,ky,kz,t,PD,etah,etav)
  ! PD is the penetration depth, measured as time of flight in units of h_bar/m.
  ! ETAH is the horizontal mosaic angle
  ! ETAV is the vertical mosaic angle
      INDX=(/1,2,3,4,5,6,7,0,0,0/)
      NX=7 ! maximum index actually used

! update CMAT by pre-process matrix
      call GetPreTrans(OBJ%FRAME,C)
      if (DBG) call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT 1',MATLEG,MATLEG,CMAT(1,1),MD)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))
      if (DBG) call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT 2',MATLEG,MATLEG,CMAT(1,1),MD)

  ! get nominal incident k-vector K0 in local coordinates and KK=|K0|^2

  ! get k-vector in local coordinates
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,1.D0/),K0)
    ! projection cosine
      cosG=-v3xv3(K0,OBJ%G)/OBJ%GTOT
    ! |k| value
      KABS=OBJ%GTOT/2.D0/abs(cosG)
  ! get nominal vectors K and KG=K0+G
      do i=1,3
        K0(i)=K0(i)*KABS
        KG(i)=K0(i)+OBJ%G(i)
      enddo
      KSQR=KABS**2
      ! KNOMINAL=KABS

! Add penetration depth, PD
! For a bent perfect crystal, it is actually not a variable
! but it will be defined through the Bragg constraint later
      NTRM=NTRM+1
      NX=8
      INDX(NX)=NTRM ! index of PD variable
      MATLEG=trim(MATLEG)//'|pd'

  ! for MODEL=0 mosaic crystals only (no multiple relections):
    ! PD is included in the set of random variables for ray-tracing => PD is not integrated in RESMAT
    ! OBJ must know where to find corresponding random variable in the XRND array:
    ! it will be obtained as XRND(IXRND(OBJ%FRAME%IRNDP))
      if (OBJ%TYP.NE.ctyp_bent) then
        if (OBJ%MODEL.EQ.0) then
          IXRND(NTRM)=NTRM
          OBJ%FRAME%IRNDP=NTRM
        endif
      endif

! update CMAT by penetration matrix along k
      call GetTransCryst(OBJ%FRAME,OBJ%G,Cin,1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cin,5,5,(/1,2,3,7,INDX(8)/))

! transmission matrix in local coordinates
      T=0.D0
      T(1,1)=1.D0/(s_step*OBJ%FRAME%SIZE(1))**2
      T(2,2)=1.D0/(s_step*OBJ%FRAME%SIZE(2))**2
      T(3,3)=1.D0/(s_step*OBJ%FRAME%SIZE(3))**2

  ! mosaic crystal: add random horizontal mosaic angle, ETA
  ! NOTE: the case of perfect undeformed crystal is not yet solved
  ! HMOS, VMOS are dispersions, not FWHM!
      if (OBJ%TYP.NE.ctyp_bent) then
        NTRM=NTRM+1
        NX=NX+1
        INDX(NX)=NTRM ! index of mosaic angle variable
        MATLEG=trim(MATLEG)//'|hmos'
        T(NX,NX)=1.D0/OBJ%HMOS**2
      endif
  ! add vertical mosaicity if any
      if (OBJ%VMOS.GT.sec) then
        NTRM=NTRM+1
        NX=NX+1
        INDX(NX)=NTRM ! index of mosaic angle variable
        MATLEG=trim(MATLEG)//'|vmos'
        T(NX,NX)=1.D0/OBJ%VMOS**2
      ! If MODEL=0, include vertical mosaicity in the set of random variables (as PD)
        if (OBJ%MODEL.EQ.0) then
          IXRND(NTRM)=NTRM
          OBJ%IRNDV=NTRM
        endif
      endif

      if (DBG) call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT 3',MATLEG,MATLEG,CMAT(1,1),MD)

! Update transmission matrix
      call AddBTAB(TRM,T,10,NX,CMAT,MD,NTRM,INDX(1))

! calculate DG=grad(G)
      DG=0.D0
  ! horizontal focusing
      if ((OBJ%TYP.eq.ctyp_gradient).or.(OBJ%TYP.eq.ctyp_bent)) then
        do i=1,3,2
          do j=1,3,2
            DG(i,j)=OBJ%DG_DR(i,j)
          enddo
        enddo
      else if (OBJ%NH.GT.1) then
        DG(1,1)=-OBJ%G(3)*OBJ%RH
        DG(3,1)=OBJ%G(1)*OBJ%RH
      endif
  ! stack misalignment along z
      if ((OBJ%NB>1).and.(OBJ%RB.NE.0.D0)) then
        DG(1,3)=DG(1,3)+OBJ%G(3)*OBJ%RB
        DG(3,3)=DG(3,3)-OBJ%G(1)*OBJ%RB
      endif
  ! vertical focusing
      if (OBJ%NV.GT.1) then
        DG(2,2)=DG(2,2)-OBJ%G(3)*OBJ%RV
      endif


! Define Bragg constraint vector
      B=0.D0
      B(1)=KG(1)*DG(1,1) + KG(3)*DG(3,1)
      B(2)=KG(2)*DG(2,2)
      B(3)=KG(1)*DG(1,3) + KG(3)*DG(3,3)
      B(4)=KABS*OBJ%G(1)
      B(5)=KABS*OBJ%G(2)
      B(6)=KABS*OBJ%G(3)
      B(9)= OBJ%G(3)*K0(1) - OBJ%G(1)*K0(3)
      B(10)=OBJ%G(2)*K0(3) - OBJ%G(3)*K0(2)
  ! create constraint matrix for this crystal
      do j=1,10
        DO i=j,10
          MB(i,j)=B(i)*B(j)/KSQR**2
          MB(j,i)=MB(i,j)
        enddo
      enddo


      if (DBG) call XML_MATRIX(IUDBG,3,3,'DG','x|y|z','x|y|z',DG(1,1),3)



! Update constraint matrix
      call AddBTAB(CON,MB,10,NX,CMAT,MD,NTRM,INDX(1))

! post-process
  ! change of beam direction due to diffraction:
  ! Cg fields correspond to the variables: x,y,z,kx,ky,kz,t,pd,hmos,vmos
      call UNIMAT(Cg,10,10)
      DO i=1,3
        do j=1,3
          Cg(i+3,j)=DG(i,j)/KABS
        enddo
      enddo
      Cg(4,9)=OBJ%G(3)/KABS
      Cg(5,10)=-OBJ%G(3)/KABS
      Cg(6,9)=-OBJ%G(1)/KABS+OBJ%G(2)/KABS
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cg,10,NX,INDX(1))
      if (DBG) call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT 4',MATLEG,MATLEG,CMAT(1,1),MD)

  ! shift by penetration depth outside
      call GetTransCryst(OBJ%FRAME,OBJ%G,Cout,-1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cout,5,5,(/1,2,3,7,INDX(8)/))
      if (DBG) call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT 5',MATLEG,MATLEG,CMAT(1,1),MD)

  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))

      if (DBG) then
      !  write(IUDBG,*) 'ADDCRYSTAL: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD CRYST','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,10,10,'Cg',MATLEG,MATLEG,Cg(1,1),10)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CON',MATLEG,MATLEG,CON(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif
      end subroutine RESMAT_ADDCRYSTAL

!-----------------------------------------------------
! Add a crystal component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDXTAL(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(XTAL),pointer :: OBJ
      TYPE(PXTAL) :: POBJ
      integer,parameter :: CD=9
      integer :: i,j,INDX(9),NX,SELNODE,iref
      real(KIND(1.D0)) :: T(9,9),C(7,7),Cin(5,5),Cout(5,5),Cg(9,9),B(9), MB(9,9), DG(3,3)
      real(KIND(1.D0)) :: K0(3),KG(3),KABS,KSQR,GABS,DBU,BDEV,AUX
      TYPE (PXREF) :: PXR
      TYPE (PRLNODE) :: PNODE
      call XTAL_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
  ! local indices correspond to (x,y,z,kx,ky,kz,t,PD)
  ! PD is the penetration depth, measured as time of flight in units of h_bar/m.
      INDX=(/1,2,3,4,5,6,7,0,0/)
      NX=7 ! maximum index actually used
  ! get root node
      call RLNODE_GET_ITEM(OBJ%ICR,1,PNODE)
  ! get selected node index
      SELNODE=RLNODE_GET(OBJ%ICR,OBJ%HKL)
  ! no reflection, or selected HKL=(0,0,0)
      ! write(*,*) 'RESMAT_ADDXTAL SELNODE=',SELNODE,' NP=',PNODE%N%NP
      if ((PNODE%N%NP.le.0).or.(SELNODE.eq.1)) then
        call RESMAT_ADDNONE(OBJ%FRAME)
        return
      endif
      iref=XREF_GET(OBJ%ICR,OBJ%HKL)
    ! in the case of forbidden reflections, use the 1st on the list
      if (iref.le.0) iref=PNODE%N%PATH(1)
      call XREF_GET_ITEM(OBJ%ICR,iref,PXR)
    ! get Darwin box unit
      DBU=2.D0*GetCrystDBU(OBJ,PXR)
    ! add misfits for all allowed primary reflections
      BDEV=0.D0
      do i=1,PNODE%N%NP
        call XREF_GET_ITEM(OBJ%ICR,PNODE%N%PATH(i),PXR)
        AUX=abs(GetCrystBragDev(OBJ,PXR))
        BDEV=max(BDEV,AUX)
        write(*,*) 'ADDXTAL: ',PXR%X%HKL,' BDEV=',NINT(AUX/DBU)
      enddo
      DBU=max(BDEV,DBU)

      KABS=2.D0*PI/OBJ%LAMBDA
  ! update CMAT by pre-process matrix
      call GetPreTrans(OBJ%FRAME,C)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))
  ! get nominal incident k-vector K0 in local coordinates and KK=|K0|^2
  ! get k-vector in local coordinates
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,KABS/),K0)
  ! |G|
      GABS=SQRT(PXR%X%G(1)**2+PXR%X%G(2)**2+PXR%X%G(3)**2)
  ! KG=K0+G
      do i=1,3
        KG(i)=K0(i)+PXR%X%G(i)
      enddo
      KSQR=KABS**2
      ! KNOMINAL=KABS
! Add penetration depth, PD
! For a bent perfect crystal, it is actually not a variable
! but it will be defined through the Bragg constraint later
      NTRM=NTRM+1
      NX=8
      INDX(NX)=NTRM ! index of PD variable
      MATLEG=trim(MATLEG)//'|pd'

! Add random deviation from Bragg law (Darwin box width)
      NTRM=NTRM+1
      NX=9
      INDX(NX)=NTRM ! index of PD variable
      MATLEG=trim(MATLEG)//'|dbw'

! update CMAT by penetration matrix along k
      call GetTransCryst(OBJ%FRAME,PXR%X%G,Cin,1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cin,5,5,(/1,2,3,7,INDX(8)/))
! transmission matrix in local coordinates

      T=0.D0
      T(1,1)=1.D0/(s_step*OBJ%A3D%SA(1))**2
      T(2,2)=1.D0/(s_step*OBJ%A3D%SA(2))**2
      T(3,3)=1.D0/(s_step*OBJ%A3D%SA(3))**2
! Darwin box width
      T(9,9)=1.D0/(s_step*DBU)**2

! Update transmission matrix
      call AddBTAB(TRM,T,9,NX,CMAT,MD,NTRM,INDX(1))

! get G-gradient
      DG=0.D0
      do i=1,3
        do j=1,3
            DG(i,j)=PXR%X%DG_DR(i,j)
        enddo
      enddo
    ! vertical focusing by segments, if not included in the crystal gradient
      if (OBJ%VFOC%NSEG.GT.1) then
        DG(2,2)=DG(2,2)-PXR%X%G(3)*OBJ%VFOC%RHO
        DG(3,2)=DG(3,2)+PXR%X%G(2)*OBJ%VFOC%RHO
      endif

! Define Bragg constraint vector
      B=0.D0
      B(1)=KG(1)*DG(1,1) + KG(2)*DG(2,1) + KG(3)*DG(3,1)
      B(2)=KG(1)*DG(1,2) + KG(2)*DG(2,2) + KG(3)*DG(3,2)
      B(3)=KG(1)*DG(1,3) + KG(2)*DG(2,3) + KG(3)*DG(3,3)
      B(4)=KABS*PXR%X%G(1)
      B(5)=KABS*PXR%X%G(2)
      B(6)=KABS*PXR%X%G(3)
! add misfit by Darwin box to Bragg condition
! this should help to avoid problems in symmetric transmission geometry
      B(9)=1.D0
      do i=1,6
        if (abs(B(i)).lt.1.D-12) B(i)=0.D0
      enddo
  ! create constraint matrix for this crystal
      MB=0.D0
      do j=1,9
        DO i=j,9
          MB(i,j)=B(i)*B(j)/KSQR**2
          MB(j,i)=MB(i,j)
        enddo
      enddo
! Update constraint matrix
      call AddBTAB(CON,MB,9,NX,CMAT,MD,NTRM,INDX(1))
! post-process
  ! change of beam direction due to diffraction:
  ! Cg fields correspond to the variables: x,y,z,kx,ky,kz,t,pd
      call UNIMAT(Cg,9,9)
      DO i=1,3
        do j=1,3
          Cg(i+3,j)=DG(i,j)/KABS
        enddo
      enddo
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cg,9,NX,INDX(1))
  ! shift by penetration depth outside
      call GetTransCryst(OBJ%FRAME,PXR%X%G,Cout,-1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cout,5,5,(/1,2,3,7,INDX(8)/))
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))

      if (DBG) then
      !  write(IUDBG,*) 'ADDXTAL: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD XTAL','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,9,9,'Cg',MATLEG,MATLEG,Cg(1,1),9)
        call XML_ARRAY(IUDBG,9,'Bragg',MATLEG,B(1))
        call XML_MATRIX(IUDBG,9,9,'MB',MATLEG,MATLEG,MB(1,1),9)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CON',MATLEG,MATLEG,CON(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif
      end subroutine RESMAT_ADDXTAL



!-----------------------------------------------------
! Add a crystal component to the beamline
!-----------------------------------------------------
      subroutine RESMAT_ADDXTAL_OLD(COMOBJ)
      TYPE(TCOMOBJ),intent(in) :: COMOBJ
      TYPE(XTAL),pointer :: OBJ
      TYPE(PXTAL) :: POBJ
      integer :: i,j,INDX(8),NX,SELNODE,iref
      real(KIND(1.D0)) :: T(8,8),C(7,7),Cin(5,5),Cout(5,5),Cg(8,8),B(8), MB(8,8), DG(3,3)
      real(KIND(1.D0)) :: K0(3),KG(3),KABS,KSQR,GABS
      TYPE (PXREF) :: PXR
      TYPE (PRLNODE) :: PNODE
      call XTAL_GET(COMOBJ%INST,POBJ)
      OBJ => POBJ%X
  ! local indices correspond to (x,y,z,kx,ky,kz,t,PD)
  ! PD is the penetration depth, measured as time of flight in units of h_bar/m.
      INDX=(/1,2,3,4,5,6,7,0/)
      NX=7 ! maximum index actually used
  ! get root node
      call RLNODE_GET_ITEM(OBJ%ICR,1,PNODE)
  ! get selected node index
      SELNODE=RLNODE_GET(OBJ%ICR,OBJ%HKL)
  ! no reflection, or selected HKL=(0,0,0)
      if ((PNODE%N%NP.le.0).or.(SELNODE.eq.1)) then
        call RESMAT_ADDNONE(OBJ%FRAME)
        return
      endif
      iref=XREF_GET(OBJ%ICR,OBJ%HKL)
    ! in the case of forbidden reflections, use the 1st on the list
      if (iref.le.0) iref=PNODE%N%PATH(1)
      call XREF_GET_ITEM(OBJ%ICR,iref,PXR)
      KABS=2.D0*PI/OBJ%LAMBDA
  ! update CMAT by pre-process matrix
      call GetPreTrans(OBJ%FRAME,C)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))
  ! get nominal incident k-vector K0 in local coordinates and KK=|K0|^2
  ! get k-vector in local coordinates
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,KABS/),K0)
  ! |G|
      GABS=SQRT(PXR%X%G(1)**2+PXR%X%G(2)**2+PXR%X%G(3)**2)
  ! KG=K0+G
      do i=1,3
        KG(i)=K0(i)+PXR%X%G(i)
      enddo
      KSQR=KABS**2
      ! KNOMINAL=KABS
! Add penetration depth, PD
! For a bent perfect crystal, it is actually not a variable
! but it will be defined through the Bragg constraint later
      NTRM=NTRM+1
      NX=8
      INDX(NX)=NTRM ! index of PD variable
      MATLEG=trim(MATLEG)//'|pd'
! update CMAT by penetration matrix along k
      call GetTransCryst(OBJ%FRAME,PXR%X%G,Cin,1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cin,5,5,(/1,2,3,7,INDX(8)/))
! transmission matrix in local coordinates

      T=0.D0
      T(1,1)=1.D0/(s_step*OBJ%A3D%SA(1))**2
      T(2,2)=1.D0/(s_step*OBJ%A3D%SA(2))**2
      T(3,3)=1.D0/(s_step*OBJ%A3D%SA(3))**2
! Update transmission matrix
      call AddBTAB(TRM,T,8,NX,CMAT,MD,NTRM,INDX(1))

! get G-gradient
      do i=1,3,2
        do j=1,3,2
            DG(i,j)=PXR%X%DG_DR(i,j)
        enddo
      enddo
    ! vertical focusing
      if (OBJ%VFOC%NSEG.GT.1) then
        DG(2,2)=-PXR%X%G(3)*OBJ%VFOC%RHO
      endif

! Define Bragg constraint vector
      B=0.D0
      B(1)=KG(1)*DG(1,1) + KG(3)*DG(3,1)
      B(3)=KG(1)*DG(1,3) + KG(3)*DG(3,3)
      B(4)=KABS*PXR%X%G(1)
      B(6)=KABS*PXR%X%G(3)
      do i=1,6
        if (abs(B(i)).lt.1.D-12) B(i)=0.D0
      enddo
  ! create constraint matrix for this crystal
      MB=0.D0
      do j=1,6
        DO i=j,6
          MB(i,j)=B(i)*B(j)/KSQR**2
          MB(j,i)=MB(i,j)
        enddo
      enddo
! Update constraint matrix
      call AddBTAB(CON,MB,8,NX,CMAT,MD,NTRM,INDX(1))
! post-process
  ! change of beam direction due to diffraction:
  ! Cg fields correspond to the variables: x,y,z,kx,ky,kz,t,pd
      call UNIMAT(Cg,8,8)
      DO i=1,3
        do j=1,3
          Cg(i+3,j)=DG(i,j)/KABS
        enddo
      enddo
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cg,8,NX,INDX(1))
  ! shift by penetration depth outside
      call GetTransCryst(OBJ%FRAME,PXR%X%G,Cout,-1)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,Cout,5,5,(/1,2,3,7,INDX(8)/))
  ! transformation from local to exit coord. system
      call GetPostTrans(OBJ%FRAME,C)
      call MltpBySubmatrix(0,CMAT,MD,NTRM,C,7,7,INDX(1))

      if (DBG) then
      !  write(IUDBG,*) 'ADDXTAL: ','class=',OBJ%FRAME%CLASS,' ',trim(OBJ%FRAME%NAME)
        call XML_VALUE(IUDBG,'ADD XTAL','string',' ',NTRM,0.D0,trim(OBJ%FRAME%NAME))
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,8,8,'Cg',MATLEG,MATLEG,Cg(1,1),8)
        call XML_ARRAY(IUDBG,6,'Bragg',MATLEG,B(1))
        call XML_MATRIX(IUDBG,8,8,'MB',MATLEG,MATLEG,MB(1,1),8)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CON',MATLEG,MATLEG,CON(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM',MATLEG,MATLEG,TRM(1,1),MD)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT',MATLEG,MATLEG,CMAT(1,1),MD)
        ! READ(*,*)
      endif
      end subroutine RESMAT_ADDXTAL_OLD

!C**************************************************************************
!C     Auxilliary subroutines:
!C**************************************************************************

!------------------------------------------------------------------------
!C Integrate over the IVAR-th variable, reducing matrices TRM,CMAT and CON
!----------------------------------------------------------------------
      subroutine IntegrateVariable(IVAR)
      integer,intent(in) :: IVAR
      integer :: i,j,k,IS,IL,LL
      real(KIND(1.D0)),parameter :: EPS=1.D-12
      real(KIND(1.D0)) :: AUX(MD,MD)
      logical :: is_integrable,isConstrained
      character(LEN=2048) TMPLEG
      character(LEN=16) :: VN
      VN=' '
      !write(*,*) 'integrate variable ',IVAR,' from ',NTRM
      if (DBG) then
        write(IUDBG,*) 'integrate variable ',IVAR,' from ',NTRM
        write(IUDBG,*) trim(MATLEG)
        ! READ(*,*)
      endif
      if ((IVAR.GT.0).and.(IVAR.LE.NTRM)) then
  ! don't integrate if the variable is used by a component
        is_integrable=.TRUE.
        call FINDSTRPAR(MATLEG,'|',IVAR,IS,IL)
        if (IL.GT.0) VN=MATLEG(IS:IS+IL-1)
        do i=1,CRND
          is_integrable=(is_integrable.and.(IXRND(i).ne.IVAR))
        enddo
        if (.not.is_integrable) then
          if (DBG) then
            write(IUDBG,*) '   ',trim(VN),' is not integrable '
          endif
          return
        endif
        if (DBG) write(IUDBG,*) '   ',trim(VN),' is integrable '
    ! remove legend item from the list
        LL=len_trim(MATLEG)
        if (IL.GT.0) then
          TMPLEG=MATLEG
          if (LL.gt.IS+IL) then
            if (IS.eq.1) then
              MATLEG=TMPLEG(IS+IL+1:LL)
            else
              MATLEG=TMPLEG(1:IS-2)//TMPLEG(IS+IL:LL)
            endif
          else
            MATLEG=TMPLEG(1:IS-2)
          endif
        endif

  ! integrate the transmission matrix
        k=IVAR

        isConstrained=(ABS(CON(k,k)).gt.0.D0)
        if (isConstrained) then
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) then
                  isConstrained = (abs(CON(i,k)*CON(k,j)).lt.CON(k,k)*1.D10)
                endif
                if (.not.isConstrained) exit
              enddo
            endif
            if (.not.isConstrained) exit
          enddo
        endif
      ! there is a non-zero constraint k-element
        if (isConstrained) then
          if (DBG) write(IUDBG,*) '   ',trim(VN),' is constrained: ',CON(k,k)
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) AUX(i,j)=TRM(i,j)-(CON(k,j)*TRM(i,k)+CON(i,k)*TRM(k,j)-CON(i,j)*TRM(k,k))/CON(k,k)
              enddo
            endif
          enddo
        ! new constraint matrix after integration
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) CON(i,j)=CON(i,j)-CON(i,k)*CON(k,j)/CON(k,k)
              enddo
            endif
          enddo
      ! no constraints
        else if (ABS(TRM(k,k)).gt.EPS) then
          if (DBG) write(IUDBG,*) '   ',trim(VN),' unconstrained, diag=',TRM(k,k)
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) AUX(i,j)=TRM(i,j)-TRM(i,k)*TRM(k,j)/TRM(k,k)
              enddo
            endif
          enddo
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) CON(i,j)=CON(i,j)-(CON(k,j)*TRM(i,k)+CON(i,k)*TRM(k,j))/TRM(k,k)
              enddo
            endif
          enddo
        endif
  ! reduce dimension
        call DeleteIndex(AUX,MD,NTRM,k,k)
        call DeleteIndex(CON,MD,NTRM,k,k)
        call DeleteIndex(CMAT,MD,NTRM,k,k)
    ! update pointers to XRND array (used by components)
        do i=1,CRND
          if (IXRND(i).gt.k) IXRND(i)=IXRND(i)-1
        enddo
        NTRM=NTRM-1
        do i=1,NTRM
          do j=1,NTRM
            TRM(i,j)=AUX(i,j)
          enddo
        enddo
        if (DBG) then
          write(IUDBG,*) '   ','new dimension: ',NTRM
        !  READ(*,*)
        endif
      endif
      !write(*,*) '   ','new dimension: ',NTRM
      end subroutine IntegrateVariable

!------------------------------------------------------------
      subroutine AddSubmatrix(A,DIMA,NA,B,DIMB,NB,INDX)
! Add a submatrix B to the matrix A
! DIMA,DIMB are leading dimensions of the arrays A,B
! NA,NB are the ranks of the matrices A,B
! NA >= NB
! INDX(NA) contains corresponding indices, INDX(1) is the index of B-th 1st row/column in A
!------------------------------------------------------------
      integer, intent(in) :: DIMA,DIMB,NA,NB,INDX(NB)
      real(KIND(1.D0)), intent(in) :: B(DIMB,NB)
      real(KIND(1.D0)) :: A(DIMA,NA)
      integer :: i,j,xi,xj
      do i=1,NB
        xi=INDX(i)
        do j=1,NB
          xj=INDX(j)
          A(xi,xj)=A(xi,xj)+B(i,j)
        enddo
      enddo
      end subroutine AddSubmatrix


!------------------------------------------------------------
      subroutine AddBTAB(C,A,DIMA,NA,B,DIMB,NB,INDX)
! Add a submatrix A to the matrix C, using conversion matrix B
! That is, C=C+B^T.A.B
! DIMA and DIMB are leading dimensions of the arrays A and B
! NA and NB are the ranks of the matrices A and B
! C has the same dimensions as B
! NB >= NA
! INDX(NA) contains corresponding indices
!------------------------------------------------------------
      integer, intent(in) :: DIMA,DIMB,NA,NB,INDX(NA)
      real(KIND(1.D0)), intent(in) :: B(DIMB,NB),A(DIMA,NA)
      real(KIND(1.D0)) :: C(DIMB,NB)
      integer :: i,j,k,m,xk,xm
c      write(*,*) 'AddBTAB: ',INDX(1:NA)
      do i=1,NB
        do j=1,NB
          do k=1,NA
            xk=INDX(k)
            do m=1,NA
              xm=INDX(m)
              if ((xm<1).or.(xk<1)) then
                 write(*,*) 'AddBTAB ',xk, xm, NA, NB, ' [', INDX,'] '
              else
                  C(i,j)=C(i,j)+B(xk,i)*A(k,m)*B(xm,j)
              endif
            enddo
          enddo
        enddo
      enddo
      end subroutine AddBTAB

!------------------------------------------------------------
      subroutine MltpBySubmatrix(ISIDE,A,DIMA,NA,B,DIMB,NB,INDX)
! Multiply matrix A by a submatrix B
! from left (A=B.A, ISIDE=0) or from right (A=A.B, ISIDE=1)
! DIMA,DIMB are leading dimensions of the arrays A,B
! NA,NB are sizes of the square matrices
! NA >= NB
! INDX(NA) contains corresponding indices, INDX(1) is the index of B-th 1st row/column in A
!------------------------------------------------------------
      integer, intent(in) :: ISIDE,DIMA,DIMB,NA,NB,INDX(NB)
      real(KIND(1.D0)), intent(in) :: B(DIMB,NB)
      real(KIND(1.D0)) :: A(DIMA,NA),C(NB)
      integer :: i,j,k,xk,xj
      select case(ISIDE)
      case(0)
  ! from left
        do i=1,NA ! scan columns
          do j=1,NB
            C(j)=0.D0
            do k=1,NB
              xk=INDX(k)
              C(j)=C(j)+B(j,k)*A(xk,i)
            enddo
          enddo
          do j=1,NB
            xj=INDX(j)
            A(xj,i)=C(j)
          enddo
        enddo
      case(1)
  ! from right
        do i=1,NA ! scan rows
          do j=1,NB
            C(j)=0.D0
            do k=1,NB
              xk=INDX(k)
              C(j)=C(j)+A(i,xk)*B(k,j)
            enddo
          enddo
          do j=1,NB
            xj=INDX(j)
            A(i,xj)=C(j)
          enddo
        enddo
      end select
      end subroutine MltpBySubmatrix

!-----------------------------------------------------
      subroutine DeleteIndex(A,DIMA,NA,IROW,ICOL)
! Delete one row and column (i.e. compact and lower matrix rank by 1)
! DIMA is the leading dimensions of A
! NA is the rank of the square matrix A
! IROW,ICOL are the indices of row and column to be deleted
!-----------------------------------------------------
      integer, intent(in) :: DIMA,NA,IROW,ICOL
      real(KIND(1.D0)) :: A(DIMA,NA)
      integer :: i,j
  ! remove row
      do i=1,NA
        do j=IROW,NA-1
          A(j,i)=A(j+1,i)
        enddo
        A(NA,i)=0.D0
      enddo
  ! remove column
      do i=1,NA-1
        do j=ICOL,NA-1
          A(i,j)=A(i,j+1)
        enddo
        A(i,NA)=0.D0
      enddo
      end subroutine DeleteIndex


!-----------------------------------------------------
      subroutine GetTransCryst(OBJ,G,C,dir)
! Get transformation matrix C which moves neutron
! to or from the point of reflection in the crystal by pd
! pd is penetration depth measured along the -G vector
! dir>0 ... along k
! dir<0 ... along k+G
! C should be rank=5 square matrix
! indices in C correspond to (x,y,z,t,pd)
!-----------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)), intent(in)  :: G(3)
      real(KIND(1.D0)), intent(out)  :: C(5,5)
      integer, intent(in) :: dir
      integer :: i
      real(KIND(1.D0)) :: K0(3),KK,sinB,vG,GABS
      integer :: ISG
      ISG=sign(1,dir)
      call UNIMAT(C,5,5)
  ! get ki-vector in local coordinates
      CALL M3XV3(OBJ%MLOC,OBJ%RLOC,(/0.D0,0.D0,1.D0/),K0)
    ! sine of Bragg angle is calculated from the current crystal orientation
      GABS=SQRT(G(1)**2+G(2)**2+G(3)**2)
      sinB=abs(v3xv3(K0,G)/GABS)
    ! |k| value
      KK=GABS/2.D0/sinB
    ! velocity along -G
      vG=ISG*GABS/2.D0
  !    write(IUDBG,*) 'GetTransCryst: vG=',vG,' sinB=',sinB,' K0=',K0
    ! add path along k or k+G vector
      do i=1,3
        K0(i)=KK*K0(i)
        if (dir.lt.0) K0(i)=K0(i)+G(i)
        C(i,5)=K0(i)/vG
      endDO
      C(4,5)=1.D0/vG
      if (DBG) then
        write(IUDBG,*) 'GetTransCryst: ',trim(OBJ%ID)
        call XML_MATRIX(IUDBG,5,5,'Cz','x|y|z|t|pd','x|y|z|t|pd',C(1,1),5)
      endif
      end subroutine GetTransCryst


!-----------------------------------------------------
      real(KIND(1.D0)) function GetCrystDBU(OBJ,PXR)
! Get darwin box unit for given reflection
! It is defined so that if the deviation from Bragg angle = fwhm of the Darwin box (y=1), then
!  DBU = (K+G)^2-K^2
! assume crystal surface along (0,0,1)
! DBU is in [A^-2]
!-----------------------------------------------------
      TYPE(XTAL) :: OBJ
      TYPE (PXREF) :: PXR
      real(KIND(1.D0)) :: K0(3),KG(3),FHKL_VC,DBU
      ! get ki-vector in local coordinates
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,2.D0*PI/OBJ%LAMBDA/),K0)
      KG=K0+PXR%X%G
      if (abs(K0(3)).gt.1e-6) then
      ! DB = 8.pi.Fhkl/vc*|cos_g/cos|
      ! gt Fhkl/vc from QML (QML is in [A^-1 cm^-1])
        FHKL_VC=1.D-4*sqrt(PXR%X%QML/4.D0/PI)/PXR%X%dhkl
      ! DB = 8.pi.Fhkl/vc*|cos_g/cos|
        DBU=8.D0*PI*FHKL_VC*sqrt(abs(KG(3)/K0(3)))
      else
        DBU=1.e-6
      endif
      GetCrystDBU=DBU
      end function GetCrystDBU


!-----------------------------------------------------
      real(KIND(1.D0)) function GetCrystBragDev(OBJ,PXR)
! Get deviation from Bragg condition as (K+G)^2-K^2
!-----------------------------------------------------
      TYPE(XTAL) :: OBJ
      TYPE (PXREF) :: PXR
      real(KIND(1.D0)) :: K0(3)
      ! get ki-vector in local coordinates
      CALL M3XV3(OBJ%FRAME%MLOC,OBJ%FRAME%RLOC,(/0.D0,0.D0,2.D0*PI/OBJ%LAMBDA/),K0)
      GetCrystBragDev=V3XV3(PXR%X%G,PXR%X%G)+2.D0*V3XV3(K0,PXR%X%G)
      end function GetCrystBragDev


!-----------------------------------------------------
      subroutine GetTransZ0(OBJ,C)
! Get transformation matrix C which moves neutron to z=0 plane
! C should be rank=4 square matrix
! indices in C correspond to (x,y,z,t)
!-----------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)), intent(out)  :: C(4,4)
      integer :: i,iref
      real(KIND(1.D0)) :: K0(3)
      call UNIMAT(C,4,4)
  ! get beam direction in local coordinates
      CALL M3XV3(OBJ%MLOC,OBJ%RLOC,(/0.D0,0.D0,KNOMINAL/),K0)
  ! move to z=0 (or x=1 or y=1, if k(3)=0), matrix equivalent of SLIT_MIDDLE
      if (abs(K0(3)/KNOMINAL).gt.1.D-6) then
        iref=3
      else if (abs(K0(1)/KNOMINAL).gt.1.D-6) then
        iref=1
      else
        iref=2
      endif
      do i=1,3
        C(i,iref)=C(i,iref)-K0(i)/K0(iref)
      enddo
      C(4,iref)=C(4,iref)-1.D0/K0(iref)
      if (DBG) then
        write(IUDBG,*) 'GetTransZ0: ',trim(OBJ%ID)
        call XML_MATRIX(IUDBG,4,4,'Cz','x|y|z|t','x|y|z|t',C(1,1),4)
      endif
      end subroutine GetTransZ0

!-----------------------------------------------------
      subroutine GetPreTrans(OBJ,C)
! Get pre-transformation matrix C for an object
! C should be rank=7 square matrix
! indices in C correspond to (x,y,z,kx,ky,kz,t)
!-----------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)) :: C(7,7)

  ! transformation from the exit system of preceding component
  ! to the incident system of OBJ
      call UNIMAT(C,7,7)
      C(1,4)=OBJ%DIST
      C(2,5)=OBJ%DIST
      C(7,6)=-OBJ%DIST/KNOMINAL
  ! transformation from incident to local coordinates, C=OBJ%RLOC.C
      call MltpBySubmatrix(0,C,7,7,OBJ%RLOC,3,3,(/1,2,3/))
      call MltpBySubmatrix(0,C,7,7,OBJ%RLOC,3,3,(/4,5,6/))
      if (DBG) then
        write(IUDBG,*) 'pre-trans: ',KNOMINAL
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,7,7,'C',MATLEG,MATLEG,C(1,1),7)
      endif
      end subroutine GetPreTrans


!-----------------------------------------------------
      subroutine GetLoc2Entry(OBJ,C)
! Get transformation matrix C for converstion
! from local to entry coordinates
! C should be rank=7 square matrix
! indices in C correspond to (x,y,z,kx,ky,kz)
!-----------------------------------------------------
      TYPE(TFRAME) :: OBJ
      real(KIND(1.D0)) :: C(6,6)
      real(KIND(1.D0)) :: U(3,3),R(3,3)
  ! transformation from the exit system of preceding component
  ! to the incident system of OBJ
      call UNIMAT(C,6,6)
    ! get transposed RLOC
      call UNIMAT(U,3,3)
      Call M3XM3(-1,OBJ%RLOC,U,R)
  ! transformation from incident to local coordinates, C=OBJ%RLOC.C
      call MltpBySubmatrix(0,C,6,6,R,3,3,(/1,2,3/))
      call MltpBySubmatrix(0,C,6,6,R,3,3,(/4,5,6/))
      end subroutine GetLoc2Entry


!-----------------------------------------------------
! Get post-transformation matrix C for an object
! C should be rank=7 square matrix
! indices in C correspond to (x,y,z,kx,ky,kz,t)
!-----------------------------------------------------
      subroutine GetPostTrans(OBJ,C)
      TYPE(TFRAME) :: OBJ
      integer :: i,j,k
      real(KIND(1.D0)) :: C(7,7),R(3,3)
  ! transformation from local to exit coordinates
  ! R=REXI*RLOC^T
      DO i=1,3
        DO j=1,3
          R(i,j)=0.D0
          DO k=1,3
            R(i,j)=R(i,j)+OBJ%REXI(i,k)*OBJ%RLOC(j,k)
          enddo
        enddo
      enddo
  ! transformation from local to exit coord. system
      call UNIMAT(C,7,7)
      call MltpBySubmatrix(0,C,7,7,R,3,3,(/1,2,3/))
      call MltpBySubmatrix(0,C,7,7,R,3,3,(/4,5,6/))
      if (DBG) then
        write(IUDBG,*) 'post-trans: '
        call XML_VALUE(IUDBG,'dimension','int',' ',NTRM,0.D0,' ')
        call XML_MATRIX(IUDBG,7,7,'C',MATLEG,MATLEG,C(1,1),7)
      endif
      end subroutine GetPostTrans


      end module RESMAT
