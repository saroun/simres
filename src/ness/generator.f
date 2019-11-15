!//////////////////////////////////////////////////////////////////////
!////  $Id: generator.f,v 1.61 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.61 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module that implements a generator of neutrons (events) for ray-tracing
!////
!//////////////////////////////////////////////////////////////////////
      MODULE GENERATOR
      use FILETOOLS
      use SIMATH
      USE TRACINGDATA
      use TRACINGOPT
      use REPORTSOPT
      use NSTORE, only: NSTORE_MAXN,NSTORE_GETEVENT
      use VCTABLE_INSTRUMENTS
      use IO
      use RNDGEN
      use mtmod
      implicit none

C GENERATOR
C generator of NEUTRON events
C uses random numbers from RNDX() array (see modules RESMAT and TRACING)
! OLDCOUNT,OLDSCOUNT register the counter values when the safety pool is switched off
      TYPE  TGENERATOR; SEQUENCE
        !  REAL*8 WIDTH,HEIGHT ! WIDTH,HEIGHT
        !  REAL*8 DIVH,DIVV  ! divergence in [rad]
          REAL(KIND(1.D0)) :: K0     ! mean K=2*pi/lambda
          REAL(KIND(1.D0)) :: KMIN   ! minimum K
          REAL(KIND(1.D0)) :: KMAX   ! maxmum K
          ! K distribution width
          REAL(KIND(1.D0)) :: SCOUNT ! sum of generated weights
          REAL(KIND(1.D0)) :: OLDSCOUNT ! old sum of generated weights
          REAL(KIND(1.D0)) :: VOL    ! samplig volume
          ! integer*4 SHAPE   ! shape: (0) none (infinite), (1) rectangle, (2) elliptic
          ! INTEGER*4 HPROF,VPROF ! divergence profile: none(0), box(1), triangle (2), gaussian (3)
          INTEGER :: KPROF ! K profile: none(0), box(1), triangle (2), gaussian (3)
          INTEGER :: COUNT ! counter of generated events
          INTEGER :: OLDCOUNT ! old counter of generated events
          INTEGER :: INDEX ! reference to STACK1 event (for separated tracing through secondary beamline)
          integer :: IMAX  ! maximum events available in STACK1
      END TYPE TGENERATOR

! global instance of event generator
      TYPE(TGENERATOR) :: EVENTGEN


! covariance matrix of random variables
      REAL(KIND(1.D0)),private :: COV(CRND,CRND),MVAL(CRND),MAXV(CRND),SCOV  ! ,MINV(CRND),CTRV(CRND)
      INTEGER,private :: NCOV
      integer,private :: NPROBE=10 ! ratio of events/probes
      REAL(KIND(1.D0)),private :: RANG=3.0D0 ! range of probes
      REAL(KIND(1.D0)),private :: LIMITS=1.0D0 ! rnd var. limits

      REAL(KIND(1.D0)) :: GEN_DLAM=0.01D0 ! default wavelength range
      REAL(KIND(1.D0)), private :: GEN_DLL=0.D0 ! GEN_DLL = SRC_PERIOD/TOF_SAMPLE

  ! debug options
      logical, private :: DBG=.false.     ! flag
      integer, private :: IUDBG=6   ! Log-file unit

      private SAMPLING_VOLUME

      contains


C--------------------------------------------------------
      SUBROUTINE SET_LIMITS(PRANGE)
C Set probe range
C--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: PRANGE
      LIMITS=PRANGE
      write(*,*) 'Rnd var. limits = ',PRANGE
      end SUBROUTINE SET_LIMITS

C--------------------------------------------------------
      SUBROUTINE SETPROBE_FREQ(PN)
C Set probe events ratio
C--------------------------------------------------------
      integer,intent(in) :: PN ! ratio of events/probes
      NPROBE=PN
      write(*,*) 'Probe frequency = ',PN
      end SUBROUTINE SETPROBE_FREQ

C--------------------------------------------------------
      SUBROUTINE SETPROBE_RANGE(PRANGE)
C Set probe range
C--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: PRANGE ! range of probes
      RANG=PRANGE
      write(*,*) 'Probe range = ',PRANGE
      end SUBROUTINE SETPROBE_RANGE

C--------------------------------------------------------
      SUBROUTINE GENERATOR_SETTMAT(EVECT,EVALS,NDIM,N)
C Set conversion matrix and sampling limits
C--------------------------------------------------------
      integer,intent(in) :: NDIM,N
      REAL(KIND(1.D0)),intent(in) :: EVECT(NDIM,NDIM),EVALS(NDIM)
      integer :: i,j
!1     format(a,7(1x,G10.4))
! get transformation matrix and limits for events from RESMAT
      RNDLIST%DIM=N
      DO I=1,RNDLIST%DIM
      IF (RNDLIST%ACTIVE(I).GT.0) THEN
        RNDLIST%LIMITS(I)=LIMITS*RNDIST_RANGE*SQRT(ABS(EVALS(I)))
        DO J=1,RNDLIST%DIM
          IF (RNDLIST%ACTIVE(J).GT.0) THEN
             TMAT(I,J)=EVECT(I,J)
          ENDIF
        ENDDO
      ENDIF
      ENDDO
      RNDIST_VOL=SAMPLING_VOLUME()
      !write(*,1) 'GENERATOR_SETTMAT vol=',RNDIST_VOL
      !write(*,1) 'limits: ',RNDLIST%LIMITS(1:7)

      end SUBROUTINE GENERATOR_SETTMAT

!--------------------------------------------------------
      SUBROUTINE GENERATOR_SETMEAN(XY,KXY,KZ)
! Set sampling svolume ceter
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: XY(2),KXY(2),KZ
      TMEAN(1:2)=XY(1:2)
      TMEAN(3:4)=KXY(1:2)
      TMEAN(5)=KZ
      end SUBROUTINE GENERATOR_SETMEAN

C-------------------------------------------------
      Subroutine GENERATOR_INIT
C initialize generator parameters
C-------------------------------------------------
      integer :: I
      type(PDATASET) :: DTS
      !REAL(KIND(1.D0)) :: x
      !integer :: j
1     format(a,': ',6(1x,G12.5))
! default rand. var. pool settings
      DO I=1,CRND
        RNDLIST%MEAN(I)=0.D0
        RNDLIST%ACTIVE(I)=1
        RNDLIST%LIMITS(I)=1.D0
      enddo
C// initialize transformation matrix to a diagonal
      call UNIMAT(TMAT,CRND,CRND)
      TMEAN=0.D0
      call GENERATOR_CLEAR
C// set generator distribution if different from default (uniform)
      call GENERATOR_SETDISTR(TROPT%SAMPLING)
C// set IMAX=number of events available in STACK1 of the module NSTORE
      select case(TROPT%IMODE)
      case(tr_secondary)
        DTS=DSET_DATA(DSET_ISEL)
        EVENTGEN%IMAX=NSTORE_MAXN(DTS%P%REG1)
      case default
        EVENTGEN%IMAX=0
      end select
      if (DBG) then
          IUDBG=OPENFILEUNIT('GENERATOR.log',.false.)
          DBG=(IUDBG.gt.0)
      endif
      EVENTGEN%VOL=RNDIST_VOL
      EVENTGEN%OLDCOUNT=0
      EVENTGEN%OLDSCOUNT=0
      isProbe=.false.
      GEN_DLL = SRC_PERIOD/TOF_SAMPLE
      !if (PULSED_SRC) then
        !write(*,1) 'Wavelengh band width', GEN_DLL
      !endif
      !do j=-1,2
      !  x = 1.D0 + j*GEN_DLL
      !  if (x>1.D-3) then
      !    write(*,1) 'frame ', j,TWOPI/EVENTGEN%K0*x
      !  endif
      !enddo
      end Subroutine GENERATOR_INIT

C-------------------------------------------------
      Subroutine GENERATOR_SETK(K0,KMIN,KMAX)
C sets generator K0 and range
C-------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: K0,KMIN,KMAX
      if ((KMAX-KMIN).GT.0.D0) then
        EVENTGEN%KPROF=1
        EVENTGEN%KMIN=KMIN
        EVENTGEN%KMAX=KMAX
        EVENTGEN%K0=K0
      else
        call GENERATOR_SETDEFAULT(K0)
      endif
      end Subroutine GENERATOR_SETK

C-------------------------------------------------
      Subroutine GENERATOR_SETDEFAULT(K0)
C sets default generator properties (all dimensions free)
C-------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: K0
      EVENTGEN%KPROF=0
      EVENTGEN%KMIN=0.D0
      EVENTGEN%KMAX=0.D0
      EVENTGEN%K0=K0
      end Subroutine GENERATOR_SETDEFAULT

C-------------------------------------------------
      Subroutine GENERATOR_SETDISTR(idist)
C set distribution for rnd numbers
C-------------------------------------------------
      integer,intent(in) :: idist
      real(KIND(1.D0)),parameter :: CSQR2PI=2.506628274631000502415765284811D0
      real(KIND(1.D0)),parameter :: CPI=3.1415926535897932384626433832795D0
      integer :: I

      RNDIST_NORM=1.D0
      RNDIST=idist
      select case(RNDIST)
      case(rndist_uni)
        do i=1,CRND
          RNDLIST%POOL(I)=1.1D0
        ENDDO
        RNDIST_RANGE=2.D0*SQRT(6.D0)
      case(rndist_gauss)
        do i=1,CRND
          RNDLIST%POOL(I)=1.D0    ! no safety pool
        ENDDO
        RNDIST_NORM=CSQR2PI
        RNDIST_RANGE=1.D0
      case(rndist_cauchy)
        do i=1,CRND
          RNDLIST%POOL(I)=1.D0    ! no safety pool
        ENDDO
        RNDIST_NORM=1.D0
        RNDIST_RANGE=1.D0
      case default
        write(*,*) 'invalid random distribution type: ',RNDIST
        RNDIST=rndist_uni
        do i=1,CRND
          RNDLIST%POOL(I)=1.1D0
        ENDDO
        RNDIST_RANGE=1.D0*SQRT(6.D0)
      end select
c      write(*,*) 'GENERATOR_SETDISTR: ',idist,RNDIST
      RNDIST_VOL=SAMPLING_VOLUME()
      end Subroutine GENERATOR_SETDISTR

C-------------------------------------------------
      Subroutine GENERATOR_CLEAR
C-------------------------------------------------
      EVENTGEN%COUNT=0
      EVENTGEN%SCOUNT=0.D0
      EVENTGEN%INDEX=0
      EVENTGEN%VOL=RNDIST_VOL
      isProbe=.false.
      end Subroutine GENERATOR_CLEAR

C----------------------------------------------------------------
      SUBROUTINE GENERATOR_APPLYLIMITS(NEU,PP)
C adjust weight PP so that constraints given by the generator are applied
C--------------------------------------------------------------
      TYPE(NEUTRON),intent(in) :: NEU
      real(KIND(1.D0)),intent(out) :: PP
      if (EVENTGEN%KPROF.NE.0) then
        if ((NEU%K0.gt.EVENTGEN%KMAX).or.(NEU%K0.lt.EVENTGEN%KMIN)) PP=0.D0
      endif
      ! don't allow backwards neutrons
      if (NEU%K(3).le.0.D0) PP=0.D0
      end SUBROUTINE GENERATOR_APPLYLIMITS

!-------------------------------------------------
      SUBROUTINE GENERATOR_GO
!-------------------------------------------------
      real(KIND(1.D0)),parameter :: PLIM=1.D-30
      real(KIND(1.D0)) :: PP,cosa,cosb, DK,x
      integer :: i,j
      type(PDATASET) :: DTS

1     format(I8,12(1x,G11.5))

      PP=0.D0
      isProbe=((.not.isProbe).and.(mod(EVENTGEN%COUNT,NPROBE).eq.0))
! repeat until an event with non-zero weight is generated
      do while (PP.LE.PLIM)
! generate uniform random numbers, XNORM
  ! probe event is generated time to time
  ! it samples larger space to check sampling limits more efficiently
        if (isProbe) then
          call GENPROBE
          NEUT%LABEL=1
        else
          call GENEVENT
          NEUT%LABEL=0
        endif
        NEUT%STATE=0

        ! if (EVENTGEN%COUNT.lt.100) NEUT%LABEL=2

! transform to physical space, XNORM -> XRND
        XRND=0.D0
        CALL M16XV16(1,RNDLIST%DIM,TMAT,XNORM,XRND)
        do i=1,RNDLIST%DIM
          XRND(i)=XRND(i)+TMEAN(I)
        enddo


  ! debug values
  !    XRND(1:2)=(/0.192572D0,0.D0/)
  !    XRND(3:5)=(/4.010368D-3,0.D0,-5.40361D-3/)
  !    XRND(6:7)=(/4.867111D0,1.660987D-3/)

! generate random event NEUT
  ! from STACK1
        if (TROPT%IMODE.eq.tr_secondary) then
          EVENTGEN%INDEX=EVENTGEN%INDEX+1
          if (EVENTGEN%INDEX.gt.EVENTGEN%IMAX) EVENTGEN%INDEX=1
          DTS=DSET_DATA(DSET_ISEL)
          call NSTORE_GETEVENT(DTS%P%REG1,EVENTGEN%INDEX,NEUT)
          NEUT%REF(1:2)=(/DTS%P%REG1,EVENTGEN%INDEX/)
          PP=NEUT%P*PNORM
          if (.not.isProbe) then
            EVENTGEN%SCOUNT=EVENTGEN%SCOUNT+NEUT%P
            EVENTGEN%COUNT=EVENTGEN%COUNT+1
          endif

        !  write(*,"(a,I6,1x,3(G10.3,1x)))") 'GENERATOR_GO IDX=',EVENTGEN%COUNT,EVENTGEN%SCOUNT
        !  write(*,*) "PP=", PP,PNORM
        !  READ(*,*)
  ! from XRND
        else
          PP=PNORM
          NEUT%T0=0.D0
          NEUT%R(1)=XRND(1)
          NEUT%R(2)=XRND(2)
          NEUT%R(3)=0.0
        ! XRND(3)=kx/k_nom, XRND(4)=ky/k_nom, XRND(5)=(kz-k_nom)/k_nom
          if (TROPT%PSMODE==tr_kspace) then
            NEUT%K(1)=XRND(3)*EVENTGEN%K0
            NEUT%K(2)=XRND(4)*EVENTGEN%K0
            NEUT%K(3)=XRND(5)*EVENTGEN%K0+EVENTGEN%K0
            NEUT%K0=SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)
        ! XRND(3)=sina, XRND(4)=sinb, XRND(5)=-delta_lammbda/lambda_nom
          else if ((abs(XRND(3))<1.D0).and.(abs(XRND(4))<1.D0).and.(XRND(5)<1.D0)) then
            NEUT%K0=EVENTGEN%K0/(1.D0-XRND(5))
            cosb=sqrt(1.D0-XRND(4)**2)
            cosa=sqrt(1.D0-XRND(3)**2)
            NEUT%K(1)=XRND(3)*NEUT%K0*cosb
            NEUT%K(2)=XRND(4)*NEUT%K0
            NEUT%K(3)=NEUT%K0*cosa*cosb
            !NEUT%K(3)=SQRT(NEUT%K0**2-NEUT%K(1)**2-NEUT%K(2)**2)
            !PP=PP*NEUT%K0**4/TWOPI/cosa
            PP=PP*NEUT%K0**4/EVENTGEN%K0**4/cosa
          else
            PP=0.D0
          endif
          NEUT%S(1:3)=0.D0
          NEUT%S(1)=2*NINT(RAN1())-1
          if (TOF_MODE) then
            NEUT%T=XRND(6)
            if (PULSED_SRC.and.SRC_OVERLAP) then
            ! shift time and wavelength by a multiple of source period and bandwidth
            ! consider the interval (-1;+2) of source and target periods
              i=NINT(4.D0*RAN1()-1.5D0)
              j=NINT(4.D0*RAN1()-1.5D0) - i
              NEUT%T=NEUT%T+i*SRC_PERIOD
              x = 1.D0 + j*GEN_DLL
              if (x>1.D-3) then
                DK = EVENTGEN%K0/x - EVENTGEN%K0
                NEUT%K(3) = NEUT%K(3) + DK
                if (NEUT%K(3)>0.D0) then
                  NEUT%K0=SQRT(NEUT%K(1)**2+NEUT%K(2)**2+NEUT%K(3)**2)
                  PP=PP*16.D0
                else
                  PP=0.D0
                endif
              else
                DK = 1.D10
                PP=0.D0
              endif
            endif
            PP=PP/HOVM  ! integration is over time*hbar/m, but should be over time only
          else
            NEUT%T=0.D0
          endif
          NEUT%REF=0
          if (TROPT%DIR.eq.tr_upstream) NEUT%T=NEUT%T+TOF_SAMPLE
          NEUT%PHI=0.D0
          if (.not.isProbe) then
            EVENTGEN%SCOUNT=EVENTGEN%SCOUNT+1.D0
            EVENTGEN%COUNT=EVENTGEN%COUNT+1
          endif

      ! "reasonable" values
        !  if (ABS(XRND(3))+ABS(XRND(4))+ABS(XRND(5)).gt.3.D0) PP=0.D0
          if (PP.GT.PLIM) call GENERATOR_APPLYLIMITS(NEUT,PP)
        endif
        NEUT%CNT=trace_cnt+1


2     format(a,': ',10(1x,G12.5))
      !if (trace_tot<18) then
      !  write(*,2) 'GENERATOR_GO ',NINT(trace_tot), mtmod_ngen, trace_cnt, PP, XNORM(1:5),isProbe
      !endif

!        if (EVENTGEN%COUNT.eq.5000) then
!2         format("GENERATOR_GO gen=",I8,2(1x,a,G10.3),1x,a,1x,I5)
!          write(*,2) EVENTGEN%COUNT,' sgen=',EVENTGEN%SCOUNT,' vol=',RNDIST_VOL,' cnt=',trace_cnt
!          read(*,*)
!        endif
!        if (EVENTGEN%COUNT.eq.1) then
!          write(*,2) EVENTGEN%COUNT,' sgen=',EVENTGEN%SCOUNT,' vol=',RNDIST_VOL,' cnt=',trace_cnt
!          read(*,*)
!        endif

!      if (isProbe) then
!3     format(a,': ',8(G10.4,1x))
!        write(*,3) 'GENERATOR_GO probe: ',NEUT%R(1),NEUT%R(3),NEUT%K(1),NEUT%K(3),PP
!        READ(*,*)
!      endif

        if (DBG) then
          write(IUDBG,1) EVENTGEN%COUNT,NEUT%R(1:3),NEUT%K(1:2),NEUT%K(3)-EVENTGEN%K0,(XRND(i),i=6,RNDLIST%DIM)
          if (EVENTGEN%COUNT.ge.REPOPT%NRAYS) then
            close(IUDBG)
            IUDBG=6
            DBG=.false.
          endif
        endif

      enddo

! weight the event
      NEUT%P=PP
! initialize gravity vector = start in lab. coord. at the generator

c      write(*,1) 'GEN: ',NEUT%R,NEUT%K
c      read(*,*)

      end SUBROUTINE GENERATOR_GO

C------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION SAMPLING_VOLUME()
C get actual sampling volume
C------------------------------------------------
      REAL(KIND(1.D0))  :: vol
      integer :: i

  ! NOTE: K0*3 factor is needed when ki vector is included in random variables
  ! because source flux distribution has to be integrated over k-space
  ! for secondary beamline, the k-vector is not involved
      if (TROPT%IMODE.eq.tr_secondary) then
        vol=1.D0
      else
        vol=EVENTGEN%K0**3
      endif
      do i=1,RNDLIST%DIM
          vol=vol*RNDLIST%LIMITS(I)*RNDIST_NORM
      enddo
      SAMPLING_VOLUME=vol
      end FUNCTION SAMPLING_VOLUME


C-------------------------------------------------------------------
      LOGICAL FUNCTION SAFETY_POOL(LASTHIT)
C  Checks, if the value of any random variable X(I) falls in the safety pool.
C  If yes, corresponding limits are relaxed.
C-------------------------------------------------------------------
      INTEGER,intent(out) :: LASTHIT
      REAL(kind(1.D0)) :: Z,NEWLIMIT
      LOGICAL :: LOG1
      INTEGER :: I
      integer :: IH(CRND),IS,IL
      logical :: db=.false.
      IH=0
      LOG1=.FALSE.
! nothing to do if there are no random variables
      if (RNDLIST%DIM<1) then
            SAFETY_POOL=LOG1
            return
      endif
  ! no safety pool for non uniform sampling
      if (RNDIST.ne.rndist_uni) goto 10
      LASTHIT=0
      DO I=1,RNDLIST%DIM
      IF (RNDLIST%ACTIVE(I).GT.0) THEN
        Z=ABS(2.D0*RNDLIST%POOL(I)*XNORM(I))
        IF (Z.GT.RNDLIST%LIMITS(I)) THEN
          LASTHIT=I
          LOG1=.TRUE.
!/// when tested, replace this:
!          NEWLIMIT=RNDLIST%LIMITS(I)*RNDLIST%POOL(I)
!          if (NEWLIMIT.lt.Z) then
!            NEWLIMIT=Z*RNDLIST%POOL(I)
!          endif
!/// with this
          NEWLIMIT=MAX(Z,RNDLIST%LIMITS(I)*RNDLIST%POOL(I))
!///
          RNDIST_VOL=RNDIST_VOL*NEWLIMIT/RNDLIST%LIMITS(I)
          RNDLIST%LIMITS(I)=NEWLIMIT
          HIT(I)=HIT(I)+1
          IH(I)=1
!1         format(a,6(1x,G12.5))
!          write(*,1) 'SAFETY_POOL HIT ', trace_tot, NEUT%LABEL,I,HIT(I), NEWLIMIT, XNORM(I)
        ENDIF
      ENDIF
      enddo

      if (db.and.LOG1) then
        write(*,*) 'SAFETY_POOL HIT ',isProbe
        DO I=1,RNDLIST%DIM
          if (IH(I).gt.0) then
            call FINDSTRPAR(TRACING_VARID,'|',I,IS,IL)
            write(*,"(a8,1x,I4,3(1x,G10.2))") TRACING_VARID(IS:IS+IL-1),HIT(I),XNORM(I),RNDLIST%LIMITS(I)/2,XRND(I)
          endif
        enddo
      endif


10    SAFETY_POOL=LOG1
      END FUNCTION SAFETY_POOL


C--------------------------------------------------
      SUBROUTINE SWPOOL(ISW)
C switch safety pool off/on
C works only for temporary switch-off,
C First switch off - limits are reduced
C then switch back on - limits are restored
C--------------------------------------------------
      INTEGER,intent(in) :: ISW
      INTEGER*4 MYPOOL(CRND)
      INTEGER :: I,ISOFF
      REAL(KIND(1.D0)) :: MYLIM(CRND)
      SAVE MYLIM,MYPOOL,ISOFF
! nothing to do if there are no random variables
      if (RNDLIST%DIM<1) return
      select case (ISW)
  ! OFF
  ! remember old limits and create new ones without safety pool
  ! set limits inactive
      case(0)
        if (SAFETY_POOL_FLAG) then
          IF (.not.REPOPT%MUTE) WRITE(*,*) 'Safety pool OFF'
          DO I=1,RNDLIST%DIM
            MYLIM(I)=RNDLIST%LIMITS(I)
            MYPOOL(I)=RNDLIST%ACTIVE(I)
            IF (MYPOOL(I).NE.0) THEN
              RNDLIST%LIMITS(I)=RNDLIST%LIMITS(I)/RNDLIST%POOL(I)
              RNDLIST%ACTIVE(I)=0
            ENDIF
          enddo
          ISOFF=1
          SAFETY_POOL_FLAG=.FALSE.
          RNDIST_VOL=SAMPLING_VOLUME()
          EVENTGEN%OLDCOUNT=EVENTGEN%COUNT
          EVENTGEN%OLDSCOUNT=EVENTGEN%SCOUNT

        endif
  ! ON
  ! restore previous limits and their activity
      case default
        if (.not.SAFETY_POOL_FLAG) then
          if (ISOFF.eq.1) then
            IF (.not.REPOPT%MUTE) WRITE(*,*) 'Safety pool ON'
            DO I=1,RNDLIST%DIM
              RNDLIST%LIMITS(I)=MYLIM(I)
              RNDLIST%ACTIVE(I)=MYPOOL(I)
            enddo
            ISOFF=0
            RNDIST_VOL=SAMPLING_VOLUME()
          !  EVENTGEN%VOL=RNDIST_VOL
          endif
        ENDIF
      END select

      END SUBROUTINE SWPOOL


C-----------------------------------------------------------------------
      SUBROUTINE MAXV_UPD(ITASK)
C record max. value history for random variables
C ONLY FOR UNIFORM RANDOM NUMBERS !!
C     ITASK=0 ... clears MAXV(I) array
C     ITASK=1 ... MAXV(I) is compared with X(I) and changed if necessary
C     ITASK=2 ... limits are changed according to MAXV(I)
C-----------------------------------------------------------------------
      integer,intent(in) :: ITASK
      integer :: I
1     format(a,10(1x, G14.7))
! nothing to do if there are no random variables
      if (RNDLIST%DIM<1) then
            if (ITASK==2) RNDIST_VOL=SAMPLING_VOLUME()
            return
      endif
      if (MAXV_UPD_FLAG) then
      select case(ITASK)
! clear max. value register
      case(0)
        DO  I=1,RNDLIST%DIM
         MAXV(I)=0.D0
         !MINV(I)=0.D0
         !CTRV(I)=0.D0
        ENDDO
! update
      case(1)
        DO  I=1,RNDLIST%DIM
          IF (RNDLIST%ACTIVE(I).gt.0) THEN
            IF(ABS(XNORM(I)).GT.MAXV(I)) then
            !  write(*,1) 'MAXV_UPD HIT ', NINT(trace_tot), NEUT%LABEL ,I, XNORM(I), MAXV(I)
              MAXV(I)=ABS(XNORM(I))
            endif
            !IF(XNORM(I).LT.MINV(I)) MINV(I)=XNORM(I)
            !IF(XNORM(I).GT.MAXV(I)) MAXV(I)=XNORM(I)
          ENDIF
        ENDDO
! evaluate
      case(2)
        DO  I=1,RNDLIST%DIM
          IF (RNDLIST%ACTIVE(I).gt.0) THEN
            RNDLIST%LIMITS(I) = 2.*MAXV(I)*RNDLIST%POOL(I)*1.001
            !  write(*,1) 'MAXV_UPD UPD ', NINT(trace_tot), trace_cnt ,I, RNDLIST%LIMITS(I)
            !RNDLIST%LIMITS(I) = 2.*(MAXV(I)-MINV(I))*RNDLIST%POOL(I)*1.001
            !CTRV(I)=(MAXV(I)+MINV(I))/2.D0
          ENDIF
        ENDDO
        RNDIST_VOL=SAMPLING_VOLUME()
      END select
      endif
      END SUBROUTINE MAXV_UPD


C------------------------------------------------
      SUBROUTINE COV_CLR
C clear covariance matrix for random variables
C------------------------------------------------
      integer :: I,J
      DO I=1,CRND
        MVAL(I)=0
        DO J=1,I
          COV(J,I)=0.D0
          COV(I,J)=0.D0
        enddo
      enddo
      NCOV=0
      SCOV=0.D0
      END SUBROUTINE COV_CLR

C-----------------------------------------------------------
      SUBROUTINE COV_UPD(P)
C accumulate event in covariance matrix for random variables
C-----------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: P
      integer :: I,J
! nothing to do if there are no random variables
      if (RNDLIST%DIM<1) return      
      DO I=1,RNDLIST%DIM
        IF (RNDLIST%ACTIVE(I).gt.0) THEN
          MVAL(I)=MVAL(I)+XRND(I)*P
          DO J=1,I
            IF (RNDLIST%ACTIVE(J).gt.0) THEN
              COV(I,J)=COV(I,J)+XRND(I)*XRND(J)*P
              COV(J,I)=COV(I,J)
            ENDIF
          enddo
        ENDIF
      enddo
      NCOV=NCOV+1
      SCOV=SCOV+P
      END SUBROUTINE COV_UPD

C-----------------------------------------------------------
      SUBROUTINE COV_NEW
C evaluate covariance matrix for random variables
C and update sampling volume transformation matrix
C-----------------------------------------------------------
      INTEGER :: I,J,NROT
      REAL(KIND(1.D0)) :: CM(CRND,CRND),AUX1(CRND,CRND),VLIM(CRND)

102   FORMAT(' COV - Old volume: ',G12.6)
103   FORMAT(' COV - New volume: ',G12.6)

! we need at least 30 events to evaluate cov. matrix
      IF (NCOV.LT.30) RETURN
! nothing to do if there are no random variables
      if (RNDLIST%DIM<1) return

! set diagonal terms for inactive variables = 1
      DO I=1,RNDLIST%DIM
        IF (RNDLIST%ACTIVE(I).lt.1) CM(I,I)=1.D0
      enddo
! normalize mean and cov. matrix
      DO I=1,RNDLIST%DIM
        IF (RNDLIST%ACTIVE(I).gt.0) THEN
          TMEAN(I)=MVAL(I)/SCOV
          DO J=1,RNDLIST%DIM
            if(RNDLIST%ACTIVE(J).gt.0) CM(I,J)=COV(I,J)/SCOV
          enddo
        ENDIF
      enddo
! subtract mean
      DO i=1,RNDLIST%DIM
        DO j=1,RNDLIST%DIM
        IF ((RNDLIST%ACTIVE(I).gt.0).and.(RNDLIST%ACTIVE(j).gt.0)) THEN
          CM(I,J)=CM(I,J)-TMEAN(I)*TMEAN(J)
        ENDIF
        ENDDO
      ENDDO
! filter covariances, neglect weak correlations
c      DO I=1,RNDLIST%DIM
c        DO J=1,RNDLIST%DIM
c         IF (ABS(CM(I,J))/SQRT(CM(I,I)*CM(J,J)).LT.0.10) CM(I,J)=0.D0
c        ENDDO
c      enddo

! diagonalize cov. matrix so that AUX1(k,i)*CM(k,m)*AUX1(m,j) = delta(i,j)*VLIM(i)
      CALL JACOBI(CM,RNDLIST%DIM,CRND,VLIM,AUX1,NROT)

    !  call XML_RSXDUMP(sout,' ',1)
    !  call XML_MATRIX(sout,5,5,'COV','x|y|kx|ky|kz','x|y|kx|ky|kz',CM(1,1),CRND)
    !  call XML_MATRIX(sout,5,5,'CMAT','x|y|kx|ky|kz','x|y|kx|ky|kz',AUX1(1,1),CRND)
    !  call XML_ARRAY(sout,5,'TMEAN','x|y|kx|ky|kz',TMEAN(1))
    !  call XML_ARRAY(sout,5,'DIA','x|y|kx|ky|kz',VLIM(1))
    !  call XML_RSXDUMP(sout,' ',0)

! get new transformation matrix and limits
  ! clear TMAT
      call UNIMAT(TMAT,CRND,CRND)
  ! report old sampling volume
      if (REPOPT%VARI) write(*,102) SAMPLING_VOLUME()
  ! calculate new TMAT and limits
      call GENERATOR_SETTMAT(AUX1,VLIM,CRND,RNDLIST%DIM)
      ! clear centers
      !CTRV=0.D0
  ! report new sampling volume
      if (REPOPT%VARI) write(*,103) SAMPLING_VOLUME()
      END SUBROUTINE COV_NEW

C----------------------------------------------------------
      SUBROUTINE GENEVENT
C     filles XNORM(I) by random numbers within specified limits
C----------------------------------------------------------
      INTEGER :: I
      REAL :: Z
1     format(a,9(1x,G12.5))

      PNORM=1.D0
      select case(RNDIST)
      case(rndist_uni)
        DO I=1,RNDLIST%DIM
          XNORM(I)=RAN1()-0.5
        END DO
      case(rndist_gauss)
        DO I=1,RNDLIST%DIM
          Z=GASDEV()
          XNORM(I)=Z
          PNORM=PNORM/EXP(-0.5D0*Z**2)
        END DO
      case(rndist_cauchy)
        call RCAUCHY(XNORM,RNDLIST%DIM,CRND,PNORM)
      end select

      DO I=1,RNDLIST%DIM
        XNORM(I)=XNORM(I)*RNDLIST%LIMITS(I) !+CTRV(I)
      END DO
      !if (trace_tot==17) then
      !  write(*,1) 'GENEVENT XNORM  ',XNORM(1:RNDLIST%DIM)
      !  write(*,1) 'GENEVENT LIMITS ',RNDLIST%LIMITS(1:RNDLIST%DIM)
      !endif
      END SUBROUTINE GENEVENT

C----------------------------------------------------------
      SUBROUTINE GENPROBE
C  generate probe = ucounted event that probes outer phase space
C----------------------------------------------------------
      INTEGER :: I
1     format(a,': ',6(G10.4,1x))
      PNORM=1
  !    DO I=1,RNDLIST%DIM
  !      XNORM(I)=GASDEV()
  !    END DO
      DO I=1,RNDLIST%DIM
        XNORM(I)=RAN1()-0.5
      END DO
      DO I=1,RNDLIST%DIM
        XNORM(I)=XNORM(I)*RANG*RNDLIST%LIMITS(I)
      END DO
      END SUBROUTINE GENPROBE

      end module GENERATOR
