!//////////////////////////////////////////////////////////////////////
!////  $Id: tracingdata.f,v 1.52 2019/01/10 20:39:01 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.52 $
!////     $Date: 2019/01/10 20:39:01 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Global ray-tracing data
!////  - defines structure of the neutron event and global instances of it
!////  - defines arrays of random variables and their generator
!////
!//////////////////////////////////////////////////////////////////////
      MODULE TRACINGDATA
      use CONSTANTS
      use TRACINGOPT
      implicit none
      SAVE

! neutron event structure
      TYPE  NEUTRON; SEQUENCE
        REAL(KIND(1.D0)) :: R(3),K(3),S(3),T,PHI,P,K0 ! r,k,spin,time,phase,probability,|k|
        REAL(KIND(1.D0)) :: T0 ! time zero (for TOF measurements)
        integer :: REF(2)  ! reference to any previously stored event (storage,event)
        INTEGER :: LABEL   ! used e.g. for probing outer phase space
        integer :: CNT     ! event number (for callback reference to events storage)
        integer :: STATE   ! (0) normal (1) absorbed
      END TYPE NEUTRON
! global records used in tracing:
  ! NEUT  ... actually tracing
  ! NEUI1 ... emitted at the source
  ! NEUI  ... incident at the sample
  ! NEUF  ... scattered at the sample
  ! NEUF1 ... incident at the detector/monitor
      TYPE(NEUTRON) :: NEUT,NEUI,NEUF,NEUI1,NEUF1


! gravity
! gravity shift is then y[mm] = -0.5*gravity*t**2, where t = z/kz, z is horizontal flight path
      real(KIND(1.D0)) :: gravity=G_EARTH
      logical :: use_gravity=.false.
! gravity vector as felt by NEUT (must be transformed together with NEUT%K)
      real(KIND(1.D0)) :: NEUT_GR(3)=0.D0

! TOF data
      REAL(KIND(1.D0)) :: TOF_ZERO,TOF_SAMPLE,TOF_DETECTOR ! reference times for T0 chopper, sample and detector
      REAL(KIND(1.D0)) :: TOF_ZERO2 ! 2nd T0 reference
      logical :: TOF_ADJ_T0=.false. ! T0 is adjusted according for the passed slit at the T0 chopper
      logical :: TOF_MODE ! true if we integrate over time (time is a random variable)
      logical :: PULSED_SRC ! true if pulsed source is defined

! global tracing counters
      integer :: trace_cnt            ! passed events (counts)
      integer :: trace_clock          ! =INT(trace_tot/1000) (used to trigger some actions)
      integer :: trace_clock_old      ! old clock value, for comaprison with trace_clock in the trigger
      real(KIND(1.D0)) :: trace_sum   ! sum of passed events (weights)
      real(KIND(1.D0)) :: trace_tot   ! total number of generated events (including repeated cycles)
      real(KIND(1.D0)) :: trace_ctot  ! total number of counts (including repeated cycles)
      real(KIND(1.D0)) :: trace_crun  ! total number of counts (including repeated cycles), not cleared in twostep tracing

! setting of RAN generator
      INTEGER :: ISEED,IGEN
  ! variable for selected rnd distribution
      integer :: RNDIST=rndist_uni
  ! normalizing constant for selected distribution
      REAL(KIND(1.D0)) :: RNDIST_NORM=1.D0  ! default for uniform dist.
  ! coefficient between the sampling range and sqrt(variance)
      REAL(KIND(1.D0)) :: RNDIST_RANGE=7.D0 ! default for uniform dist.

! variables controling event flow
      INTEGER :: TCOUNT_REQ            ! count target requested by user
      integer :: TCOUNT                ! count target current loop
      LOGICAL :: TIMEOUT               ! flag for timeout
      integer :: CNT_VARI              ! count target for variance reduction
      real(KIND(1.D0)) :: CNT_TOUT     ! count target for timeout
! flag indicates current up-stream tracing, set automatically before tracing
      logical :: TRACING_UP
      logical :: NEUT_DBG=.false.

! source multiplication
      logical :: SRC_OVERLAP=.false.
      integer :: SRC_NPLS=3
      REAL(KIND(1.D0)) :: SRC_PERIOD   ! source period in us*h_bar/m

! dimension of arrays with random variables
      INTEGER,parameter :: CRND=64
! array of random numbers and limits setting for NESS
      TYPE RANDFIELD; SEQUENCE
         REAL(KIND(1.D0)) :: MEAN(CRND),LIMITS(CRND),POOL(CRND)
         INTEGER ::  ACTIVE(CRND),DIM
      END TYPE RANDFIELD
      INTEGER :: NSEED,HIT(CRND),IXRND(CRND)
      REAL(KIND(1.D0)) :: XNORM(CRND),XRND(CRND),PNORM
      TYPE(RANDFIELD) :: RNDLIST
  ! actual sampling volume
      real(KIND(1.D0)) :: RNDIST_VOL

!  transformation arrays for sampling optimization
      REAL(KIND(1.D0)) :: TMAT(CRND,CRND),TMEAN(CRND)

! safety pool flag
      LOGICAL :: SAFETY_POOL_FLAG
      logical :: isProbe=.false. ! probe flag
! switch and target count for limits update based on maximum values
      logical :: MAXV_UPD_FLAG
      integer :: MAXV_UPD_COUNT

! results
      real(KIND(1.D0)) :: TRACING_INT             ! intensity, using factor 10^14/100 (source flux is in [10^14/cm^2/s])
      real(KIND(1.D0)) :: TRACING_DINT            ! error of INT
      real(KIND(1.D0)) :: TRACING_TOF             ! mean flight time
      real(KIND(1.D0)) :: TRACING_DTOF            ! error of TOF
      real(KIND(1.D0)) :: TRACING_DE,TRACING_DDE  ! energy spread and its error
      real(KIND(1.D0)) :: TRACING_PRECISION       ! dInt/Int
      real(KIND(1.D0)) :: FLUX_NORM               ! events normalization factor
      real(KIND(1.D0)) :: PRIMARY_INT             ! intensity for the primary beam

      real(KIND(1.D0)) :: TRACING_VOL             ! sampling volume per event
      real(KIND(1.D0)) :: TRACING_DVOL            ! error of TRACING_VOL

! during tray-tracing, we need to know where to store events
! i.e. we need to know current subset and dataset represented by the beamline

      integer :: CSUBSET(2)    ! current subset index
      integer :: CDATASET(2)   ! current dataset index

      logical :: STOPTRACING=.false. ! signal to stop tracing, used by GUI
      logical :: STOPTIME=.false. ! signal to stop tracing, when DT<0

      character(128) :: TRACING_VARID  ! ID strings for random variables


      real(KIND(1.D0)) :: modulation_dt,modulation_d0


      contains

!-------------------------------------------------------------
      SUBROUTINE incCounter(cnt)
!-------------------------------------------------------------
        integer,intent(out) :: cnt
        if (.not.isProbe) cnt=cnt+1
      end SUBROUTINE incCounter

!-------------------------------------------------------------
      subroutine TransportK(K,DT)
! transport neutron as a particle with velocity hk/m by DT
! DT is always in units of m_n/h_bar
!-------------------------------------------------------------
        real(kind(1.D0)) :: K(3),DT
        integer :: I
        if (DT.eq.0.D0) return
        do i=1,3
          NEUT%R(i)=NEUT%R(i)+DT*K(i)
        enddo
        if (TRACING_UP) then
          NEUT%T=NEUT%T-DT! /HOVM
        else
          NEUT%T=NEUT%T+DT! /HOVM
        endif

      end SUBROUTINE TransportK

!-------------------------------------------------------------
      subroutine Transport(DT)
! transport neutron as free particle by DT
! DT is always in units of m_n/h_bar
!-------------------------------------------------------------
        real(kind(1.D0)) :: DT
        integer :: I
        if (DT.eq.0.D0) return
        do i=1,3
          NEUT%R(i)=NEUT%R(i)+DT*NEUT%K(i)
        enddo
        if (TRACING_UP) then
          NEUT%T=NEUT%T-DT! /HOVM
        else
          NEUT%T=NEUT%T+DT! /HOVM
        endif

      end SUBROUTINE Transport

!-------------------------------------------------------------
      subroutine TransportG(DT)
! transport neutron in gravity field
! DT is always in units of m_n/h_bar
!-------------------------------------------------------------
        real(kind(1.D0)) :: DT
        integer :: i
        if (DT.eq.0.D0) return
        if (.not.use_gravity) then
          call Transport(DT)
          return
        endif
        do i=1,3
          NEUT%R(i)=NEUT%R(i)+DT*NEUT%K(i)+0.5D0*NEUT_GR(i)*DT**2
          NEUT%K(i)=NEUT%K(i)+NEUT_GR(i)*DT
        enddo
        if (TRACING_UP) then
          NEUT%T=NEUT%T-DT! /HOVM
        else
          NEUT%T=NEUT%T+DT! /HOVM
        endif

      end SUBROUTINE TransportG



!------------------------------------------------------
      SUBROUTINE GENERATE_NOM(K0,NEU)
! generate nominal event at the source with given |K|
! for tracing beam axis
!-------------------------------------------------------
      real(kind(1.D0)),intent(in) :: K0
      TYPE(NEUTRON),intent(out) :: NEU
      NEU%R=0.D0
      NEU%K=0.D0
      NEU%K(3)=K0
      NEU%K0=K0
      NEU%S=0.D0
      NEU%S(1)=1.D0
      NEU%T=0.D0
      NEU%PHI=0.D0
      NEU%P=1.D0
      end SUBROUTINE GENERATE_NOM

!---------------------------------------------------
      real(KIND(1.D0)) function GETVAR(IVAR,NEU,K0)
! Get value of the neutron phase-space variable defined by index=IVAR
!---------------------------------------------------
      integer,intent(in) :: IVAR
      type(NEUTRON),intent(in) :: NEU
      real(KIND(1.D0)),intent(in) :: K0
      real(KIND(1.D0)) :: XVAL
        select case(IVAR)
        case(1,2,3)
          XVAL=NEU%R(IVAR)
        case(4,5)
          XVAL=NEU%K(IVAR-3)/NEU%K0
        case(6)
        ! dkz/|k| refers to the nominal |k|, not the current event
          XVAL=NEU%K(IVAR-3)/K0
        case(7)
          XVAL=NEU%T/HOVM
        case(8)
          XVAL=HSQOV2M*NEU%K0**2
        case(9)
          XVAL=2.D0*PI/NEU%K0
        case default
          XVAL=NEU%R(1)
        end select
        GETVAR=XVAL
      end function GETVAR


!---------------------------------------------------
      subroutine GETVARSTR(IVAR,SARG)
! get variable name by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      character(*),intent(out) :: SARG
      character(1),parameter :: CH(6)=(/'X','Y','Z','x','y','z'/)
      if (len(SARG).ge.16) then
        select case(IVAR)
        case(1,2,3)
          SARG=CH(IVAR)//' [mm]'
        case(4,5,6)
          SARG='d_'//CH(IVAR)//'/k'
        case(7)
          SARG='TOF [us]'
        case(8)
          SARG='E [meV]'
        case(9)
          SARG='lambda [A]'
        case default
          SARG='undefined'
        end select
      endif
      end subroutine GETVARSTR



      end module TRACINGDATA

