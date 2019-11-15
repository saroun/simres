C------------------------------------------------------------------------
      LOGICAL*4 FUNCTION GUIDE_RUN(OBJ,NCNT,tof0)
C run a single event for OBJ=BENDER
C NCNT: counter of transmitted events
C tof0: nominal time-of-flight (constant for all events)
C------------------------------------------------------------------------
      USE COMPONENTS
      USE TRACING
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'randvars.inc'
      INCLUDE 'ness_common.inc'
! arguments
      TYPE(BENDER) :: OBJ
      INTEGER :: NCNT
      REAL(KIND(1.D0)),intent(in) :: tof0
! functions
      LOGICAL :: SAFETY_POOL,BENDER_GO
! locals
      LOGICAL :: LOG1
      REAL :: RAN1
      REAL(KIND(1.D0)) :: DT
      INTEGER :: i,idir,ihit
      TYPE(NEUTRON) :: NEU

! generate array of random numbers
      call RANDFILL
      CALL M16XV16(-1,RNDLIST.DIM,TMAT,XNORM,XRND)
      do i=1,RNDLIST.DIM
        XRND(i)=XRND(i)+TMEAN(I)
      enddo

! fill NEU
      NEU.R(1)=XRND(1)
      NEU.R(2)=XRND(2)
      NEU.R(3)=0.0

c        NEU.K0=K0+(RAN1()-0.5)*DK
c        C1=SQRT(ABS(1-XRND(1)**2))
c        C2=SQRT(ABS(1-XRND(2)**2))
c        NEU.K(1)=NEU.K0*XRND(1)*C2
c        NEU.K(2)=NEU.K0*XRND(2)
c        NEU.K(3)=NEU.K0*C1*C2

      NEU.K(1)=XRND(3)
      NEU.K(2)=XRND(4)
      NEU.K(3)=XRND(5)
      NEU.K0=SQRT(NEU.K(1)**2+NEU.K(2)**2+NEU.K(3)**2)

      NEU.S(1:3)=0.D0
      NEU.S(1)=2*NINT(RAN1())-1
      NEU.T=0.D0
      NEU.P=1.D0

! transform to previous object ref. frame
! move back by DIST, because xxx_GO procedure does the reverse
c        NEU.R(3)=NEU.R(3)+OBJ.FRAME.DIST
c        dt=-NEU.R(3)/NEU.K(3)
c        DO i=1,3
c          NEU.R(i)=NEU.R(i)+dt*NEU.K(i)
c        enddo
c        NEU.T=NEU.T+dt/HOVM

! pass through OBJ
      NEUI=NEU

      LOG1=BENDER_GO(OBJ,NEU)

      if (.NOT.LOG1) goto 10

! check pool
      IF (SAFETY_POOL(ihit)) THEN
        NCNT=0
        OBJ.FRAME.COUNT=0
        goto 10
      endif

c      call WrtNeu(NEUI)
c      write(*,*) OBJ.FRAME.COUNT,LOG1
!      read(*,*)


! transfer to the next time frame
      IDIR=1
      if (TRACING_UP) IDIR=-1
      DT=tof0-NEU.T*HOVM
      DO i=1,3
        NEU.R(i)=NEU.R(i)+IDIR*DT*NEU.K(i)
      enddo
      NEU.T=NEU.T+IDIR*DT/HOVM
      NEUI.P=NEU.P
      NEUF=NEU
      NCNT=NCNT+1
      LOG1=.TRUE.
c      call WrtNeu(NEUF)
c110   format(a20,4(G12.5,2x))
c      write(*,110) 'NEUI: ',NEUI.R(1),NEUI.K(1),NEUI.T
c      write(*,110) 'NEUF: ',NEUF.R(1),NEUF.K(1),NEUF.T,tof0/hovm
c      write(*,110) 'DIFF: ',NEUF.R(1)-NEUF.K(1)*tof0,NEUF.R(3)-NEUF.K(3)*tof0,NEUF.R(3)
c      write(*,*)

10    GUIDE_RUN=LOG1

      END FUNCTION

C------------------------------------------------------------------------
      SUBROUTINE GUIDE_MAP(OBJ,K0,DK,DIM,MCOVAR,CTRANS,MEAN)
C get resolution matrix, transport matrix and mean of random variables for a guide
C random variables are:
C XRND(1) ... SIN(AH) ... AH is horizontal divergence
C XRND(2) ... SIN(AV) ... AV is horizontal divergence
C XRND(3) ... X, horizontal position
C XRND(4) ... Y, vertical position
C------------------------------------------------------------------------
      USE COMPONENTS
      USE COVARIANCES,only:STAT_INI,STAT_INP,STAT_COVAR,STAT_CONV
      USE TRACING
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'randvars.inc'
      INCLUDE 'ness_common.inc'

! arguments
      TYPE(BENDER) :: OBJ
      integer,intent(in) :: DIM
      REAL(KIND(1.D0)),intent(in) :: K0,DK
      REAL(KIND(1.D0)),intent(out) :: MCOVAR(DIM,DIM),CTRANS(DIM,DIM),MEAN(DIM)
! dimension
      integer, parameter :: ND = 4
! number of trials
      integer, parameter :: NC = 1000

      integer :: i,ic
      REAL(KIND(1.D0)) :: X(ND),NTOT,tof0
      REAL(KIND(1.D0)) :: BENDER_DIV
      LOGICAL :: GUIDE_RUN

! default rand. var. pool settings
      RNDLIST.DIM=5
      DO I=1,RNDLIST.DIM
        RNDLIST.MEAN(I)=0.D0
        RNDLIST.POOL(I)=1.1D0
        RNDLIST.ACTIVE(I)=1
        RNDLIST.LIMITS(I)=1.D0
      enddo

! set sampling limits
      RNDLIST.LIMITS(1)=OBJ.FRAME.SIZE(1)*RNDLIST.POOL(1)
      RNDLIST.LIMITS(2)=OBJ.FRAME.SIZE(2)*RNDLIST.POOL(2)
      RNDLIST.LIMITS(3)=k0*BENDER_DIV(OBJ,2*PI/K0,1)*RNDLIST.POOL(1)
      RNDLIST.LIMITS(4)=k0*BENDER_DIV(OBJ,2*PI/K0,2)*RNDLIST.POOL(2)
      RNDLIST.LIMITS(5)=DK

      RNDLIST.ACTIVE(5)=0
      RNDLIST.POOL(5)=1.D0

! set transformation matrix
      DO I=1,RNDLIST.DIM
        TMAT(1:CRND,I)=0.D0
        TMAT(I,I)=1.D0
      enddo
      TMEAN(1:CRND)=0.D0
      TMEAN(5)=K0

! nominal TOF from previous component:
      tof0=OBJ.FRAME.DIST/K0

! reset event counters
      ic=0
      NTOT=0

!-------------- TRACING CYCLE -------------------
      do while (ic.lt.NC)
        if (ic.eq.0) then
      ! initialize covariance data
          call STAT_INI(1,ND)
          call STAT_INI(2,ND)
        endif
        NTOT=NTOT+1
100     format('.',$)
        if (MOD(INT(NTOT),100000).eq.0) write(*,100)
        if (GUIDE_RUN(OBJ,ic,tof0)) then

! get array with random numbers for transmitted event
c        X(2)=NEUF.K(2)/NEUF.K0
c        C2=SQRT(ABS(1-X(2)**2))
c        X(1)=NEUF.K(1)/NEUF.K0/C2
c        X(3)=NEUF.R(1)
c        X(4)=NEUF.R(2)
          X(1)=NEUF.K(1)
          X(2)=NEUF.K(2)
          X(3)=NEUF.R(1)
          X(4)=NEUF.R(2)

! accumulate covariance and mean
c100       format(a20,4(G12.5,2x))
c          write(*,100) ' ',ic
c          write(*,100) 'XRND ',XRND(1:ND)
c          write(*,100) 'X    ',X(1:ND)
          call STAT_INP(1,ND,XRND(1),NEUF.P)
          call STAT_INP(2,ND,X(1),NEUF.P)
        endif
      enddo
!-------------- END OF TRACING CYCLE ---------------
      write(*,*)

! evaluate cov. matrix 1
      call STAT_COVAR(1,1,DIM,MCOVAR,MEAN)

! evaluate transport matrix
      call STAT_CONV(CTRANS,DIM)

1     format('conts: ',I8,'  total: ',G12.5,'  efficiency: ',G12.5)
      write(*,1) ic,ntot,ic/ntot
      end

C------------------------------------------------------------------------
      SUBROUTINE GUIDE_BENCH(OBJ,K0,DK,DIM,MCOVAR,CTRANS,MEAN)
C run tracing through a single component (GUIDE)
C random variables are:
C XRND(1) ... SIN(AH) ... AH is horizontal divergence
C XRND(2) ... SIN(AV) ... AV is horizontal divergence
C XRND(3) ... X, horizontal position
C XRND(4) ... Y, vertical position
C------------------------------------------------------------------------
      USE COMPONENTS
      USE COVARIANCES,only:STAT_INI,STAT_INP,STAT_COVAR,STAT_CONV
      USE TRACING
      implicit none
      INCLUDE 'const.inc'
      INCLUDE 'randvars.inc'
      INCLUDE 'ness_common.inc'

! arguments
      TYPE(BENDER) :: OBJ
      integer,intent(in) :: DIM
      REAL(KIND(1.D0)),intent(in) :: K0,DK
      REAL(KIND(1.D0)),intent(in) :: MCOVAR(DIM,DIM),CTRANS(DIM,DIM),MEAN(DIM)
! dimension
      integer, parameter :: ND = 4
! number of trials
      integer, parameter :: NC = 10000

      integer :: i,j,ic
      REAL(KIND(1.D0)) :: NTOT,tof0
      LOGICAL :: GUIDE_RUN,LOG1

      REAL(KIND(1.D0)) :: AUX(ND,ND),U(ND,ND),D(ND)

! default rand. var. pool settings
      RNDLIST.DIM=5
      DO I=1,RNDLIST.DIM
        RNDLIST.MEAN(I)=0.D0
        RNDLIST.POOL(I)=1.1D0
        RNDLIST.ACTIVE(I)=1
        RNDLIST.LIMITS(I)=1.D0
      enddo

      TMEAN(5)=K0
      RNDLIST.LIMITS(5)=DK
      RNDLIST.ACTIVE(5)=0
      RNDLIST.POOL(5)=1.D0

! calculate new TMAT and limits
      CALL JACOBI(MCOVAR,AUX,DIM,ND,D,U,I)

      DO I=1,RNDLIST.DIM
      IF (RNDLIST.ACTIVE(I).GT.0) THEN
        RNDLIST.LIMITS(I)=2*SQRT(6*ABS(D(I)))
        DO J=1,RNDLIST.DIM
          IF (RNDLIST.ACTIVE(J).GT.0) THEN
             TMAT(I,J)=U(J,I)
          ENDIF
        ENDDO
        TMEAN(I)=MEAN(I)
      ENDIF
      ENDDO

! nominal TOF from previous component:
      tof0=OBJ.FRAME.DIST/K0

! reset event counters
      ic=0
      NTOT=0

! TRACING CYCLE
      do while (ic.lt.NC)
        NTOT=NTOT+1
100     format('.',$)
        if (MOD(INT(NTOT),100000).eq.0) write(*,100)
        LOG1 = GUIDE_RUN(OBJ,ic,tof0)
      enddo
! END OF TRACING CYCLE
      write(*,*)
1     format('conts: ',I8,'  total: ',G12.5,'  efficiency: ',G12.5)
      write(*,1) ic,ntot,ic/ntot

      end SUBROUTINE

