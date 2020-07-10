!//////////////////////////////////////////////////////////////////////
!////  $Id: rndgen.f,v 1.8 2019/08/13 10:02:33 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2013, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.8 $
!////     $Date: 2019/08/13 10:02:33 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Random generators
!////
!//////////////////////////////////////////////////////////////////////
      MODULE RNDGEN
      use SIMATH
      USE TRACINGDATA
      USE mtmod
      implicit none
      SAVE
      PRIVATE

      integer :: gasdev_set=0
      integer :: RAN1NR_INI=0
      integer :: RAN1NR_IX1,RAN1NR_IX2,RAN1NR_IX3
      REAL :: RAN1NR_R(97)
      REAL(KIND(1.D0)) :: gset

      public RCAUCHY,GASDEV1,GASDEV,GASDEV2,GENINT
      public RAN1SEED,RAN1,RAN1NR,RAN1TEST,M16XV16

      logical :: rndgen_dbg=.false.
      integer :: dbg_cnt=0
      public rndgen_dbg

      contains

!-----------------------------------
      SUBROUTINE M16XV16(IT,M,A,B,C)
!-----------------------------------
      REAL(KIND(1.D0)),intent(in) :: A(CRND,CRND),B(CRND)
      REAL(KIND(1.D0)),intent(out) :: C(CRND)
      INTEGER,intent(in) :: M,IT
      INTEGER :: J,I
      DO 10 J=1,M
         C(J)=0.
         IF (IT.GT.0) THEN
           DO 20 I=1,M
20         C(J)=C(J)+A(J,I)*B(I)
         ELSE
           DO 30 I=1,M
30         C(J)=C(J)+A(I,J)*B(I)
         ENDIF
10    CONTINUE
      END SUBROUTINE M16XV16


!----------------------------------------------
! random points on n-dim. unit sphere surface
      SUBROUTINE RSPHERE(X,NX,MX)
!----------------------------------------------
      integer, intent(in) :: NX,MX
      REAL(KIND(1.D0)),intent(out) :: X(MX)
      real,parameter :: EPS=0.001
      real,parameter :: CPI=3.14159265358979
      REAL(KIND(1.D0)) :: Z
      integer :: i
      Z=0.0
      do i=1,NX
        X(i)=GASDEV()
        Z=Z+X(i)**2
      enddo
      do i=1,NX
        X(i)=X(i)/SQRT(Z)
      enddo
      END SUBROUTINE RSPHERE

!----------------------------------------------
! rnd numbers with Cauchy distribution ***
      REAL(KIND(1.D0)) FUNCTION GENCAUCHY()
! limited to ~ (-1/EPS,+1/EPS)
!----------------------------------------------
      REAL(KIND(1.D0)),parameter :: EPS=0.001D0
      REAL(KIND(1.D0)),parameter :: CPI=3.14159265358979D0
      REAL(KIND(1.D0)) :: Z
      Z=0.0
      DO WHILE ((Z.lt.EPS).or.(1.0-Z.lt.EPS))
        Z=RAN1()
      enddo
      GENCAUCHY=1.0/TAN(Z*CPI)
      END FUNCTION GENCAUCHY

!----------------------------------------------
! generate random integer in given range
      integer FUNCTION GENINT(IMIN,IMAX)
!----------------------------------------------
      integer, intent(in) :: IMIN,IMAX
      integer :: I
      REAL(kind(1.D0)) :: Z
      I=IMIN-1
      if (IMAX.GE.IMIN) then
      DO WHILE ((I.LT.IMIN).or.(I.GT.IMAX))
        Z=RAN1()
        I=NINT(IMIN+Z*(IMAX-IMIN+1)-0.5D0)
      enddo
      endif
      GENINT=I
      END FUNCTION GENINT

!----------------------------------------------
! multidimensional Cauchy distribution ***
      SUBROUTINE RCAUCHY(X,NX,MX,PX)
! isotropic nx-dimensional vectors with Cauchy |X| distribution
! generator works with REAL, but X is KIND(1.D0) !!
! limited to ~ (-1/EPS,+1/EPS)
!----------------------------------------------
      integer, intent(in) :: NX,MX
      REAL(KIND(1.D0)),intent(out) :: X(MX),PX
      real,parameter :: EPS=1.0E-8
      real,parameter :: CPI=3.14159265358979
      REAL(KIND(1.D0)),parameter :: CSQRTPI=1.772453850905516D0
      REAL(KIND(1.D0)) :: SN,R
      integer :: i
      call RSPHERE(X,NX,MX)
      R=0.0
      do while (R.lt.EPS)
        R=ABS(GENCAUCHY())
      ENDDO
      do i=1,NX
        X(i)=X(i)*R
      enddo
      PX=CPI*(1.D0+R**2)
    ! n-dim sphere surface
!      write(*,*) 'R,NX,GAMMAN2(NX)',R,NX,GAMMAN2(NX)
      SN=(CSQRTPI*R)**NX/R/GAMMAN2(NX)
!      write(*,*) 'PX,SN',PX,SN
      PX=PX/SN
      end SUBROUTINE RCAUCHY


!----------------------------------------------
! GASDEV with centre and limits
! scaled in units sigma=1
      REAL(kind(1.D0)) FUNCTION GASDEV1(centre,limits)
!----------------------------------------------
      REAL(kind(1.D0)), intent(in) :: limits,centre
      REAL(kind(1.D0)) :: Z, zlim
      REAL(kind(1.D0)), parameter :: eps=1.D-8
      zlim = abs(limits)
      if (zlim<eps) then
        Z = 0.D0
      else
        Z=GASDEV()
        DO WHILE (ABS(Z).GT.zlim)
          Z=GASDEV()
        enddo
      endif
      GASDEV1=Z+centre
      END FUNCTION GASDEV1

!----------------------------------------------
! GASDEV on interval (xmin, INF)
! scaled in units sigma=1
      REAL(kind(1.D0)) FUNCTION GASDEV2(xmin)
!----------------------------------------------
      REAL(kind(1.D0)), intent(in) :: xmin
      REAL(kind(1.D0)) :: Z
      Z=GASDEV()
      DO WHILE (Z<xmin)
        Z=GASDEV()
      enddo
      GASDEV2=Z
      END FUNCTION GASDEV2

!----------------------------------------------
! Generation of Gaussian deviates by GASDEV from Numerical Recipes ***
      REAL(kind(1.D0)) FUNCTION GASDEV()
!----------------------------------------------
      REAL(kind(1.D0)) :: R, V1,V2,FAC
2     format(a,': ', 8(1x, G14.7))
      IF (gasdev_set.EQ.0) THEN
1       V1=2.D0*RAN1()-1.D0
        V2=2.D0*RAN1()-1.D0
        R=V1**2+V2**2
        IF(R.GE.1.)GO TO 1
        FAC=SQRT(-2.*LOG(R)/R)
        GSET=V1*FAC
        GASDEV=V2*FAC
        gasdev_set=1
        if (rndgen_dbg) write(*,2) 'GASDEV gen',V1, V2, FAC, V2*FAC, GSET
      ELSE
        GASDEV=GSET
        gasdev_set=0
        if (rndgen_dbg) write(*,2) 'GASDEV GSET',GSET
      ENDIF

      END FUNCTION GASDEV

!----------------------------------------------
      SUBROUTINE RAN1SEED(IDUM)
! Initialize random number generator
!----------------------------------------------
      INTEGER,intent(in) :: IDUM
      REAL :: DUM,Z
      REAL :: RAND,SECNDS

      if (rndgen_dbg) write(*,*) 'RAN1SEED ',IDUM
! Generate ISEED from the system time on the first call
      IF (IDUM.EQ.0)  then
         Z=SECNDS(0.0)
         ISEED=2*INT(10000.+Z)+1
      ENDIF
! If argument<>0, use it as a new seed
      IF (IDUM.NE.0) ISEED=ABS(IDUM)
      write(*,*) 'SEED = ',ISEED,' ',IGEN
! Initialize required generator
      IF (IGEN.EQ.0) call sgrnd(iseed)  ! Mersenne Twister
      IF (IGEN.EQ.1) DUM=RAN1NR(-iseed) ! Numerical Recipes RAN1
      IF (IGEN.EQ.2) DUM=RAND(iseed)     ! System generator
      gasdev_set=0
      END SUBROUTINE RAN1SEED


!----------------------------------------------
      REAL(kind(1.D0)) FUNCTION RAN1()
! Call random number generator
! generate uniform random numbers in the interval [1e-8..1-1e-8]
!----------------------------------------------
      REAL(kind(1.D0)) :: Z
      REAL :: RAND
      REAL(kind(1.D0)),PARAMETER :: EPS=1.0D-8
      REAL(kind(1.D0)),PARAMETER :: EPS1=1.D0-1.0D-8

      dbg_cnt=dbg_cnt+1
      if (rndgen_dbg) write(*,*) 'RAN1 IGEN=',IGEN,dbg_cnt
10    IF (IGEN.EQ.0) THEN ! Mersenne Twister
        Z=GRND()
        if (rndgen_dbg) write(*,*) 'RAN1 ',Z
      ELSE IF (IGEN.EQ.1) THEN  ! Numerical Recipes RAN1
        Z=1.D0*RAN1NR(ISEED)
      ELSE IF (IGEN.EQ.2) THEN ! System generator
        Z=1.D0*RAND(ISEED)
      ENDIF


!      write(*,*) ISEED,' ',Z
!      pause
      IF (Z.LT.EPS) GOTO 10
      IF (Z.GT.EPS1) GOTO 10
      RAN1=Z
      END FUNCTION RAN1


!----------------------------------------------
      REAL FUNCTION RAN1NR(IDUM)
! *** Random number generator from Numerical Recipes (RAN1): ***
!----------------------------------------------
!      implicit real*4 (a-h,o-z)
!      implicit integer*4 (i-n)

      INTEGER,intent(in) :: IDUM
      INTEGER :: J
      integer,parameter :: M1=259200
      integer,parameter :: M2=134456
      integer,parameter :: M3=243000
      integer,parameter :: IA1=7141
      integer,parameter :: IA2=8121
      integer,parameter :: IA3=4561
      integer,parameter :: IC1=54773
      integer,parameter :: IC2=28411
      integer,parameter :: IC3=51349
      REAL,parameter :: RM1=3.8580247E-6
      REAL,parameter :: RM2=7.4373773E-6

      IF ((IDUM.LT.0).OR.(RAN1NR_INI.EQ.0)) THEN
        ISEED=IDUM
        RAN1NR_INI=-1
        RAN1NR_IX1=MOD(IC1-ISEED,M1)
        RAN1NR_IX1=MOD(IA1*RAN1NR_IX1+IC1,M1)
        RAN1NR_IX2=MOD(RAN1NR_IX1,M2)
        RAN1NR_IX1=MOD(IA1*RAN1NR_IX1,M2)
        RAN1NR_IX1=MOD(IA1*RAN1NR_IX1+IC1,M1)
        RAN1NR_IX3=MOD(RAN1NR_IX1,M3)
        DO 11 J=1,97
          RAN1NR_IX1=MOD(IA1*RAN1NR_IX1+IC1,M1)
          RAN1NR_IX2=MOD(IA2*RAN1NR_IX2+IC2,M2)
          RAN1NR_R(J)=(FLOAT(RAN1NR_IX1)+FLOAT(RAN1NR_IX2)*RM2)*RM1
11      CONTINUE
        ISEED=1
      ENDIF
12    RAN1NR_IX1=MOD(IA1*RAN1NR_IX1+IC1,M1)
      RAN1NR_IX2=MOD(IA2*RAN1NR_IX2+IC2,M2)
      RAN1NR_IX3=MOD(IA3*RAN1NR_IX3+IC3,M3)

      J=1+(97*RAN1NR_IX3)/M3
      IF(J.GT.97.OR.J.LT.1) THEN
        write(*,*) 'RAN1 error:'
        WRITE(*,*) J, 'IX = ',RAN1NR_IX1, RAN1NR_IX2, RAN1NR_IX3, 'ISEED = ', ISEED
      END IF
      RAN1NR=RAN1NR_R(J)
      RAN1NR_R(J)=(FLOAT(RAN1NR_IX1)+FLOAT(RAN1NR_IX2)*RM2)*RM1
      RETURN
      END FUNCTION RAN1NR

!----------------------------------------------
!  Test Random number covariances
!----------------------------------------------
      SUBROUTINE RAN1TEST(N,NOE)
      integer,intent(in) :: N,NOE
      real :: SECNDS
      integer :: i,j,M,nc,NMAX,IM,JM
      real :: C(128,128),V(128),MEAN(128),RW,C0,MAX,DISP,S,S2,Z,T1,T2
      real :: W,W2,Z0
      M=N
      if (M.LT.2) M=2
      if (M.GT.128) M=128
      DO i=1,M
        MEAN(i)=0
        DO j=1,M
          C(I,J)=0
        ENDDO
      ENDDO
      NMAX=NOE/M
      DO nc=1,NMAX
         DO I=1,M
               V(I)=RAN1()-0.5
         END DO
         DO I=1,M
         MEAN(I)=MEAN(I)+V(I)
         DO J=1,M
            C(I,J)=C(I,J)+V(I)*V(J)
         END DO
         END DO
      END DO
1     format(E8.2,' ',$)
      MAX=0.
      write(*,*)
      write(*,*) 'Covariances:'
      DO I=1,M
      DO J=1,M
         C0=MEAN(I)*MEAN(J)/NMAX/NMAX
         if (I.EQ.J) C0=C0+1./12
         write(*,1) (C(I,J)/NMAX-C0)*12.
         if(ABS(C(I,J)/NMAX-C0).GT.MAX) THEN
            MAX=ABS(C(I,J)/NMAX-C0)
            IM=I
            JM=J
         endif
      END DO
      write(*,*)
      END DO
      write(*,*)
      write(*,*) 'Mean:'
      DO I=1,M
        write(*,1) MEAN(I)/NMAX
      END DO
      write(*,*)
      S=0
      S2=0
      W=0
      W2=0
      T1=SECNDS(0.0)
      DO I=1,10000
        RW=0.
        Z0=0.
        DO J=1,10000
          Z=RAN1()-0.5
                 IF (Z.GT.0) THEN
             RW=RW+1.
          ELSE
             RW=RW-1.
          ENDIF
          Z0=Z0+Z
        ENDDO
        S=S+RW
        S2=S2+RW**2
        W=W+Z0/10000
        W2=W2+(Z0/10000)**2
      ENDDO
      T2=SECNDS(0.0)
      DISP=W2/10000-(W/10000)**2
      write(*,*)
4     format('Variance of mean value: ',G12.6)
      write(*,4) SQRT(DISP)/SQRT(1./12/(10000-1))
2     format('Variance of discrete random walk : ',G12.6)
      write(*,2) (S2/10000-(S/10000)**2)/10000
3     format('Speed: ',G12.5,'/msec')
      write(*,3) 1.e8/(T2-T1)/1000
      END SUBROUTINE RAN1TEST


      end module RNDGEN









