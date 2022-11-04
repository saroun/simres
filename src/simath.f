!//////////////////////////////////////////////////////////////////////
!////  $Id: simath.f,v 1.35 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.35 $
!////     $Date: 2019/08/16 17:16:25 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Math subroutines
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SIMATH
      use CONSTANTS

      implicit none

      interface INVERT
        module procedure INVERT_R4
        module procedure INVERT_R8
      end interface

      private DIAG
      contains

!------------------------------------------
      SUBROUTINE PRN_VECTOR(TITLE,V)
! print 3-vector
!------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: V(3)
      character*(*),intent(in) :: TITLE
1     format(a,':')
2     format(4x,G12.5)
      write(*,1) trim(TITLE)
      write(*,2) V(1:3)
      end SUBROUTINE PRN_VECTOR

!------------------------------------------
      SUBROUTINE PRN_MATRIX(TITLE,R)
! print 3x3 matrix
!------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: R(3,3)
      character*(*),intent(in) :: TITLE
      integer :: i
1     format(a,':')
2     format(3x,3(1x,G12.5))
      write(*,1) trim(TITLE)
      do i=1,3
        write(*,2) R(i,1:3)
      enddo
      end SUBROUTINE PRN_MATRIX

!------------------------------------------
      SUBROUTINE NORM(R)
! normalize given 3-vector
!------------------------------------------
      REAL(KIND(1.D0)),intent(inout) :: R(3)
      REAL(KIND(1.D0)) :: Z
      integer :: i
      Z=sqrt(R(1)**2+R(2)**2+R(3)**2)
      if (Z>1.D-12) then
        do i=1,3
          R(i)=R(i)/Z
        enddo
      endif
      end SUBROUTINE NORM


!------------------------------------------
      real(kind(1.D0)) function  DW_PHI1(x)
! calculates the integral used in Debye-Waller factor
! Integral[z/(exp(z)-1)]dz from 0 to x
! Simpson/3 rule
!------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: x
      REAL(KIND(1.D0)) :: Z,ksi,dx
      integer :: i
      REAL(KIND(1.D0)),parameter :: P(2)=(/4.D0,2.D0/)
      if (x.gt.5) then
    ! use asymptotic approximation
        DW_PHI1=(PI**2)/6.D0 - exp(-x)/(x+1)
      else
    ! integrate by Simpson/3
        Z=1.D0+x/(exp(x)-1.D0)
        dx=x/100
        DO i=2,100
          ksi=(i-1)*dx
          Z=Z+P(MOD(i,2)+1)*ksi/(exp(ksi)-1.D0)
        ENDDO
        DW_PHI1=Z*dx/3.D0
      endif
      end function DW_PHI1

!------------------------------------------
      real(kind(1.D0)) function  DW_PHI(x)
! calculates the integral used in Debye-Waller factor
! Integral[z/(exp(z)-1)]dz from 0 to x
! Bernoulli formula
!------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: x
      REAL(KIND(1.D0)) :: Z,B(0:30),Xn
      integer :: i
      REAL(KIND(1.D0)),parameter :: P(2)=(/4.D0,2.D0/)
      if (x.gt.6.D0) then
    ! use asymptotic approximation
        Z=(PI**2)/6.D0 - exp(-x)/(x+1)
      else
    ! integral using Bernoulli numbers
        call BERNOULLI(30,B)
        Z=0.D0
        Xn=x
        DO I=0,30
          Z=Z+B(i)*Xn
          Xn=Xn*x/(I+2)
        enddo
      endif
      DW_PHI=Z
      end function DW_PHI

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION FACT(N)
C N factorial
C------------------------------------------
      integer, intent(in) :: N
      integer :: i
      REAL(KIND(1.D0)) :: Z
      if (N.le.1) then
        Z=1.D0
      else
        Z=1.D0
        do i=1,N
          Z=Z*i
        enddo
      endif
      FACT=Z
      end FUNCTION FACT

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DOUBLEFACT(N)
C N double factorial
C------------------------------------------
      integer, intent(in) :: N
      integer :: N2
      REAL(KIND(1.D0)) :: Z
      if (N.le.1) then
        Z=1.D0
      else
        if (MOD(N*1.0,2.0).eq.0) then
          N2=N/2
          Z=FACT(N2)*2**N2
c          write(*,*) 'DOUBLEFACT EVEN: ',N,N2
c          write(*,*) 'DOUBLEFACT N2!,Z: ',FACT(N2),Z
        else
          N2=(N-1)/2
          Z=FACT(N)/FACT(N2)/2**N2
c          write(*,*) 'DOUBLEFACT ODD: ',N,N2
c          write(*,*) 'DOUBLEFACT N!,N2!,Z: ',FACT(N),FACT(N2),Z
        endif
      endif
      DOUBLEFACT=Z
      end FUNCTION DOUBLEFACT

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION GAMMAN2(N)
C GAMMA(N/2)
C------------------------------------------
      integer, intent(in) :: N
      REAL(KIND(1.D0)),parameter :: CSQRTPI=1.772453850905516D0
      integer :: N2
      REAL(KIND(1.D0)) :: Z
      Z=0.D0
      if (N.EQ.1) then
        Z=CSQRTPI
      else if (N.EQ.2) then
        Z=1.D0
      else
        if (MOD(N*1.0,2.0).eq.0) then
          N2=(N-2)/2
          Z=FACT(N2)
c          write(*,*) 'GAMMAN2 EVEN: ',N,N2
c          write(*,*) 'GAMMAN2 FACT(N2): ',Z
        else
          N2=(N-1)/2
          Z=CSQRTPI*DOUBLEFACT(N-2)/2**N2
c          write(*,*) 'GAMMAN2 ODD: ',N,N2
c          write(*,*) 'GAMMAN2 DOUBLEFACT(N-2): ',DOUBLEFACT(N-2),Z
        endif
      endif
      GAMMAN2=Z
      end FUNCTION GAMMAN2

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DETJAC(A,MD,N)
C Calculate matrix determinant using Jacobi diagonalization
C------------------------------------------
      integer, intent(in) :: MD,N
      REAL(KIND(1.D0)),intent(in) :: A(MD,MD)
      REAL(KIND(1.D0)) :: DIA(MD),VT(MD,MD),D
      integer :: NROT,I
      call JACOBI(A,N,MD,DIA,VT,NROT)
      D=1.D0
      do i=1,N
        D=D*DIA(i)
      enddo
      DETJAC=D
      end FUNCTION DETJAC

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DETNAP(A,MD,N)
C Calculate matrix determinant using NAPACK library
C preserves A
C------------------------------------------
      integer, intent(in) :: MD,N
      REAL(KIND(1.D0)),intent(in) :: A(MD,MD)
      REAL(KIND(1.D0)) :: AUX(2+N*(N+2)),D,KDET
      integer :: i,j,k
      k=0
      do j=1,N
        do i=1,N
          k=k+1
          AUX(k)=A(i,j)
!          AUX(i,j)=A(i,j)
        enddo
      enddo
      call KFACT(AUX(1),MD,N)
      D=KDET(I,AUX(1),2+N*(N+2))
      DETNAP=D*10.D0**I
      end FUNCTION DETNAP


c------------------------------------------
      integer function INVERT_R4(N,A,NA,B,NB)
C     Inverts matrix A (A is not destroyed)
C returns 1 if not invertible, otherwise 0
C------------------------------------------
      INTEGER,intent(in) ::  N,NA,NB
      real,intent(in) :: A(NA,NA)
      real,intent(out) :: B(NB,NB)
      integer,PARAMETER :: NMAX=256
      real(KIND(1.D0)) :: A1(NMAX,NMAX),WK(2*NMAX)
      INTEGER :: I,J,IRES,KVERTD
      if (N<1) then
            INVERT_R4=0
            return
      endif      
      DO I=1,N
       DO J=1,N
          A1(I,J)=A(I,J)
        enddo
      enddo
      IRES=KVERTD(A1,NMAX,N,WK)
      IF (IRES.EQ.0) then
        DO I=1,N
          DO J=1,N
            B(I,J)=A1(I,J)
          enddo
        enddo
      else
        write(*,*) 'INVERT: cannot invert matrix, dim=',N
      endif
      INVERT_R4 = IRES
      END function INVERT_R4

c------------------------------------------
      integer function INVERT_R8(N,A,NA,B,NB)
C     Inverts matrix A (A is not destroyed)
C returns 1 if not invertible, otherwise 0
C------------------------------------------
      INTEGER,intent(in) ::  N,NA,NB
      real(KIND(1.D0)),intent(in) :: A(NA,NA)
      real(KIND(1.D0)),intent(out) :: B(NB,NB)
      integer,PARAMETER :: NMAX=256
      real(KIND(1.D0)) :: A1(NMAX,NMAX),WK(2*NMAX)
      INTEGER :: I,J,IRES,KVERTD
      if (N<1) then
            INVERT_R8=0
            return
      endif
      DO I=1,N
       DO J=1,N
          A1(I,J)=A(I,J)
        enddo
      enddo
      IRES=KVERTD(A1,NMAX,N,WK)
      IF (IRES.EQ.0) then
        DO I=1,N
          DO J=1,N
            B(I,J)=A1(I,J)
          enddo
        enddo
      else
        write(*,*) 'INVERT: cannot invert matrix, dim=',N
      endif
      INVERT_R8 = IRES
      END function INVERT_R8

c------------------------------------------
      integer function INVERT2(N,A,NA)
C     Inverts matrix A (A is inout)
C returns 1 if not invertible, otherwise 0
C------------------------------------------
      INTEGER,intent(in) ::  N,NA
      real(KIND(1.D0)),intent(inout) :: A(NA,NA)
      integer,PARAMETER :: NMAX=256
      real(KIND(1.D0)) :: A1(NMAX,NMAX),WK(2*NMAX)
      INTEGER :: I,J,IRES,KVERTD
      if (N<1) then
            INVERT2=0
            return
      endif      
      DO I=1,N
       DO J=1,N
          A1(I,J)=A(I,J)
        enddo
      enddo
      IRES=KVERTD(A1,NMAX,N,WK)
      IF (IRES.EQ.0) then
        DO I=1,N
          DO J=1,N
            A(I,J)=A1(I,J)
          enddo
        enddo
      else
        write(*,*) 'INVERT: cannot invert matrix, dim=',N
      endif
      INVERT2 = IRES
      END function INVERT2

C--------------------------------------
      SUBROUTINE BABT(A,B,C,DIM,N)
C calculate square matrix product C=B*A*B^T
C--------------------------------------
      integer, intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(in) :: A(DIM,DIM),B(DIM,DIM)
      REAL(KIND(1.D0)),intent(out) :: C(DIM,DIM)
      INTEGER :: i,j,k,m
      do i=1,N
        do j=1,N
          C(i,j)=0.D0
          do k=1,N
            do m=1,N
              C(i,j)=C(i,j)+B(i,k)*B(j,m)*A(k,m)
            enddo
          enddo
        enddo
      enddo
      end SUBROUTINE BABT

!--------------------------------------
      SUBROUTINE BTAB(A,B,C,N)
! calculate square matrix product C=B^T*A*B
!--------------------------------------
      integer, intent(in) :: N
      REAL(KIND(1.D0)),intent(in) :: A(:,:),B(:,:)
      REAL(KIND(1.D0)),intent(out) :: C(:,:)
      INTEGER :: i,j,k,m
      do i=1,N
        do j=1,N
          C(i,j)=0.D0
          do k=1,N
            do m=1,N
              C(i,j)=C(i,j)+B(k,i)*B(m,j)*A(k,m)
            enddo
          enddo
        enddo
      enddo
      end SUBROUTINE BTAB

!--------------------------------------
      SUBROUTINE BTAC(B,A,C,D,N)
! calculate vector D=B^T*A*C
! B,A ... matrices NxN
! C   ... vector N
!--------------------------------------
      integer, intent(in) :: N
      REAL(KIND(1.D0)),intent(in) :: A(:,:),B(:,:),C(:)
      REAL(KIND(1.D0)),intent(out) :: D(:)
      INTEGER :: i,j,k,m
      do i=1,N
        do j=1,N
          D(i)=0.D0
          do k=1,N
            do m=1,N
              D(i)=D(i)+B(k,i)*A(k,m)*C(m)
            enddo
          enddo
        enddo
      enddo
      end SUBROUTINE BTAC

!------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION XTMX(X,M,NM)
! return X^.M.X
! M is a square matrix, dim>=NM
! X is vector, DIM >= NM
!------------------------------------------------------------
      integer, intent(in) :: NM
      real(KIND(1.D0)), intent(in) :: M(:,:),X(:)
      real(KIND(1.D0)) :: suma,aux
      integer :: i,j
      suma=0.D0
      do j=1,NM
        aux=0.D0
        do i=1,NM
          aux=aux+M(i,j)*X(i)
        enddo
        suma=suma+aux*X(j)
      enddo
      XTMX=suma
      end FUNCTION XTMX

C--------------------------------------
      SUBROUTINE M3XV3(IT,M,B,C)
C matrix.vector, 3-dimensional
C use transposed matrix if IT=-1
C just copy C=B, if IT=0
C--------------------------------------
      INTEGER,intent(in) :: IT
      REAL(KIND(1.D0)),intent(in) :: M(3,3),B(3)
      REAL(KIND(1.D0)),intent(out) :: C(3)
      INTEGER :: I,J
      select case(IT)
      case(0)
        do i=1,3
          C(i)=B(i)
        enddo
      case(-1)
        C=(/(0.D0,i=1,3)/)
        DO J=1,3
          DO I=1,3
            C(J)=C(J)+M(I,J)*B(I)
          enddo
        enddo
      case default
        C=(/(0.D0,i=1,3)/)
        DO I=1,3
          DO J=1,3
            C(J)=C(J)+M(J,I)*B(I)
          enddo
        enddo
      end select
      END SUBROUTINE M3XV3

C     ------------------------------
      SUBROUTINE MXV(IT,N,NP,A,B,C)
C     ------------------------------
      INTEGER*4 IT,N,NP,I,J
      REAL*8 A(NP,NP),B(NP),C(NP)
      DO 10 J=1,N
         C(J)=0.
         IF (IT.GT.0) THEN
           DO 20 I=1,N
20         C(J)=C(J)+A(J,I)*B(I)
         ELSE
           DO 30 I=1,N
30         C(J)=C(J)+A(I,J)*B(I)
         ENDIF
10    CONTINUE
      RETURN
      END SUBROUTINE MXV

C-----------------------------------------
      SUBROUTINE MXM(IT,N,ND,A,B,C)
C multiply two square matrices
C-----------------------------------------
      integer, intent(in) :: IT,N,ND
      REAL(KIND(1.D0)), intent(in) :: A(ND,ND),B(ND,ND)
      REAL(KIND(1.D0)), intent(out) :: C(ND,ND)
      INTEGER :: I,J,K

      DO J=1,N
        DO K=1,N
          C(J,K)=0.
          IF (IT.GT.0) THEN
            DO I=1,N
              C(J,K)=C(J,K)+A(J,I)*B(I,K)
            enddo
          ELSE
            DO I=1,N
              C(J,K)=C(J,K)+A(I,J)*B(I,K)
            enddo
          ENDIF
        enddo
      enddo
      END SUBROUTINE MXM


C     ------------------------------
      SUBROUTINE M3XM3(IT,A,B,C)
C     ------------------------------
      INTEGER*4 IT,I,J,K
      REAL*8 A(3,3),B(3,3),C(3,3)
      DO 10 J=1,3
      DO 10 K=1,3
         C(J,K)=0.
         IF (IT.GT.0) THEN
           DO 20 I=1,3
20         C(J,K)=C(J,K)+A(J,I)*B(I,K)
         ELSE
           DO 30 I=1,3
30         C(J,K)=C(J,K)+A(I,J)*B(I,K)
         ENDIF
10    CONTINUE
      RETURN
      END SUBROUTINE M3XM3

C     ------------------------------
      SUBROUTINE V3AV3(IT,A,B,C)
C     ------------------------------
      INTEGER*4 IT,I
      REAL*8 A(3),B(3),C(3)

      DO 10 I=1,3
10      C(I)=A(I)+IT*B(I)
      RETURN
      END SUBROUTINE V3AV3

C     ------------------------------
      REAL*8 FUNCTION ABSV3(A)
C     ------------------------------
      REAL*8 A(3)
      ABSV3=SQRT(V3XV3(A,A))
      END FUNCTION ABSV3

C     ------------------------------
      REAL*8   FUNCTION V3XV3(A,B)
C     ------------------------------
      INTEGER*4 I
      REAL*8 A(3),B(3),Z
      Z=0
      DO 10 I=1,3
10      Z=Z+A(I)*B(I)
      V3XV3=Z
      RETURN
      END FUNCTION V3XV3


!----------------------------------
! cross product of 3-dim real vectors
!----------------------------------
      SUBROUTINE V3cV3(A,B,C)
      real(KIND(1.D0)) :: A(3),B(3),C(3)
      C(1)=A(2)*B(3)-A(3)*B(2)
      C(2)=A(3)*B(1)-A(1)*B(3)
      C(3)=A(1)*B(2)-A(2)*B(1)
      end SUBROUTINE V3cV3

C----------------------------------------------------
      SUBROUTINE GENROT(IAX,PHI,ROT,LD)
C generate rotation matrix ROT(3,3) by angle PHI
C IAX determines the axis index and direction
C LD is the dimension of ROT
C rotates axes 1..3, other dimensions are invariant
C------------------------------------------------------
      implicit none
      integer,intent(in) :: IAX,LD
      real(kind(1.D0)),intent(in) :: PHI
      real(kind(1.D0)),intent(out) :: ROT(LD,LD)
      real(kind(1.D0)) :: SI,CO
      integer :: I,J
      call UNIMAT(ROT,LD,LD)
      SI=SIN(PHI)
      CO=COS(PHI)
      if (IAX.eq.2) SI=-SI
      DO I=1,3
        if (I.NE.IAX) then
          DO J=1,3
            if (J.NE.IAX) then
              IF(I.EQ.J) THEN
                ROT(I,J)=CO
              ELSE
                IF (I.GT.J) THEN
                  ROT(I,J)=SI
                ELSE
                  ROT(I,J)=-SI
                endif
              endif
            ENDIF
          enddo
        ENDIF
      ENDDO
      END SUBROUTINE GENROT

c**************************************************************
c
      REAL*8 FUNCTION DETERM1(B,MD,N,A)
C     COMPUTES THE DETERMINANT OF THE MATRIX B
      INTEGER*4 I,J,K,N1,N,I1,J1,K1,J2,K2,MD
      real*8 A(MD,N),B(MD,N),X
      DO 55 I=1,N
      DO 55 J=1,N
   55 A(I,J)=B(I,J)
      N1=N-1
      DETERM1=1.
      DO 1 I=1,N1
      J1=I
      K1=I
      DO 10 J2=I,N
      DO 10 K2=I,N
      IF(ABS(A(J1,K1)).GE.ABS(A(J2,K2)))GO TO 10
      J1=J2
      K1=K2
10    CONTINUE
      IF(ABS(A(J1,K1)).GT.1.E-30)GO TO 11
      DETERM1=0.
      RETURN
11    CONTINUE
      IF(J1.EQ.I)GO TO 12
      DO 5 K=I,N
      X=A(I,K)
      A(I,K)=A(J1,K)
5     A(J1,K)=-X
12    IF(K1.EQ.I)GO TO 13
      DO 6 J=1,N
      X=A(J,I)
      A(J,I)=A(J,K1)
6     A(J,K1)=-X
13    I1=I+1
      DO 30 J=I1,N
      IF(A(J,I).EQ.0.)GO TO 30
      X=A(J,I)/A(I,I)
      DO 7 K=I,N
7     A(J,K)=A(J,K)-X*A(I,K)
30    CONTINUE
1     DETERM1=DETERM1*A(I,I)
      DETERM1=DETERM1*A(N,N)
      RETURN
      END FUNCTION DETERM1
C
C

c**************************************************************
c
      REAL*8 FUNCTION DETERM(B,N,A)
C     COMPUTES THE DETERMINANT OF THE MATRIX B
      INTEGER*4 I,J,K,N1,N,I1,J1,K1,J2,K2
      real*8 A(N,N),B(N,N),X
      DO 55 I=1,N
      DO 55 J=1,N
   55 A(I,J)=B(I,J)
      N1=N-1
      DETERM=1.
      DO 1 I=1,N1
      J1=I
      K1=I
      DO 10 J2=I,N
      DO 10 K2=I,N
      IF(ABS(A(J1,K1)).GE.ABS(A(J2,K2)))GO TO 10
      J1=J2
      K1=K2
10    CONTINUE
      IF(ABS(A(J1,K1)).GT.1.E-30)GO TO 11
      DETERM=0.
      RETURN
11    CONTINUE
      IF(J1.EQ.I)GO TO 12
      DO 5 K=I,N
      X=A(I,K)
      A(I,K)=A(J1,K)
5     A(J1,K)=-X
12    IF(K1.EQ.I)GO TO 13
      DO 6 J=1,N
      X=A(J,I)
      A(J,I)=A(J,K1)
6     A(J,K1)=-X
13    I1=I+1
      DO 30 J=I1,N
      IF(A(J,I).EQ.0.)GO TO 30
      X=A(J,I)/A(I,I)
      DO 7 K=I,N
7     A(J,K)=A(J,K)-X*A(I,K)
30    CONTINUE
1     DETERM=DETERM*A(I,I)
      DETERM=DETERM*A(N,N)
      RETURN
      END FUNCTION DETERM
C
C
C
        SUBROUTINE DIAG(A,ADA,B)
C***********************************************************************
C   diagonalizes real*4 matrix A(4,4), B(4,4) is corresponding rotation matrix
C***********************************************************************
        INTEGER*4 I,J,K,L,N,ND,KK,KI,KJ,JR,KL,NDK,NDI,NDJ,JTES,NDN,II
        INTEGER*4 IJ,JI,JJ,JK,IK,ITES,JRMAX(16)
        REAL*4 A(16),ADA(16),B(16),ARMAX(16)
        REAL*4 E,Y,X,T,TY,TSQ,C,S,CSQ,AMAX,AII,AJJ,AIJ
        DATA N,ND,E/4,4,1.E-24/
        NDN=ND*N
        DO 1 K=1,NDN
        ADA(K)=A(K)
        B(K)=0.
    1        CONTINUE
        DO 2 K=1,N
        KK=K*(ND+1)-ND
        ARMAX(K)=0.
        B(KK)=1.
        DO 3 L=K,N
        IF(L-K)4,3,4
    4        KL=K+ND*(L-1)
        Y=ABS(ADA(KL))
        IF(ARMAX(K)-Y)5,3,3
    5        ARMAX(K)=Y
        JRMAX(K)=L
    3        CONTINUE
    2        CONTINUE
   11        AMAX=0.
        DO 6 K=1,N
        Y=ABS(ARMAX(K))
        IF(AMAX-Y)7,6,6
    7        AMAX=Y
        I=K
    6        CONTINUE
        J=JRMAX(I)
        IF(E-AMAX)8,9,9
    8        NDI=ND*(I-1)
        NDJ=ND*(J-1)
        II=I+NDI
        JJ=J+NDJ
        IJ=I+NDJ
        JI=J+NDI
        AII=ADA(II)
        AJJ=ADA(JJ)
        AIJ=ADA(IJ)
        Y=2.*AIJ
        X=AII-AJJ
        T=SIGN(1.E0,X)*Y/(ABS(X)+SQRT(X**2+Y**2))
        TSQ=T**2
        C=1./SQRT(ABS(1.+TSQ))
        TY=T*Y
        S=T*C
        CSQ=C**2
        ADA(II)=CSQ*(AII+TY+AJJ*TSQ)
        ADA(JJ)=CSQ*(AJJ-TY+AII*TSQ)
        ADA(IJ)=0.
        ADA(JI)=0.
        DO 10 K=1,N
        JTES=(K-I)*(K-J)
        NDK=ND*(K-1)
        KI=K+NDI
        KJ=K+NDJ
        IF(JTES)13,12,13
   13        JK=J+NDK
        IK=I+NDK
        ADA(KI)=C*ADA(IK)+S*ADA(JK)
        ADA(KJ)=-S*ADA(IK)+C*ADA(JK)
        ADA(JK)=ADA(KJ)
        ADA(IK)=ADA(KI)
   12        X=B(KI)
        B(KI)=C*X+S*B(KJ)
        B(KJ)=-S*X+C*B(KJ)
   10        CONTINUE
        ARMAX(I)=0.
        DO 14 K=1,N
        IF(K-I)15,14,15
   15        IK=I+ND*(K-1)
        Y=ABS(ADA(IK))
        IF(ARMAX(I)-Y)16,14,14
   16        ARMAX(I)=Y
        JRMAX(I)=K
   14        CONTINUE
        ARMAX(J)=0.
        DO 17 K=1,N
        IF(K-J)18,17,18
   18        JK=J+ND*(K-1)
        Y=ABS(ADA(JK))
        IF(ARMAX(J)-Y)19,17,17
   19        ARMAX(J)=Y
        JRMAX(J)=K
   17        CONTINUE
        DO 20 K=1,N
        ITES=(K-I)*(K-J)
        KI=K+NDI
        KJ=K+NDJ
        IF(ITES)21,20,21
   21        X=ABS(ADA(KI))
        Y=ABS(ADA(KJ))
        JR=J
        IF(X-Y)22,22,23
   23        Y=X
        JR=I
   22        IF(ARMAX(K)-Y)24,20,20
   24        ARMAX(K)=Y
        JRMAX(K)=JR

   20        CONTINUE
        GOTO 11
9        CONTINUE
        RETURN
        END SUBROUTINE DIAG

C-----------------------------------------------------------------
      SUBROUTINE LINFIT(X,Y,N,DX,DY,DZ,ND,AMPL,BACK,CHISQ)
C     Linear fit of the function (X,Y,N) to the data (DX,DY,DZ,ND)
C in REAL*4 !!
C-----------------------------------------------------------------
      INTEGER*4 I,K,N,ND,KK
      REAL*4 AMPL,BACK,CHISQ
      REAL*4 X(N),Y(N),DX(ND),DY(ND),DZ(ND)
      REAL*4 C1,C2,C3,C4,C5,C6,Z,YY,W,DDX,ZMIN

      C1=0
      C2=0
      C3=0
      C4=0
      C5=0
      C6=0
      KK=0
      zmin=0
      DO i=1,N
         if(abs(Y(i)).gt.zmin) zmin=abs(Y(i))
      enddo
      zmin=abs(zmin/10.)
      DO I=1,ND
        DDX=X(2)-X(1)
        Z=(DX(I)-X(1))/DDX
        IF(Z.GE.0) THEN
           K=INT(Z)+1
        ELSE
           K=INT(Z)
        ENDIF
        IF((K.GT.0).AND.(K.LT.N)) THEN
           KK=KK+1
           YY=Y(K)+(Y(K+1)-Y(K))*(Z+1-K)  ! linear interpolation
C           SIG2=DZ(I)**2
C           IF(SIG2.EQ.0) SIG2=1.
C           W=1/SIG2

            IF (abs(YY).LE.zmin) THEN
               W=SQRT(ZMIN)
            ELSE
               W=SQRT(abs(YY))   ! weighted by SQRT(Y) (for RESTRAX only)
            ENDIF
           C1=C1+DY(I)*W
           C2=C2+YY*W
           C3=C3+W
           C4=C4+DY(I)*YY*W
           C5=C5+YY*YY*W
           C6=C6+DY(I)*DY(I)*W
        ENDIF
      END DO
      IF(KK.GT.0) THEN
         IF((C2*C2-C3*C5).EQ.0) THEN
            AMPL=0
            BACK=0
            CHISQ=C6
         ELSE
            AMPL= (C1*C2-C3*C4)/(C2*C2-C3*C5)
            BACK=(C4*C2-C5*C1)/(C2*C2-C3*C5)
            CHISQ=(AMPL**2)*C5+(BACK**2)*C3+C6+2*AMPL*BACK*C2-2*AMPL*C4-
     1         2*BACK*C1
            CHISQ=CHISQ/KK/C3
         ENDIF
      ELSE
         WRITE(*,*) 'Cannot fit data ! '
         AMPL=0.
         BACK=0.
      ENDIF
      RETURN
      END SUBROUTINE LINFIT

C------------------------------------------------
      REAL*4 FUNCTION LINTERP4(X,Y,N,Z)
C linear interpolation on equidistant data
C in REAL*4 !!
C------------------------------------------------
      INTEGER*4 I0,N
      REAL*4 Z0,X(N),Y(N),Z,DX

      IF (Z.LE.X(1)) THEN
         LINTERP4=Y(1)
         RETURN
      ELSE IF (Z.GE.X(N)) THEN
         LINTERP4=Y(N)
         RETURN
      ELSE
         DX=(X(N)-X(1))/(N-1)
         Z0=(Z-X(1))/DX
         I0=INT(Z0)+1
         LINTERP4=Y(I0)+(Y(I0+1)-Y(I0))*(Z0-I0+1)
         RETURN
      ENDIF
      END FUNCTION LINTERP4

C------------------------------------------------
      REAL*8 FUNCTION LINTERP8(X,Y,N,Z)
C linear interpolation on equidistant data
C in REAL*8!!
C------------------------------------------------
      INTEGER*4 I0,N
      REAL*8 Z0,X(N),Y(N),Z,DX

      IF (Z.LE.X(1)) THEN
         LINTERP8=Y(1)
         RETURN
      ELSE IF (Z.GE.X(N)) THEN
         LINTERP8=Y(N)
         RETURN
      ELSE
         DX=(X(N)-X(1))/(N-1)
         Z0=(Z-X(1))/DX
         I0=INT(Z0)+1
         LINTERP8=Y(I0)+(Y(I0+1)-Y(I0))*(Z0-I0+1)
         RETURN
      ENDIF
      END FUNCTION LINTERP8

C------------------------------------------------
      REAL*8 FUNCTION LINTERP(Y,N,X0,DX,Z)
C linear interpolation on equidistant data
C x=X0+(i-1)*DX, i=1..N
C extrapolation=0
C------------------------------------------------
      INTEGER,intent(in) :: N
      REAL(KIND(1.D0)),intent(in) :: Y(N),X0,DX,Z
      REAL(KIND(1.D0)) :: Z0
      INTEGER :: I0
      IF ((Z.LE.X0).OR.(Z.GE.X0+(N-1)*DX)) THEN
         LINTERP=0.D0
      ELSE
         Z0=(Z-X0)/DX
         I0=INT(Z0)+1
         LINTERP=Y(I0)+(Y(I0+1)-Y(I0))*(Z0-I0+1)
      ENDIF
      END FUNCTION LINTERP

C------------------------------------------------
      REAL*8 FUNCTION LOGINTERP(Y,N,X0,DX,Z)
C linear interpolation on logaritmic scale
C x=X0*exp((i-1)*DX), i=1..N
C extrapolation=0
C------------------------------------------------
      INTEGER,intent(in) :: N
      REAL(KIND(1.D0)),intent(in) :: Y(N),X0,DX,Z
      REAL(KIND(1.D0)) :: Z0
      INTEGER :: I0
      IF ((Z.LE.X0).or.(X0.LE.0.D0)) then
        LOGINTERP=0.D0
      ELSE
        Z0=LOG(Z/X0)/DX
        I0=INT(Z0)+1
        if (I0.GE.N) then
          LOGINTERP=0.D0
        else
          LOGINTERP=Y(I0)+(Y(I0+1)-Y(I0))*(Z0-I0+1)
        endif
      endif
      END FUNCTION LOGINTERP

!-----------------------------------------------------
      REAL*8 FUNCTION LINTERP2D(Y,NI,NJ,ND,ZI,ZJ,DI,DJ,Z)
! linear interpolation on equidistant data in 2D
! x=ZI+(i-1)*DI, i=1..NI
! y=ZJ+(j-1)*DJ, j=1..NJ
! extrapolation=0
!-----------------------------------------------------
      INTEGER*4 NI,NJ,ND
      REAL*8 Y(ND,NJ),ZI,ZJ,DI,DJ,Z(2)
      REAL*8 JJ,YY1,YY2
      INTEGER*4 J0
      IF (Z(2).LE.ZJ.OR.Z(2).GE.ZJ+(NJ-1)*DJ) THEN
         LINTERP2D=0.D0
         RETURN
      ELSE
         JJ=(Z(2)-ZJ)/DJ
         J0=INT(JJ)+1
         YY1=LINTERP(Y(1,J0),NI,ZI,DI,Z(1))
         YY2=LINTERP(Y(1,J0+1),NI,ZI,DI,Z(1))
         LINTERP2D=YY1+(YY2-YY1)*(JJ-J0+1)
         RETURN
      ENDIF
      END FUNCTION LINTERP2D

!-----------------------------------------------------
      REAL*8 FUNCTION LOGLINTERP2D(Y,NI,NJ,ND,ZI,ZJ,DI,DJ,Z)
! 2D interpolation in a table, which has logarithmic scale in the 1st indexes
! and linear in the 2nd.
! x=ZI*exp((i-1)*DI), i=1..NI
! y=ZJ+(j-1)*DJ, j=1..NJ
! extrapolation=0
!-----------------------------------------------------
      INTEGER*4 NI,NJ,ND
      REAL*8 Y(ND,NJ),ZI,ZJ,DI,DJ,Z(2)
      REAL*8 JJ,YY1,YY2
      INTEGER*4 J0
      IF (Z(2).LE.ZJ.OR.Z(2).GE.ZJ+(NJ-1)*DJ) THEN
        LOGLINTERP2D=0.D0
      ELSE
        JJ=(Z(2)-ZJ)/DJ
        J0=INT(JJ)+1
        YY1=LOGINTERP(Y(1,J0),NI,ZI,DI,Z(1))
        YY2=LOGINTERP(Y(1,J0+1),NI,ZI,DI,Z(1))
        LOGLINTERP2D=YY1+(YY2-YY1)*(JJ-J0+1)
      ENDIF
      END FUNCTION LOGLINTERP2D

C------------------------------------------------
      REAL*4 FUNCTION QINTERP4(X,Y,N,Z)
C quadratic interpolation (X monotonous)
C in REAL*4 !!
C------------------------------------------------
      INTEGER*4 I0,N,I
      REAL*4 X(N),Y(N),Z,A(3,3),B(3),C(3)

      IF (Z.LE.X(1)) THEN
         QINTERP4=Y(1)
         RETURN
      ELSE IF (Z.GE.X(N)) THEN
         QINTERP4=Y(N)
         RETURN
      ELSE
        I0=1
        DO WHILE(X(I0).LT.Z)
          I0=I0+1
        ENDDO
        IF(I0.GT.N-1) I0=N-1
        IF(I0.LT.2) I0=2
        IF((X(I0).EQ.X(I0-1)).OR.(X(I0).EQ.X(I0+1))) THEN
           IF(X(I0-1).EQ.X(I0+1)) THEN
              QINTERP4=Y(I0)
           ELSE
              C(1)=(Y(I0+1)-Y(I0-1))/(X(I0+1)-X(I0-1))
              QINTERP4=Y(I0-1)+C(1)*(Z-X(I0-1))
           ENDIF
           RETURN
        ENDIF
        DO I=1,3
          A(I,1)=X(I+I0-2)**2
          A(I,2)=X(I+I0-2)
          A(I,3)=1
          B(I)=Y(I+I0-2)
        ENDDO
        CALL GAUSS(A,B,C,3)
        QINTERP4=C(1)*Z**2+C(2)*Z+C(3)
      ENDIF
      END FUNCTION QINTERP4

!------------------------------------------------
      REAL(kind(1.D0)) FUNCTION QINTERP8(X,Y,N,Z)
! quadratic interpolation (X monotonous)
! in DOUBLE !!
!------------------------------------------------
      INTEGER,intent(in) :: N
      REAL(kind(1.D0)),intent(in) :: X(N),Y(N),Z
      INTEGER*4 I
      REAL(kind(1.D0)) :: aa,bb,cc,d(3),RES,SS
      if (N<3) then
        RES=LINTERP8(X,Y,N,Z)
      ELSE IF (Z.LE.X(1)) THEN
        RES=Y(1)
      ELSE IF (Z.GE.X(N)) THEN
        RES=Y(N)
      ELSE
        i=1
        SS=SIGN(1.D0,X(2)-X(1))
        DO WHILE((SS*X(i).LE.SS*Z).and.(i<N))
          i=i+1
        ENDDO
        i=i-1
        IF (i.GT.N-2) i=N-2
        call quadinterp3(X(i:i+2),Y(i:i+2),d)
        aa=d(1)
        bb=d(2)-2*d(1)*x(i)
        cc=d(3)-d(2)*x(i)+d(1)*x(i)**2
        RES=aa*z**2+bb*z+cc
      endif
      QINTERP8=RES
      END FUNCTION QINTERP8

C------------------------------------------------
      SUBROUTINE GAUSS(A,B,C,N)
C// solve linear equations by Gauss elimination
C in REAL*4 !!
C------------------------------------------------
      INTEGER, PARAMETER :: DMAX=64
      INTEGER*4 N,I,J,K
      REAL*4 A(N,N),C(N),B(N),AUX(DMAX,DMAX+1),V(DMAX+1),M,sum

      DO I=1,N
      DO J=1,N
        AUX(I,J)=A(I,J)
      ENDDO
      ENDDO

      DO I=1,N
        AUX(I,N+1)=B(I)
      ENDDO

      DO K=1,N-1
        I=K
        DO WHILE(AUX(I,K).EQ.0)
          I=I+1
          IF (I.GT.N) GOTO 10
        ENDDO
        DO J=1,N+1
          V(j)=AUX(i,j)
          AUX(i,j)=AUX(k,j)
          AUX(k,j)=V(j)
        ENDDO
        DO I=k+1,N
          M=AUX(I,K)/AUX(k,k)
          do j=k,N+1
            AUX(i,j)=AUX(i,j)-m*AUX(k,j)
          enddo
        enddo
      enddo
      if (AUX(N,N).EQ.0) GOTO 10
      C(N)=AUX(N,N+1)/AUX(N,N)
      do k=1,n-1
        sum=0
        do j=n-k+1,n
          sum=sum+AUX(n-k,j)*C(j)
        enddo
        C(n-k)=(AUX(n-k,n+1)-sum)/AUX(n-k,n-k)
      enddo
      RETURN
10    write(*,*) 'Matrix not invertible'
      DO I=1,N
        C(I)=0
      ENDDO
      return
      end SUBROUTINE GAUSS

!-----------------------------------------------------------------------------------
      REAL*8 FUNCTION ROT3(I,J,K,ALFA)
! Return element(J,K) of a rotation matrix with rot. axis I
! ALPHA is the rotation angle in [rad]
! Assume right-handed system, i.e. ROT(1,2,3,A)=+sin(A), ROT(1,3,2,A)=-sin(A)
! If I<0, use oposite rotation direction
! NOTE: ROT3 does not rotate vectors - it transforms vector coordinates into rotated basis
! to rotate vectors, use inverse of ROT3 (use I<0)
!------------------------------------------------------------------------------------
      INTEGER,intent(in) :: I,J,K
      REAL(KIND(1.D0)),intent(in) :: ALFA
      integer :: IA,ISG

      ISG=SIGN(1,I)
      IA=ABS(I)
      IF (IA.EQ.2) ISG=-ISG
      IF((IA.EQ.J).OR.(IA.EQ.K)) THEN
        IF(J.EQ.K) THEN
          ROT3=1.D0
        ELSE
          ROT3=0.D0
        ENDIF
      ELSE IF (J.EQ.K) THEN
        ROT3=COS(alfa)
      ELSE
        IF(J.GT.K) THEN
           ROT3=-ISG*sin(alfa)
        ELSE
           ROT3=ISG*sin(alfa)
        ENDIF
      ENDIF
      END FUNCTION ROT3

!--------------------------------------------------------
      SUBROUTINE MK_ROT3(I,ALFA,RT)
! generate 3-dim rotation matrix with rot. axis I
! see subroutine ROT3 for rules
!--------------------------------------------------------
      integer,intent(in) :: I
      REAL(KIND(1.D0)),intent(in) :: ALFA
      REAL(KIND(1.D0)),intent(out) :: RT(3,3)
      INTEGER*4 J,K

      DO J=1,3
      DO K=1,3
         RT(J,K)=ROT3(I,J,K,ALFA)
      ENDDO
      ENDDO
      END SUBROUTINE MK_ROT3

!--------------------------------------------------------
      SUBROUTINE Givens(a,b,co,si,r)
! Givens rotation
! a,b are the diagonal and off-diagonal matrix elements
! return cos, sin and r
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: a,b
      REAL(KIND(1.D0)),intent(out) :: co,si,r
      REAL(KIND(1.D0)) :: t,u
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      if (abs(b).lt.EPS) then
        co=sign(1.D0,a)
        si=0.D0
        r=abs(a)
      else if (abs(a).lt.EPS) then
        co=0.D0
        si=sign(1.D0,b)
        r=abs(b)
      else
        if (abs(b).gt.abs(a)) then
          t=a/b
          u=sign(1.D0,b)*sqrt(1+t**2)
          si=1/u
          co=si*t
          r=b*u
        else
          t=b/a
          u=sign(1.D0,a)*sqrt(1+t**2)
          co=1/u
          si=co*t
          r=a*u
        endif
      endif
      END SUBROUTINE Givens


!--------------------------------------------------------
      SUBROUTINE getEulerYZX(R,e)
! Return (y,z,x) Euler angles, given 3x3 rotation matrix
! transforming coordinates from lab to local frame
! Suitable for A3 (omega)  and  GL, GU (lower and upper craddles)
! of ILL TAS instruments
! Inverse to getRotYZX
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: R(3,3)
      REAL(KIND(1.D0)),intent(out) :: e(3)
      REAL(KIND(1.D0)) :: e1(3)
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      REAL(KIND(1.D0)) :: si(3),co(3),c22,RZ(3,3),R2(3,3),ss2
      integer :: i
      c22=SQRT(R(1,1)**2 + R(1,3)**2)
      if (c22.lt.1.D-6) then
        ss2=sign(1.D0,R(1,2))
        call MK_ROT3(3,-ss2*PI/2,RZ)
        call M3xM3(1,RZ,R,R2)
        e=(/atan2(-R2(1,3),R2(1,1)),ss2*PI/2,0.D0/)
      else
        call Givens(R(1,1),-R(1,3),co(1),si(1),co(2))
        si(3)=si(1)*R(2,1)+co(1)*R(2,3)
        co(3)=si(1)*R(3,1)+co(1)*R(3,3)
        si(2)=R(1,2)
        do i=1,3
          e(i)=atan2(si(i),co(i))
        enddo
    ! minimize e(2),e(3)
        e1(2)=-e(2)
        e1(1)=e(1)+PI
        e1(3)=e(3)+PI
        if (e1(1).gt.PI) e1(1)=e1(1)-2*PI
        if (e1(3).gt.PI) e1(3)=e1(3)-2*PI
        if ((abs(e(2))+abs(e(3))).gt.(abs(e1(2))+abs(e1(3)))) e=e1
      endif
    ! round zeros
      do i=1,3
        if (abs(e(i)).le.EPS) e(i)=0.D0
      enddo
      end SUBROUTINE getEulerYZX

!--------------------------------------------------------
      SUBROUTINE getEulerYXY(R,e)
! Return (y,x,y) Euler angles, given 3x3 rotation matrix
! transforming coordinates from lab to local frame
! Standard gonio angles for all SIMRES components
! Inverse to getRotYXY
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: R(3,3)
      REAL(KIND(1.D0)),intent(out) :: e(3)
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      REAL(KIND(1.D0)) :: si(3),co(3),s22,sc2
      integer :: i
      s22=SQRT(R(1,1)**2 + R(1,3)**2)
      if (s22.lt.1.D-6) then
        sc2=sign(1.D0,R(2,2))
        e=(/atan2(-R(1,3),R(1,1)),(1.D0-sc2)*PI/2,0.D0/)
      else
        call Givens(R(2,3),R(2,1),co(1),si(1),si(2))
        si(3)=co(1)*R(3,1)-si(1)*R(3,3)
        co(3)=co(1)*R(1,1)-si(1)*R(1,3)
        co(2)=R(2,2)
        do i=1,3
          e(i)=atan2(si(i),co(i))
        enddo
      endif
    ! round zeros
      do i=1,3
        if (abs(e(i)).le.EPS) e(i)=0.D0
      enddo
      end SUBROUTINE getEulerYXY


!--------------------------------------------------------
      SUBROUTINE getEulerYZY(R,e)
! Return (y,z,y) Euler angles, given 3x3 rotation matrix
! transforming coordinates from lab to local frame
! Suitable for positioning the Kf axis (AX sample angles in SIMRES and RESTRAX)
! Inverse to getRotYZY
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: R(3,3)
      REAL(KIND(1.D0)),intent(out) :: e(3)
      REAL(KIND(1.D0)) :: R1(3,3),R2(3,3),AUX(3,3)
      call MK_ROT3(2,PI/2,R1)
      call M3XM3(1,R,R1,AUX)
      call M3XM3(-1,R1,AUX,R2)
      call getEulerYXY(R2,e)
      end SUBROUTINE getEulerYZY

!--------------------------------------------------------
      SUBROUTINE getRotYXY(e,R)
! Return Rotation matrix from (y,x,y) Euler angles
! Inverse to getEulerYXY
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: e(3)
      REAL(KIND(1.D0)),intent(out) :: R(3,3)
      REAL(KIND(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),AUX(3,3)
        call MK_ROT3(2,e(1),R1)
        call MK_ROT3(1,e(2),R2)
        call MK_ROT3(2,e(3),R3)
        call M3XM3(1,R2,R1,AUX)
        call M3XM3(1,R3,AUX,R)
        call roundMatrix(R)
      end SUBROUTINE getRotYXY

!--------------------------------------------------------
      SUBROUTINE getRotYZX(e,R)
! Return Rotation matrix from (y,x,y) Euler angles
! Inverse to getEulerYXY
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: e(3)
      REAL(KIND(1.D0)),intent(out) :: R(3,3)
      REAL(KIND(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),AUX(3,3)
        call MK_ROT3(2,e(1),R1)
        call MK_ROT3(3,e(2),R2)
        call MK_ROT3(1,e(3),R3)
        call M3XM3(1,R2,R1,AUX)
        call M3XM3(1,R3,AUX,R)
        call roundMatrix(R)
      end SUBROUTINE getRotYZX

!--------------------------------------------------------
      SUBROUTINE getRotYZY(e,R)
! Return Rotation matrix from (y,z,y) Euler angles
! Inverse to getEulerYZY
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: e(3)
      REAL(KIND(1.D0)),intent(out) :: R(3,3)
      REAL(KIND(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),AUX(3,3)
        call MK_ROT3(2,e(1),R1)
        call MK_ROT3(3,e(2),R2)
        call MK_ROT3(2,e(3),R3)
        call M3XM3(1,R2,R1,AUX)
        call M3XM3(1,R3,AUX,R)
        call roundMatrix(R)
      end SUBROUTINE getRotYZY

!--------------------------------------------------------
      SUBROUTINE roundMatrix(R)
! Round mateix elements to remove numerical uncertainity
!--------------------------------------------------------
      REAL(KIND(1.D0)),intent(inout) :: R(3,3)
      REAL(KIND(1.D0)),parameter :: EPS=1.D-14
      REAL(KIND(1.D0)) :: a
      integer :: i,j
        do i=1,3
        do j=1,3
          a=R(i,j)
          if (abs(a)<EPS) a=0.D0
          if (abs(a-1.D0)<EPS) a=1.D0
          if (abs(a+1.D0)<EPS) a=-1.D0
          R(i,j)=a
        enddo
        enddo
      end SUBROUTINE roundMatrix



c--------------------------------------------------------
      SUBROUTINE UNIMAT(U,DIM,N)
c create unit square matrix
C N ..  matrix rank
C DIM .. leading dimension
c--------------------------------------------------------
      integer,intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(out) :: U(DIM,DIM)
      INTEGER :: I,J
      DO I=1,N
        DO J=I+1,N
          U(I,J)=0.D0
          U(J,I)=0.D0
        ENDDO
        U(I,I)=1.D0
      ENDDO
      END SUBROUTINE UNIMAT


c--------------------------------------------------------
      SUBROUTINE getRotZX(Z,X,R)
! Return rotation matrix to the coordinates defined by
! z'//Z, y'//ZxX
c--------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: Z(3),X(3)
      REAL(KIND(1.D0)),intent(out) :: R(3,3)
      REAL(KIND(1.D0)) :: U(3,3),W
      integer :: i,j
      W=SQRT(Z(1)**2+Z(2)**2+Z(3)**2)
      U(1:3,3)=Z(1:3)/W
      call V3cV3(Z,X,U(1,2))
      W=SQRT(U(1,2)**2+U(2,2)**2+U(3,2)**2)
      U(1:3,2)=U(1:3,2)/W
      call V3cV3(U(1,2),U(1,3),U(1,1))
      do i=1,3
        do j=1,3
          R(i,j)=U(j,i)
        enddo
      enddo
      end SUBROUTINE getRotZX


C------------------------------------------
      SUBROUTINE JACOBI(AA,N,NP,D,V,NROT)
C  Diagonalize matrix AA,
C returns the diagonal, D
C and the transformation matrix, V
C  D=V^T*AA*V
C  (Numerical Recipes)
C------------------------------------------
      integer,intent(in) :: N,NP
      integer,intent(out) :: NROT
      real(KIND(1.D0)),intent(in) :: AA(NP,NP)
      real(KIND(1.D0)),intent(out) :: D(NP),V(NP,NP)
      real(KIND(1.D0)) :: A(NP,NP)

    !  INTEGER,PARAMETER :: NMAX=100
      real(KIND(1.D0)) :: B(NP),Z(NP),G
      integer :: IQ,IP,I
      real(KIND(1.D0)) :: SM,TRESH
      Z=0.D0
      call UNIMAT(V,NP,N)
      DO IP=1,N
        DO IQ=1,N
          A(IP,IQ)=AA(IP,IQ)
        enddo
      enddo
      DO IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
      enddo
      NROT=0
      DO I=1,50
        SM=0.D0
        DO IP=1,N-1
          DO IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
          enddo
        enddo
      !  write(*,*) 'JACOBI iter=',I,SM,N,NP
        IF(SM.LE.1.D-30) GOTO 99
        IF(I.LT.4) THEN
          TRESH=0.2D0*SM/N**2
        ELSE
          TRESH=0.D0
        ENDIF
        DO IP=1,N-1
          DO IQ=IP+1,N
            G=100.D0*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ)))) THEN
              A(IP,IQ)=0.D0
            ELSE IF (ABS(A(IP,IQ)).GT.TRESH) THEN
              call JACOBI_ROT(IP,IQ,A,V,D,Z,N,NP)
              NROT=NROT+1
            ENDIF
          enddo
        enddo
        DO IP=1,N
          B(IP)=B(IP)+Z(IP)
          D(IP)=B(IP)
          Z(IP)=0.D0
        enddo
      enddo
      WRITE(*,*) 'JACOBI : 50 ITERATION SHOULD NEVER HAPPEN'
      RETURN
99    CONTINUE
      END SUBROUTINE JACOBI


C------------------------------------------------------------
      SUBROUTINE JACOBI_ROT(IP,IQ,A,V,D,Z,N,NP)
C  Make one Jacobi rotation for IP,IQ element
C------------------------------------------------------------
      integer,intent(in) :: IP,IQ,N,NP
      real(KIND(1.D0)),intent(inout) :: A(NP,NP),V(NP,NP),D(NP),Z(NP)
      real(KIND(1.D0)) :: G,H,S,C,T,THETA,TAU
      integer :: J
      !  write(*,*) 'JACOBI_ROT ',IP,IQ
        G=100.D0*ABS(A(IP,IQ))
        H=D(IQ)-D(IP)
        IF (ABS(H)+G.EQ.ABS(H)) THEN
          T=A(IP,IQ)/H
        ELSE
          THETA=0.5D0*H/A(IP,IQ)
          T=1.D0/(ABS(THETA)+SQRT(1.+THETA**2))
          IF(THETA.LT.0.D0) T=-T
        ENDIF
        S=SQRT(1.D0+T**2)
        C=1.D0/S
        S=C*T
        TAU=S/(1.D0+C)
        H=T*A(IP,IQ)
        Z(IP)=Z(IP)-H
        Z(IQ)=Z(IQ)+H
        D(IP)=D(IP)-H
        D(IQ)=D(IQ)+H
        A(IP,IQ)=0.D0
        DO J=1,IP-1
          G=A(J,IP)
          H=A(J,IQ)
          A(J,IP)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
        enddo
        DO J=IP+1,IQ-1
          G=A(IP,J)
          H=A(J,IQ)
          A(IP,J)=G-S*(H+G*TAU)
          A(J,IQ)=H+S*(G-H*TAU)
        enddo
        DO J=IQ+1,N
          G=A(IP,J)
          H=A(IQ,J)
          A(IP,J)=G-S*(H+G*TAU)
          A(IQ,J)=H+S*(G-H*TAU)
        enddo
        DO J=1,N
          G=V(J,IP)
          H=V(J,IQ)
          V(J,IP)=G-S*(H+G*TAU)
          V(J,IQ)=H+S*(G-H*TAU)
        enddo
      end SUBROUTINE JACOBI_ROT


C------------------------------------------------------------
      SUBROUTINE EXCHANGE_INDEX(A,I1,I2)
C  exchange cols and rows of a 4x4 matrix for indexes I1 and I2
C------------------------------------------------------------
      INTEGER,intent(in) :: I1,I2
      REAL(kind(1.D0)) :: A(4,4)
      REAL(kind(1.D0)) :: C(4),R(4)
      integer :: I
      if ((I1.ne.I2).and.(I1.gt.0).and.(I2.gt.0).and.(I1.le.4).and.(I2.le.4)) then
        C(1:4)=A(1:4,I1)
        A(1:4,I1)=A(1:4,I2)
        A(1:4,I2)=C(1:4)
        do i=1,4
          R(i)=A(I1,i)
          A(I1,i)=A(I2,i)
          A(I2,i)=R(i)
        enddo
      endif
      end SUBROUTINE EXCHANGE_INDEX

C------------------------------------------------------------
      SUBROUTINE EXCHANGE_INDEX_S(A,I1,I2)
C Single precission version of  EXCHANGE_INDEX
C------------------------------------------------------------
      INTEGER,intent(in) :: I1,I2
      REAL :: A(4,4)
      REAL :: C(4),R(4)
      integer :: I
      if ((I1.ne.I2).and.(I1.gt.0).and.(I2.gt.0).and.(I1.le.4).and.(I2.le.4)) then
        C(1:4)=A(1:4,I1)
        A(1:4,I1)=A(1:4,I2)
        A(1:4,I2)=C(1:4)
        do i=1,4
          R(i)=A(I1,i)
          A(I1,i)=A(I2,i)
          A(I2,i)=R(i)
        enddo
      endif
      end SUBROUTINE EXCHANGE_INDEX_S

!--------------------------------------------------------------------
      SUBROUTINE EllipsoidView4(A,IX,IY,PROJ,SECT)
! Works with single precission REAL !!!
! Given A(4,4) matrix describing a 4D ellipsoid:
! Return projection and section 2x2 matrices for the plane (IX,IY)
!--------------------------------------------------------------------
      integer, intent(in) :: IX,IY
      REAL, intent(in) :: A(4,4)
      REAL, intent(out) ::PROJ(2,2),SECT(2,2)
      REAL :: M(4,4,4)
      integer :: I,J,K,IDX(2)

  ! projection indices
      IDX=(/IX,IY/)

      M(1:4,1:4,4)=A(1:4,1:4)
      call EXCHANGE_INDEX_S(M(1,1,4),IX,1)
      call EXCHANGE_INDEX_S(M(1,1,4),IY,2)
  ! projection
      DO 30 I=4,3,-1
      DO 30 J=1,I
      DO 30 K=1,I
        M(J,K,I-1)=M(J,K,I)-M(J,I,I)*M(K,I,I)/M(I,I,I)
30    CONTINUE
      PROJ(1:2,1:2)=M(1:2,1:2,2)
  ! section
      do I=1,2
      DO J=1,2
        SECT(I,J)=A(IDX(I),IDX(J))
      enddo
      enddo
      end SUBROUTINE EllipsoidView4

!--------------------------------------------------------------------
      integer function SGN2INT(SIG)
! converts sign (-1,+1) to ordinal value (0,1)
!--------------------------------------------------------------------
      integer, intent(in) :: SIG
        if (SIG<0) then
          SGN2INT=0
        else
          SGN2INT=1
        endif
      end function SGN2INT

!--------------------------------------------------------------------
      integer function INT2SGN(NUM)
! converts ordinal value (0,1) to sign (-1,+1)
!--------------------------------------------------------------------
      integer, intent(in) :: NUM
        if (NUM.le.0.D0) then
          INT2SGN=-1
        else
          INT2SGN=1
        endif
      end function INT2SGN

!--------------------------------------------------------------------
      integer function LOG2INT(LVAL)
! converts logical (false,true) to ordinal value (0,1)
!--------------------------------------------------------------------
      logical, intent(in) :: LVAL
        if (LVAL) then
          LOG2INT=1
        else
          LOG2INT=0
        endif
      end function LOG2INT

!--------------------------------------------------------------------
      SUBROUTINE HEAPSORT(ARG,NARG,DIR,ORD)
! Make an index to given array using Heapsort algorithm.
! NOTE: the array is not sorted, only the index is created.
! DIR>=0 for ascending, DIR<0 for descending order
!--------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: ARG(*)
      integer, intent(in) :: DIR,NARG
      integer, intent(out) :: ORD(*)
      integer :: iroot,i,j,L,IR
      do i=1,NARG
        ORD(i)=i
      enddo
      if (NARG.le.1) return
      L=NARG/2+1
      IR=NARG
10    if (L > 1) then
        L=L-1
        iroot=ORD(L)
      else
        iroot=ORD(IR)
        ORD(IR)=ORD(1)
        IR=IR-1
        if (IR.eq.1) then
          ORD(1)=iroot
          goto 30
        end if
      end if
      I=L
      J=L+L
20    if (J.le.IR) then
        if (J < IR) then
          if (ARG(ORD(J)) < ARG(ORD(J+1)))  J=J+1
        endif
        if (ARG(iroot) < ARG(ORD(J))) then
          ORD(I)=ORD(J)
          I=J
          J=J+J
        else
          J=IR+1
        end if
        goto 20
      end if
      ORD(I)=iroot
      goto 10
! revert order for descending sort
30    if (DIR.lt.0) then
        i=1
        j=NARG
        do while (i.lt.j)
          L=ORD(i)
          ORD(i)=ORD(j)
          ORD(j)=L
          i=i+1
          j=j-1
        enddo
      endif
      END SUBROUTINE HEAPSORT

!--------------------------------------------------------------------
      SUBROUTINE BERNOULLI(NMAX,B)
! return a series of Bernoulli numbers from 0 to NMAX.
! NMAX<=30
!--------------------------------------------------------------------
      integer,intent(in) :: NMAX
      real(kind(1.D0)),intent(out) :: B(0:)
      integer :: i
      real(kind(1.D0)),parameter :: BER2(0:15)= (/
     & 1.D0,
     & 1.D0/6,
     & -1.D0/30,
     & 1.D0/42,
     & -1.D0/30,
     & 5.D0/66,
     & -6.91D2/2730,
     & 7.D0/6,
     & -3.617D3/510,
     & 4.3867D4/798,
     & -1.74611D5/330,
     & 8.54513D5/138,
     & -2.36364091D8/2730,
     & 8.553103D6/6,
     & -2.3749461029D10/870,
     & 8.615841276005D12/14322
     &/)

      do i=0,min(NMAX,30)
        if (i == 1) then
          B(i)=-0.5D0
        else if (MOD(i,2) == 0) then
          B(i)=BER2(i/2)
        else
          B(i)=0.D0
        endif
      enddo
      end SUBROUTINE BERNOULLI

!--------------------------------------------------------------------
      real(kind(1.D0)) function INTEG1D(y,dx,n)
! integrate given array using Simpson/3 method
! data must be equidistant with step dx
! n must be odd integer, otherwise the data are extrapolated by 0
!--------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: y(*),dx
      integer,intent(in) :: n
      real(kind(1.D0)) :: Z
      INTEGER,parameter :: P(2)=(/2,4/)
      integer :: i,m
      Z=0.D0
      ! get the nearest even integer m<=n
      m=n
      if (mod(n,2).eq.1) m=m-1
      if (n.gt.2) then
        Z=y(1)+P(2)*y(2)
        do i=2,m/2
          Z=Z+P(1)*y(i*2-1)+P(2)*y(2*i)
        enddo
        if (n.gt.m) Z=Z+y(n)
      endif
      INTEG1D=Z*dx/3.D0
      end function INTEG1D

!--------------------------------------------------------------------
      real(kind(1.D0)) function INTEG1DLOG(y,x0,dx,n)
! integrate given array using Simpson/3 method
! data must be on equidistant logarithmic x-scale with log-step dx
! x(i)=x0*exp((i-1)*dx), i=1..n
! n must be odd integer, otherwise the data are extrapolated by 0
!--------------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: y(*),x0,dx
      integer,intent(in) :: n
      real(kind(1.D0)) :: Z,XX(2)
      INTEGER,parameter :: P(2)=(/2,4/)
      integer :: i,m
      Z=0.D0
      ! get the nearest even integer m<=n
      m=n
      if (mod(n,2).eq.1) m=m-1
      if (n.gt.2) then
        XX(1)=x0
        XX(2)=x0*exp(dx)
        Z=XX(1)*y(1)+P(2)*XX(2)*y(2)
        do i=2,m/2
          XX(1)=x0*exp((i*2-2)*dx)
          XX(2)=x0*exp((i*2-1)*dx)
          Z=Z+P(1)*XX(1)*y(i*2-1)+P(2)*XX(2)*y(2*i)
        enddo
        if (n.gt.m) Z=Z+y(n)*x0*exp((n-1)*dx)
      endif
      INTEG1DLOG=Z*dx/3.D0
      end function INTEG1DLOG

!--------------------------------------------------------------------
      SUBROUTINE spline(x,y,n,yp1,ypn,y2)
! from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
! rewritten to f90 double precission
!Given arrays x(1:n) and y(1:n) containing a tabulated function, i.e., yi = f(xi), with
!x1 < x2 < .. . < xN, and given values yp1 and ypn for the first derivative of the interpolating
!function at points 1 and n, respectively, this routine returns an array y2(1:n) of
!length n which contains the second derivatives of the interpolating function at the tabulated
!points xi. If yp1 and/or ypn are equal to 1 × 1030 or larger, the routine is signaled to set
!the corresponding boundary condition for a natural spline, with zero second derivative on
!that boundary.
!Parameter: NMAX is the largest anticipated value of n.
!--------------------------------------------------------------------
      INTEGER,intent(in) :: n
      REAL(kind(1.D0)),intent(in) :: yp1,ypn,x(n),y(n)
      REAL(kind(1.D0)),intent(out) :: y2(n)
      INTEGER :: i,k
      integer,parameter :: NMAX=500
      REAL(kind(1.D0)) :: p,qn,sig,un,u(NMAX)
      if (yp1>0.99D30) then
        y2(1)=0.D0
        u(1)=0.D0
      else
        y2(1)=-0.5D0
        u(1)=(3.D0/(x(2)-x(1)))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif
      do i=2,n-1
        sig=(x(i)-x(i-1))/(x(i+1)-x(i-1))
        p=sig*y2(i-1)+2.
        y2(i)=(sig-1.D0)/p
        u(i)=(6.D0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
      enddo
      if (ypn>0.99D30) then
        qn=0.D0
        un=0.D0
      else
        qn=0.5D0
        un=(3.D0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif
      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.D0)
      do k=n-1,1,-1
        y2(k)=y2(k)*y2(k+1)+u(k)
      enddo
      END SUBROUTINE spline

!--------------------------------------------------------------------
      SUBROUTINE splint(xa,ya,y2a,n,x,y)
! from NUMERICAL RECIPES IN FORTRAN 77: THE ART OF SCIENTIFIC COMPUTING (ISBN 0-521-43064-X)
! rewritten to f90 double precission
!Given the arrays xa(1:n) and ya(1:n) of length n, which tabulate a function (with the
!xai’s in order), and given the array y2a(1:n), which is the output from spline above,
!and given a value of x, this routine returns a cubic-spline interpolated value y.
!--------------------------------------------------------------------
      INTEGER,intent(in) :: n
      REAL(kind(1.D0)),intent(in) :: x,xa(n),ya(n),y2a(n)
      REAL(kind(1.D0)),intent(out) :: y
      INTEGER :: k,khi,klo
      REAL(kind(1.D0)) :: a,b,h
      klo=1
      khi=n
1     if (khi-klo.gt.1) then
        k=(khi+klo)/2
        if(xa(k).gt.x)then
          khi=k
        else
          klo=k
        endif
        goto 1
      endif ! klo and khi now bracket the input value of x.
      h=xa(khi)-xa(klo)
      if (h.eq.0.D0) then
        y=ya(khi)
        return ! ’bad xa input in splint’ The xa’s must be distinct.
      endif
      a=(xa(khi)-x)/h !Cubic spline polynomial is now evaluated.
      b=(x-xa(klo))/h
      y=a*ya(klo)+b*ya(khi)+((a**3-a)*y2a(klo)+(b**3-b)*y2a(khi))*(h**2)/6.D0
      END SUBROUTINE splint

!--------------------------------------------------------------------
      SUBROUTINE quadinterp3(x,y,d)
! quadratic interpolation of 3 points
! with respect to x(1)
! given x,y, return derivatives d so that y(x)=d(1)*(x-x(1))**2 + d(2)*(x-x(1)) + d(3)
!--------------------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: x(3),y(3)
      REAL(kind(1.D0)),intent(out) :: d(3)
      REAL(kind(1.D0)) :: h1,h2,dh
      h1=x(2)-x(1)
      h2=x(3)-x(1)
      dh=h2-h1
      d(1)=y(1)/(h1*h2) - y(2)/(h1*dh) + y(3)/(h2*dh)
      d(2)=-y(1)*(h1+h2)/(h1*h2) + y(2)*h2/(h1*dh) - y(3)*h1/(h2*dh)
      d(3)=y(1)

      end SUBROUTINE quadinterp3


!--------------------------------------------------------------------
      REAL(kind(1.D0)) function SERF(x)
! Error function from Numerical Recipes.
! erf(x) = 2/sqrt(pi) INT{exp(-x^2)}dx from 0 to x
!--------------------------------------------------------------------
      REAL(kind(1.D0)), intent(in) :: x
      REAL(kind(1.D0)) :: t, z, arg
      integer :: i
      real(kind(1.D0)),parameter :: P(1:10)= (/
     & 0.17087277D0,-0.82215223D0,1.48851587D0,-1.13520398D0,0.27886807D0,
     &  -0.18628806D0,0.09678418D0,0.37409196D0,1.00002368D0,-1.26551223D0/)
      z = abs(x)
      t = 1.D0 / ( 1.D0 + 0.5D0 * z )
      arg=P(1)*t+P(2)
      do i=3,10
        arg = P(i) + t*arg
      enddo
      arg=arg-z*z
      if (x.lt.0.0) then
        SERF=t*exp(arg)-1.D0
      else
        SERF=1.D0 - t*exp(arg)
      endif
      end function SERF


      end module SIMATH

