C//////////////////////////////////////////////////////////////////////
C////
C////  R E S T R A X   4.1
C////
C////  Subroutines for simple matrix operations:
C////
C////
C//////////////////////////////////////////////////////////////////////

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION FACT(N)
C N factorial
C------------------------------------------
      implicit none
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
      implicit none
      integer, intent(in) :: N
      integer :: N2
      REAL(KIND(1.D0)) :: Z,FACT
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
      implicit none
      integer, intent(in) :: N
      REAL(KIND(1.D0)),parameter :: CSQRTPI=1.772453850905516D0
      integer :: N2
      REAL(KIND(1.D0)) :: Z,FACT,DOUBLEFACT
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
      REAL(KIND(1.D0)) FUNCTION DETJAC(A,DIM,N)
C Calculate matrix determinant using Jacobi diagonalization
C------------------------------------------
      implicit none

      integer, intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(in) :: A(DIM,DIM)
      REAL(KIND(1.D0)) :: AUX(DIM,DIM),DIA(DIM),VT(DIM,DIM),D
      integer :: NROT,I

      call JACOBI(A,AUX,N,DIM,DIA,VT,NROT)
      D=1.D0
      do i=1,N
        D=D*DIA(i)
      enddo
      DETJAC=D
      end

C------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DETNAP(A,DIM,N)
C Calculate matrix determinant using NAPACK library
C preserves A
C------------------------------------------
      implicit none

      integer, intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(in) :: A(DIM,DIM)
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
      call KFACT(AUX(1),DIM,N)
      D=KDET(I,AUX(1),2+N*(N+2))
      DETNAP=D*10.D0**I
      end


c------------------------------------------
      SUBROUTINE INVERT(N,A,NA,B,NB)
C     Inverts matrix A (A is not destroyed)
C returns 1 if not invertible, otherwise 0
C------------------------------------------
      IMPLICIT NONE
      INTEGER*4 N,NA,NB,NMAX
      PARAMETER(NMAX=64)
      REAL*8 A(NA,NA),B(NB,NB),A1(NMAX,NMAX),WK(2*NMAX)
      INTEGER*4 I,J,IRES,KVERTD

      DO 5 I=1,N
      DO 5 J=1,N
5        A1(I,J)=A(I,J)
!       write(*,*) 'res_mat.f, INVERT, N,NA,NB: ',N,NA,NB
      IRES=KVERTD(A1,NMAX,N,WK)
      IF (IRES.EQ.0) then
        DO 10 I=1,N
        DO 10 J=1,N
10            B(I,J)=A1(I,J)
      else
        write(*,*) 'INVERT: cannot invert matrix, dim=',N
        READ(*,*)
      endif
      END


C--------------------------------------
      SUBROUTINE BABT(A,B,C,DIM,N)
C calculate square matrix product C=B*A*B^T
C--------------------------------------
      IMPLICIT NONE
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
      end

C--------------------------------------
      SUBROUTINE BTAB(A,B,C,DIM,N)
C calculate square matrix product C=B^T*A*B
C--------------------------------------
      IMPLICIT NONE
      integer, intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(in) :: A(DIM,DIM),B(DIM,DIM)
      REAL(KIND(1.D0)),intent(out) :: C(DIM,DIM)
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
      end





C--------------------------------------
      SUBROUTINE M3XV3(IT,M,B,C)
C matrix.vector, 3-dimensional
C use transposed matrix if IT=-1
C just copy C=B, if IT=0
C--------------------------------------
      IMPLICIT NONE
      INTEGER*4 IT,I,J
      REAL*8 M(3,3),B(3),C(3)
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
      END

C     ------------------------------
      SUBROUTINE MXV(IT,N,NP,A,B,C)
C     ------------------------------
      IMPLICIT NONE
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
      END

C-----------------------------------------
      SUBROUTINE MXM(IT,N,ND,A,B,C)
C multiply two square matrices
C-----------------------------------------
      IMPLICIT NONE
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
      END


C     ------------------------------
      SUBROUTINE M3XM3(IT,A,B,C)
C     ------------------------------
      IMPLICIT NONE
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
      END

C     ------------------------------
      SUBROUTINE V3AV3(IT,A,B,C)
C     ------------------------------
      IMPLICIT NONE
      INTEGER*4 IT,I
      REAL*8 A(3),B(3),C(3)

      DO 10 I=1,3
10      C(I)=A(I)+IT*B(I)
      RETURN
      END

C     ------------------------------
      REAL*8 FUNCTION ABSV3(A)
C     ------------------------------
      IMPLICIT NONE
      REAL*8 A(3),V3XV3
      ABSV3=SQRT(V3XV3(A,A))
      RETURN
      END



C     ------------------------------
      REAL*8 FUNCTION V3XV3(A,B)
C     ------------------------------
      IMPLICIT NONE
      INTEGER*4 I
      REAL*8 A(3),B(3),Z
      Z=0
      DO 10 I=1,3
10      Z=Z+A(I)*B(I)
      V3XV3=Z
      RETURN
      END


!----------------------------------
! cross product of 3-dim real vectors
!----------------------------------
      SUBROUTINE V3cV3(A,B,C)
      IMPLICIT NONE
      real(KIND(1.D0)) :: A(3),B(3),C(3)
      integer :: i,j,k
      integer,parameter :: ix(6)=(/2,3,3,1,1,2/)
      do i=1,3
        j=ix(2*i-1)
        k=ix(2*i)
        C(i)=A(j)*B(k)-A(k)*B(j)
      enddo
      end SUBROUTINE V3cV3

C     ------------------------------
      SUBROUTINE GENROT(IAX,PHI,AUX)
C     ------------------------------
      IMPLICIT NONE
      INTEGER*4 I,J,IAX
      REAL*8 PHI,CO,SI,AUX(3,3)

      SI=SIN(PHI)
      CO=SQRT(1-SI**2)
      DO 20 I=1,3
      DO 20 J=1,3
      IF(I.EQ.J) THEN
         IF (I.EQ.IAX) THEN
            AUX(I,J)=1.
         ELSE
            AUX(I,J)=CO
         ENDIF
      ELSE
         IF((I.EQ.IAX).OR.(J.EQ.IAX)) THEN
           AUX(I,J)=0.
         ELSE IF (I.GT.J) THEN
           AUX(I,J)=SI
         ELSE
           AUX(I,J)=-SI
         ENDIF
      ENDIF
20    CONTINUE
      RETURN
      END

c**************************************************************
c
      REAL*8 FUNCTION DETERM1(B,MD,N,A)
C     COMPUTES THE DETERMINANT OF THE MATRIX B
      implicit none
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
      END
C
C

c**************************************************************
c
      REAL*8 FUNCTION DETERM(B,N,A)
C     COMPUTES THE DETERMINANT OF THE MATRIX B
      implicit none
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
      END
C
C
C
        SUBROUTINE DIAG(A,ADA,B)
C***********************************************************************
C   diagonalizes real*4 matrix A(4,4), B(4,4) is corresponding rotation matrix
C***********************************************************************
        IMPLICIT NONE
        INTEGER*4 I,J,K,L,N,ND,KK,KI,KJ,JR,KL,NDK,NDI,NDJ,JTES,NDN,II
        INTEGER*4 IJ,JI,JJ,JK,IK,ITES
        REAL*4 A(16),ADA(16),B(16),ARMAX(16),JRMAX(16)
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
        END

C-----------------------------------------------------------------
      SUBROUTINE LINFIT(X,Y,N,DX,DY,DZ,ND,AMPL,BACK,CHISQ)
C     Linear fit of the function (X,Y,N) to the data (DX,DY,DZ,ND)
C in REAL*4 !!
C-----------------------------------------------------------------
      IMPLICIT NONE
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
      END

C------------------------------------------------
      REAL*4 FUNCTION LINTERP4(X,Y,N,Z)
C linear interpolation on equidistant data
C in REAL*4 !!
C------------------------------------------------
      IMPLICIT NONE

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
      END

C------------------------------------------------
      REAL*8 FUNCTION LINTERP8(X,Y,N,Z)
C linear interpolation on equidistant data
C in REAL*8!!
C------------------------------------------------
      IMPLICIT NONE

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
      END

C------------------------------------------------
      REAL*8 FUNCTION LINTERP(Y,N,X0,DX,Z)
C linear interpolation on equidistant data
C x=X0+(i-1)*DX, i=1..N
C extrapolation=0
C------------------------------------------------
      IMPLICIT NONE
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
      IMPLICIT NONE
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

C-----------------------------------------------------
      REAL*8 FUNCTION LINTERP2D(Y,NI,NJ,ND,ZI,ZJ,DI,DJ,Z)
C linear interpolation on equidistant data in 2D
C x=ZI+(i-1)*DI, i=1..NI
C y=ZJ+(j-1)*DJ, j=1..NJ
C extrapolation=0
C-----------------------------------------------------
      IMPLICIT NONE

      INTEGER*4 NI,NJ,ND
      REAL*8 Y(ND,NJ),ZI,ZJ,DI,DJ,Z(2)
      REAL*8 JJ,YY1,YY2
      INTEGER*4 J0
      REAL*8 LINTERP

c      write(*,*)  '2D: ',ZJ,Z(2),ZJ+(NJ-1)*DJ

      IF (Z(2).LE.ZJ.OR.Z(2).GE.ZJ+(NJ-1)*DJ) THEN
         LINTERP2D=0.D0
         RETURN
      ELSE
         JJ=(Z(2)-ZJ)/DJ
         J0=INT(JJ)+1
         YY1=LINTERP(Y(1,J0),NI,ZI,DI,Z(1))
         YY2=LINTERP(Y(1,J0+1),NI,ZI,DI,Z(1))
         LINTERP2D=YY1+(YY2-YY1)*(JJ-J0+1)
c       write(*,*) '2D: ', JJ,J0,YY1,YY2
         RETURN
      ENDIF

      END
C------------------------------------------------
      REAL*4 FUNCTION QINTERP4(X,Y,N,Z)
C quadratic interpolation (X monotonous)
C in REAL*4 !!
C------------------------------------------------
      IMPLICIT NONE

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
      END




C------------------------------------------------
      SUBROUTINE GAUSS(A,B,C,N)
C// solve linear equations by Gauss elimination
C in REAL*4 !!
C------------------------------------------------
      IMPLICIT NONE
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
      pause
      return

      end

c------------------------------------------------------------------
      REAL*8 FUNCTION ROT3(I,J,K,ALFA)
c Return element(J,K) of a rotation matrix with rot. axis I
C ALPHA is the rotation angle in [rad]
C Assume right-handed system, i.e. ROT(1,2,3,A)=+sin(A), ROT(1,3,2,A)=-sin(A)
C If I<0, use oposite rotation direction
c------------------------------------------------------------------
      IMPLICIT NONE

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
      END

c--------------------------------------------------------
      SUBROUTINE MK_ROT3(I,ALFA,RT)
c generate 3-dim rotation matrix with rot. axis I
C see subroutine ROT3 for rules
c--------------------------------------------------------
      IMPLICIT NONE
      integer,intent(in) :: I
      REAL(KIND(1.D0)),intent(in) :: ALFA
      REAL(KIND(1.D0)),intent(out) :: RT(3,3)
      INTEGER*4 J,K
      REAL(KIND(1.D0)) :: ROT3

      DO J=1,3
      DO K=1,3
         RT(J,K)=ROT3(I,J,K,ALFA)
      ENDDO
      ENDDO
      END

c--------------------------------------------------------
      SUBROUTINE Givens(a,b,co,si,r)
c Givens rotation
c a,b are the diagonal and off-diagonal matrix elements
c return cos, sin and r
c--------------------------------------------------------
      IMPLICIT NONE
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

c--------------------------------------------------------
      SUBROUTINE getEulerYXY(R,e)
! Return (y,x,y) Euler angles
! given 3x3 rotation matrix
c--------------------------------------------------------
      use CONSTANTS
      IMPLICIT NONE
      REAL(KIND(1.D0)),intent(in) :: R(3,3)
      REAL(KIND(1.D0)),intent(out) :: e(3)
      REAL(KIND(1.D0)) :: e1(3)
      REAL(KIND(1.D0)),parameter :: EPS=1.D-10
      REAL(KIND(1.D0)) :: si(3),co(3)
      integer :: i
      if (abs(abs(R(2,2))-1.D0).lt.EPS) then
        e=(/atan2(-R(1,3),R(1,1)),(1.D0-sign(1.D0,R(2,2)))*PI/2,0.D0/)
      else
        call Givens(R(2,3),R(2,1),co(1),si(1),si(2))
        si(3)=co(1)*R(3,1)-si(1)*R(3,3)
        co(3)=co(1)*R(1,1)-si(1)*R(1,3)
        co(2)=R(2,2)
        do i=1,3
          e(i)=atan2(si(i),co(i))
        enddo
    ! minimize e(1),e(3)
        e1(2)=-e(2)
        e1(1)=e(1)+PI
        e1(3)=e(3)+PI
        if (e1(1).gt.PI) e1(1)=e1(1)-2*PI
        if (e1(3).gt.PI) e1(3)=e1(3)-2*PI
        if ((abs(e(1))+abs(e(3))).gt.(abs(e1(1))+abs(e1(3)))) e=e1
      endif
    ! round zeros
      do i=1,3
        if (abs(e(i)).le.EPS) e(i)=0.D0
      enddo
      end SUBROUTINE getEulerYXY

c--------------------------------------------------------
      SUBROUTINE getRotYXY(e,R)
! Return Rotation matrix from (y,x,y) Euler angles
! Inverse to getEulerYXY
c--------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND(1.D0)),intent(in) :: e(3)
      REAL(KIND(1.D0)),intent(out) :: R(3,3)
      REAL(KIND(1.D0)) :: R1(3,3),R2(3,3),R3(3,3),AUX(3,3)
        call MK_ROT3(2,e(1),R1)
        call MK_ROT3(1,e(2),R2)
        call MK_ROT3(2,e(3),R3)
        call M3XM3(1,R2,R1,AUX)
        call M3XM3(1,R3,AUX,R)
      end SUBROUTINE getRotYXY

c--------------------------------------------------------
      SUBROUTINE UNIMAT(U,DIM,N)
c create unit square matrix
C N ..  matrix rank
C DIM .. leading dimension
c--------------------------------------------------------
      IMPLICIT NONE
      integer,intent(in) :: DIM,N
      REAL(KIND(1.D0)),intent(out) :: U(DIM,DIM)
      INTEGER*4 I,J
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
      IMPLICIT NONE
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
      SUBROUTINE JACOBI(AA,A,N,NP,D,V,NROT)
C  Diagonalize matrix AA,
C returns the diagonal, D
C and the transformation matrix, V
C  D=V^T*AA*V
C  (Numerical Recipes)
C------------------------------------------
      implicit none
      integer,intent(in) :: N,NP
      integer,intent(out) :: NROT
      real(KIND(1.D0)),intent(in) :: AA(NP,NP)
      real(KIND(1.D0)),intent(out) :: A(NP,NP),D(NP),V(NP,NP)
      INTEGER,PARAMETER :: NMAX=100
      real(KIND(1.D0)) :: B(NMAX),Z(NMAX)
      integer :: IQ,IP,I,J
      real(KIND(1.D0)) :: SM,G,S,C,T,TAU,H,TRESH,THETA
      DO IP=1,N
        DO IQ=1,N
          V(IP,IQ)=0.D0
          A(IP,IQ)=AA(IP,IQ)
        enddo
        V(IP,IP)=1.D0
      enddo
      DO IP=1,N
        B(IP)=A(IP,IP)
        D(IP)=B(IP)
        Z(IP)=0.D0
      enddo
      NROT=0
      DO I=1,50
        SM=0.D0
        DO IP=1,N-1
          DO IQ=IP+1,N
            SM=SM+ABS(A(IP,IQ))
          enddo
        enddo
        IF(SM.EQ.0) GOTO 99
        IF(I.LT.4) THEN
          TRESH=0.2*SM/N**2
        ELSE
          TRESH=0.D0
        ENDIF
        DO IP=1,N-1
          DO IQ=IP+1,N
            G=100.D0*ABS(A(IP,IQ))
            IF((I.GT.4).AND.(ABS(D(IP))+G.EQ.ABS(D(IP))).AND.(ABS(D(IQ))+G.EQ.ABS(D(IQ)))) THEN
              A(IP,IQ)=0.D0
            ELSE IF (ABS(A(IP,IQ)).GT.TRESH) THEN
              H=D(IQ)-D(IP)
              IF (ABS(H)+G.EQ.ABS(H)) THEN
                T=A(IP,IQ)/H
              ELSE
                THETA=0.5*H/A(IP,IQ)
                T=1./(ABS(THETA)+SQRT(1.+THETA**2))
                IF(THETA.LT.0.D0) T=-T
              ENDIF
              C=1.D0/SQRT(1.D0+T**2)
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
      END


C------------------------------------------------------------
      SUBROUTINE EXCHANGE_INDEX(A,I1,I2)
C  exchange cols and rows of a 4x4 matrix for indexes I1 and I2
C------------------------------------------------------------
      IMPLICIT NONE
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
      IMPLICIT NONE
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


C--------------------------------------------------------------------
      SUBROUTINE EllipsoidView4(A,IX,IY,PROJ,SECT)
c Works with single precission REAL !!!
C Given A(4,4) matrix describing a 4D ellipsoid:
C Return projection and section 2x2 matrices for the plane (IX,IY)
C--------------------------------------------------------------------
      IMPLICIT NONE
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






