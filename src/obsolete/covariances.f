!//////////////////////////////////////////////////////////////////////
!////  $Id: covariances.f,v 1.1 2007/10/08 08:35:37 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.1 $
!////     $Date: 2007/10/08 08:35:37 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module with data for variance sampling optimization
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COVARIANCES
      implicit none

      integer, parameter :: MD=16
      TYPE TSTAT; SEQUENCE
        REAL(KIND(1.D0)) :: SUM2(MD,MD),C(MD,MD)
        REAL(KIND(1.D0)) :: SUM1(MD),M(MD)
        REAL(KIND(1.D0)) :: SUMN,P
        INTEGER :: ND,NC
      END TYPE TSTAT
      TYPE(TSTAT) :: CV(2)

      CONTAINS

C---------------------------------------------
      SUBROUTINE STAT_INI(ID,N)
C initialize covariance matrix data for dim=N
C---------------------------------------------
      integer, intent(in) :: ID,N
      INTEGER :: I,J
      CV(ID).ND=N
      CV(ID).SUMN=0
      CV(ID).NC=0
      DO I=1,CV(ID).ND
        CV(ID).SUM1(I)=0
        DO J=1,CV(ID).ND
          CV(ID).SUM2(I,J)=0
        enddo
      enddo
      END SUBROUTINE

C--------------------------------------------------------------------
      SUBROUTINE STAT_INP(ID,N,X,P)
C accumulates covariantes matrix of a vector X with weight P
C--------------------------------------------------------------------
      integer, intent(in) :: ID,N
      REAL(KIND(1.D0)), intent(in) :: X(N),P
      INTEGER :: I,J

      CV(ID).SUMN=CV(ID).SUMN+P

c      write(*,*) 'STAT_INP ',CV(ID).NC,X(1:N)

      CV(ID).NC=CV(ID).NC+1
      DO I=1,N
        CV(ID).SUM1(I)=CV(ID).SUM1(I)+X(I)*P
        DO J=1,N
          CV(ID).SUM2(I,J)=CV(ID).SUM2(I,J)+X(I)*X(J)*P
        enddo
      enddo
      END SUBROUTINE

C--------------------------------------
      SUBROUTINE STAT_EVAL(ID)
C evaluate covariance matrix and means
C--------------------------------------
      integer, intent(in) :: ID
      INTEGER :: I,J

      if (CV(ID).NC.LE.10) then
        call MSG_ERROR('STAT_EVAL','not enough events for covariance estimation',0,1)
        return
      endif
c      write(*,*) 'STAT_EVAL 1 ',ID
! can't evaluate if NC<10
      IF ((CV(ID).NC.LT.10).OR.(CV(ID).SUMN.LE.0)) RETURN
! normalize
      CV(ID).P=CV(ID).SUMN/CV(ID).NC
      DO I=1,CV(ID).ND
        CV(ID).M(I)=CV(ID).SUM1(I)/CV(ID).SUMN
        DO J=1,CV(ID).ND
          CV(ID).C(I,J)=CV(ID).SUM2(I,J)/CV(ID).SUMN
        enddo
      enddo
c      write(*,*) 'STAT_EVAL 2 ',ID
! subtract mean
      DO I=1,CV(ID).ND
        DO J=1,CV(ID).ND
          CV(ID).C(I,J)=CV(ID).C(I,J)-CV(ID).M(I)*CV(ID).M(J)
        enddo
      enddo
      END SUBROUTINE

C--------------------------------------------------------------------
      SUBROUTINE STAT_DIAG(ID,DIM,U,D)
C diagonalize covariantes matrix so that U(k,j)*CV(k,m)*U(m,j) = D(i)*delta(i,j)
C RETURN U(M,M) and D(M)
C--------------------------------------------------------------------
      integer, intent(in) :: ID,DIM
      REAL(KIND(1.D0)), intent(out) :: U(DIM,DIM),D(DIM)
      REAL(KIND(1.D0)) :: AUX(MD,MD),VT(MD,MD),DIA(MD)
      INTEGER :: I,J,N

      if (CV(ID).NC.LE.10) then
        call MSG_ERROR('STAT_DIAG','not enough events for covariance estimation',0,1)
        return
      endif

      N=CV(ID).ND
      call STAT_EVAL(ID)
c      write(*,*) 'STAT_DIAG 1 ',ID
      CALL JACOBI(CV(ID).C(1,1),AUX(1,1),N,MD,DIA(1),VT,J)
c      write(*,*) 'STAT_DIAG 2 ',ID
      do j=1,N
        do i=1,N
          U(i,j)=VT(i,j)
        enddo
        D(j)=DIA(j)
      enddo
      end SUBROUTINE

C--------------------------------------------------------------------
      SUBROUTINE STAT_COVAR(ID,IT,DIM,COV,MEAN)
C return covariance matrix or its inverse (IT<0) and mean
C--------------------------------------------------------------------
      integer, intent(in) :: ID,IT,DIM
      REAL(KIND(1.D0)), intent(out) :: COV(DIM,DIM),MEAN(DIM)
      REAL(KIND(1.D0)) :: AUX(MD,MD),VT(MD,MD),DIA(MD)
      INTEGER :: I,J,N

      if (CV(ID).NC.LE.10) then
        call MSG_ERROR('STAT_COVAR','not enough events for covariance estimation',0,1)
        return
      endif

      N=CV(ID).ND
      call STAT_EVAL(ID)
      if (IT.lt.0) then
        call INVERT(N,CV(ID).C,MD,COV,DIM)
      else
        do j=1,N
          do i=1,N
            COV(i,j)=CV(ID).C(i,j)
          enddo
        enddo
      endif
      do i=1,N
        MEAN(i)=CV(ID).M(i)
      enddo
      end SUBROUTINE

C--------------------------------------------------------------------
      SUBROUTINE STAT_CONV(C,DIM)
C Calculate conversion matrix, C(DIM,DIM) for the two covariance matrices CV
C Task:
C ============
C Define: A(i,j)=Cov(X), B(i,j)=Cov(Y), X,Y being vectors (i=1..N)
C Let Y=C*X, then B=C*A*C^T
C We need to find C for given A and B
C Solution:
C =========
C Let's find U and W such that U^T*B*U = W^T*A*W = I (unit matrix)
C Then B=U-1^T*U-1 and hence C = W*U^-1
C--------------------------------------------------------------------
      include 'inout.inc'
      integer, intent(in) :: DIM
      REAL(KIND(1.D0)), intent(out) :: C(DIM,DIM)
      REAL(KIND(1.D0)) :: A(DIM,DIM),W(DIM,DIM),U(DIM,DIM),D(DIM),UINV(DIM,DIM)
      REAL(KIND(1.D0)) :: DET1,DET2,DETNAP,AUX(DIM,DIM),AUX1(DIM,DIM)
      INTEGER :: i,j,k,N



1     format(a20,6(2x,G10.4))
      call XML_RSXDUMP(sout,' ',1)

      N=CV(1).ND
      if (N.NE.CV(2).ND) then
        call MSG_ERROR('STAT_CONV','unequal matrix dimensions',0,1)
        return
      endif
      if ((CV(1).NC.LE.10).or.(CV(1).NC.LE.10)) then
        call MSG_ERROR('STAT_CONV','not enough events for covariance estimation',0,1)
        return
      endif

! decompose CV(1) to D=W^T*CV*W and normalize W
      call STAT_DIAG(1,DIM,W,D)
      do j=1,N
        do i=1,N
          W(i,j)=W(i,j)/SQRT(D(j))
        enddo
      enddo

      CALL XML_ARRAY(sout,N,'D(i) ','kx|ky|x|y',D(1))
  ! test W^T*COV1*W = 1
      do j=1,N
        do i=1,N
          U(i,j)=CV(1).C(i,j)
        enddo
      enddo
      call XML_MATRIX(sout,N,N,'COV1 ','kx|ky|x|y','kx|ky|x|y',U(1,1),DIM)
      DET1=DETNAP(U,DIM,N)
      write(*,*) 'det(Cov1)=',DET1
      call BTAB(U,W,AUX,DIM,N)
      write(*,*)
      do i=1,N
        write(*,1) 'W^T*COV1*W ',(AUX(i,j),j=1,N)
      enddo

! decompose CV(2) to D=A^T*CV*A, normalize A and get U=A^-1
! use A^-1=A^T
      call STAT_DIAG(2,DIM,A,D)
c      call INVERT(N,A,DIM,UINV,DIM)

c      write(*,*) 'STAT_DIAG OK ',N,DIM
      do j=1,N
        do i=1,N
          UINV(i,j)=A(j,i)*SQRT(D(i))
        enddo
      enddo
c      write(*,*) 'UINV OK '
      call XML_ARRAY(sout,N,'D(i) ','kx|ky|x|y',D(1))
c      call XML_MATRIX(sout,DIM,DIM,'U(i,j) ','kx|ky|x|y','kx|ky|x|y',U(1,1))


  ! test U^T*COV2*U = 1
    ! U has been inverted, therefore we have to invert it back
c      write(*,*) 'before INVERT '
      call INVERT(N,UINV,DIM,A,DIM)
      do j=1,N
        do i=1,N
          AUX1(i,j)=CV(2).C(i,j)
        enddo
      enddo
      call XML_MATRIX(sout,N,N,'COV2 ','kx|ky|x|y','kx|ky|x|y',AUX1(1,1),DIM)
c      write(*,*) 'before DETNAP '
       DET2=DETNAP(AUX1,DIM,N)

      write(*,*) 'det(Cov2)=',DET2
      write(*,*) 'det(Cov2)/det(Cov1)=',DET2/DET1

      call BTAB(AUX1,A,AUX,DIM,N)
      write(*,*)
      do i=1,N
        write(*,1) 'U^T*COV2*U ',(AUX(i,j),j=1,N)
      enddo


! construct C=W*UINV
      call MXM(1,N,DIM,W,UINV,AUX)
      do j=1,N
        do i=1,N
          C(i,j)=AUX(j,i)
        enddo
      enddo

c     DO j=1,N
c       do i=1,N
c         C(i,j)=0.D0
c         DO k=1,N
c            C(i,j)=C(i,j)+W(i,k)*U(k,j)
c          enddo
c        enddo
c      enddo
      write(*,*) 'det(C)=',DETNAP(C,DIM,N)


! check transport matrix, C*COV1*C^T = COV2

      do j=1,N
        do i=1,N
          U(i,j)=CV(1).C(i,j)
        enddo
      enddo
      call BABT(U,C,A,DIM,N)
      do j=1,N
        do i=1,N
          A(i,j)=A(i,j)-CV(2).C(i,j)
        enddo
      enddo

      write(*,*)
      do i=1,N
        write(*,1) 'C*COV1*C^T - COV2 ',(A(i,j),j=1,N)
      enddo


c      write(*,*)
c      do i=1,N
c        write(*,1) 'C1 ',(CV(1).C(i,j),j=1,N)
c      enddo
c      write(*,*)
c      do i=1,N
c        write(*,1) 'C2 ',(CV(2).C(i,j),j=1,N)
c      enddo

c      call XML_MATRIX(sout,DIM,DIM,'check transformation: ','kx|ky|x|y','kx|ky|x|y',A(1,1))

      call XML_RSXDUMP(sout,' ',0)
      END SUBROUTINE


      END MODULE COVARIANCES
