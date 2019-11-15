!//////////////////////////////////////////////////////////////////////
!////  $Id: optimization.f,v 1.16 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.16 $
!////     $Date: 2019/08/16 17:16:25 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Numerical optimization, basic subroutines and interfaces
!////
!//////////////////////////////////////////////////////////////////////
      MODULE OPTIMIZATION
      use XMLINFO
      implicit none
      SAVE
      PRIVATE

      public LMOPT,LMFIT,MGAUSS,GETCHI2

      contains

!------------------------------------------------------------
      SUBROUTINE GETDERIV1(A,DELTA,FUNCTN,DER1,DER2)
!    1st and 2nd derivative for 1 parameter
!------------------------------------------------------------
      REAL,intent(inout) :: A,DELTA
      REAL(kind(1.D0)),intent(out) :: DER1,DER2
      INTERFACE
        REAL FUNCTION FUNCTN(PAR)
          REAL,intent(inout) :: PAR(:)
        END FUNCTION FUNCTN
      END INTERFACE
      REAL(kind(1.D0)) :: Y1,Y2,Y0
      REAL :: B(1)
      DER1=0
      DER2=0
      IF (DELTA.NE.0) THEN
10       B(1) = A+DELTA
         Y1 = functn(B)
         B(1) = A-DELTA
         Y2 = functn(B)
         IF(ABS((Y1-Y2)/(Y1+Y2)).GT.0.4) THEN
            DELTA=DELTA/2.
c            write(*,*) DELTA,Y1,Y2,ABS((Y1-Y2)/(Y1+Y2))
            GOTO 10
         ENDIF
         B(1) = A
         Y0 = functn(B)
         DER1=(Y1-Y2)/2/DELTA
         DER2=(Y1+Y2-2*Y0)/DELTA/DELTA
      ENDIF
      END SUBROUTINE GETDERIV1

!-------------------------------------------------------------------
      REAL(kind(1.D0)) FUNCTION GETDERIV(IDER,J,K,A,DELTAA,N,FUNCTN)
!-------------------------------------------------------------------
      INTEGER,intent(in) :: IDER,J,K,N
      REAL,intent(inout) :: A(:),DELTAA(:)
      INTERFACE           ! for EXTERNAL procedures
        REAL FUNCTION FUNCTN(PAR)
          REAL,intent(inout) :: PAR(:)
        END FUNCTION FUNCTN
      END INTERFACE
      REAL :: B(N)
      REAL(kind(1.D0)) :: Y1,Y2,Y0,W1,W2,DELTA,DELTB,RES
      RES=0.D0
      IF (IDER.EQ.0) THEN
        RES=functn(A)
      ELSE IF (IDER.EQ.1) THEN
        B(1:N)=A(1:N)
        IF(DELTAA(J).NE.0) THEN
10        DELTA = DELTAA(J)
          B(J) = A(J)+DELTA
          Y1 = functn(B)
          B(J) = A(J)-DELTA
          Y2 = functn(B)
!1         format(a,I3,2x,5(G12.5,1x))
!          write(*,1) 'GETDERIV ',J,Y1,Y2,A(J)+DELTA,A(J)-DELTA,(Y1-Y2)/2/DELTA
          IF (ABS((Y1-Y2)/(Y1+Y2)).GT.0.4) THEN
            DELTAA(J)=DELTAA(J)/2.
            GOTO 10
          ENDIF
          RES=(Y1-Y2)/2/DELTA
        ENDIF
      ELSE IF (IDER.EQ.2) THEN
        IF ((DELTAA(J).NE.0).AND.(DELTAA(K).NE.0)) THEN
          B(1:N)=A(1:N)
          DELTA = DELTAA(J)
          DELTB = DELTAA(K)
          IF (J.EQ.K) THEN
            B(J) = A(J)+DELTA
            Y1 = functn(B)
            B(J) = A(J)-DELTA
            Y2 = functn(B)
            B(J) = A(J)
            Y0 = functn(B)
           RES=(Y1+Y2-2*Y0)/DELTA/DELTA
         ELSE
           B(K) = A(K)+DELTB
           B(J) = A(J)+DELTA
           Y1 = functn(B)
           B(J) = A(J)-DELTA
           Y2 = functn(B)
           W1=(Y1-Y2)/2/DELTA
           B(K) = A(K)-DELTB
           B(J) = A(J)+DELTA
           Y1 = functn(B)
           B(J) = A(J)-DELTA
           Y2 = functn(B)
           W2=(Y1-Y2)/2/DELTA
           RES=(W1-W2)/2/DELTB
         ENDIF
        ENDIF
      ENDIF
      GETDERIV=RES
      END FUNCTION GETDERIV

!---------------------------------------------------------------------
      SUBROUTINE LMOPT(FUNCTN,A,AMI,AMA,INCA,NP,TOL,ISIL)
!  find minimum of FUNCTN function by Levenberg-Marquardt algorithm
!  A(i)        ...  parameters
!  AMI,AMA(i)  ...  min/max parameter values
!  INCA(i)     ...  minimum increments for numerical derivatives of A
!  TOL         ...  tolerance indicator
!  ISIL        ...  silence level (no output for ISIL>1)
!  result is indicated through the  COMMON /ERROR/IERR
!  IERR=-2 ... max. iteration reached
!  IERR=-1 ... lambda diverges
!  IERR=1  ... O.K. exact match (linear fit)
!  IERR=2  ... O.K. TOL reached
!----------------------------------------------------------------------
      use SIMATH
      use IO
      INTEGER,intent(in) :: NP,ISIL
      REAL,intent (inout) :: A(:)
      REAL,intent (in) :: AMI(:),AMA(:),INCA(:),TOL
      integer,PARAMETER :: M=192

      INTEGER :: MAXIT,J,K,INS(M),IT,IERR
      REAL :: FLAMDA
      REAL(kind(1.D0)) :: ARRAY(M,M),ALPHA(M,M),BETA(M),AN(M,M)
      REAL :: B(M),DELTAA(M)
      REAL :: CHISQR,CHISQ1,DB,KSI
    !  REAL*8 GETDERIV
    !  REAL*4 FUNCTN
    !  external functn
      INTERFACE           ! for EXTERNAL procedures
        REAL FUNCTION FUNCTN(PAR)
          REAL,intent(inout) :: PAR(:)
        END FUNCTION FUNCTN
      END INTERFACE
1     format(/,'CHISQR1: ',G12.5,2x,'CHISQR: ',G12.5,2x,'PAR: ',6(G12.5,2x))
2     format('ITER: ',I5,2x,'FLAMDA: ',G12.5)
       if (isil<1) call XML_PROGRESS(SMES,0,10,1,'Fitting in progress')

C/// initializations
      IERR=0
      IT=0      ! iter. counter
      MAXIT=100 ! max. 100 iterations
      DO J=1,NP
        DELTAA(J)=ABS(A(J)*TOL) ! derive increments from TOL
        IF(DELTAA(J).LT.INCA(J)) DELTAA(J)=INCA(J)
        IF (A(J).GT.AMA(J).OR.A(J).LT.AMI(J)) THEN
          RETURN ! initial values outside limits !!
        ENDIF
      ENDDO
      FLAMDA=0.1
C/// Get initial Chi^2
      CHISQR = FUNCTN(A)
      !write(SOUT,1) 0.0D0,CHISQR,(A(J),J=1,MIN(NP,6))
C/// Start main cycle
41    if (isil.lt.1)  write(SOUT,1) CHISQ1,CHISQR,(A(J),J=1,MIN(NP,6))
      CHISQ1 = CHISQR
      !write(*,*) 'TOL=',TOL

C/// Calculate derivatives:
      IF(NP.EQ.1) THEN
        CALL GETDERIV1(A(1),DELTAA(1),FUNCTN,BETA(1),ALPHA(1,1))
      ELSE
        DO J=1,NP
            BETA(J)=GETDERIV(1,J,J,A,DELTAA,NP,FUNCTN)
            DO K=1,J
              ALPHA(J,K)=GETDERIV(2,J,K,A,DELTAA,NP,FUNCTN)
            END DO
        END DO
      ENDIF

C/// symetrize ALPHA
      DO J=1,NP
        DO K=1,J
          ALPHA(K,J) = ALPHA(J,K)
        END DO
      END DO

C/// check zeros on digonal
      DO  J=1,NP
          INS(J) = 0
          IF(ABS(ALPHA(J,J)).LT.1.E-19) INS(J)=1
      END DO

C/// create Hessian matrix
71    DO J=1,NP
        DO K=1,NP
          IF((INS(J)+INS(K)).EQ.0) THEN
            AN(J,K)=SQRT(ABS(ALPHA(J,J))*ABS(ALPHA(K,K)))
            ARRAY(J,K) = ALPHA(J,K)/AN(J,K)
          ELSE
            AN(J,K)=1.D0
            ARRAY(J,K) = 0.D0
          ENDIF
        END DO
        ARRAY(J,J) = 1.D0+FLAMDA
      END DO

C/// invert Hessian matrix
      J=INVERT2(NP,ARRAY,M)
      !write(*,*) 'LMOPT invert=',J

C/// increment with limits checking
      KSI=1.0
      DO J=1,NP
10      DB = 0.0
        DO K=1,NP
          IF((INS(J)+INS(K)).EQ.0) THEN
            DB=DB-KSI*BETA(K)*ARRAY(J,K)/AN(J,K)
          ENDIF
        END DO
        IF (A(J)+DB.GT.AMA(J).OR.A(J)+DB.LT.AMI(J)) THEN
          KSI=KSI/2.0
          GOTO 10
        ENDIF
!3       format('J,KSI,A,DA: ',5(G13.5,2x))
!4       format(a,5(G13.5,2x))
!5       format(a,5(I5,2x))
!        if (J==1) write(*,3) J,KSI,A(J),DB,FLAMDA
!        if (J==1) write(*,4) 'BETA: ',BETA(1:NP)
!        if (J==1) write(*,4) 'DELTAA: ',DELTAA(1:NP)
        B(J)=A(J)+DB
      END DO

C// get new CHISQR
      CHISQR = FUNCTN(B)
      !write(SOUT,1) CHISQ1,CHISQR,(B(J),J=1,MIN(NP,6))
      !write(SOUT,2) IT,FLAMDA


C// check stop conditions
      IF(CHISQR.LT.1E-20) GOTO 110 ! exact match ... speed-up linear fit
      IF ((CHISQ1-CHISQR).LE.0.D0) THEN
        FLAMDA = 5.*FLAMDA
        IF(FLAMDA.LE.1000.) THEN ! try again with larger FLAMBDA
           GOTO 71
        ELSE
           IERR=-1
           GOTO 111  ! exit if FLAMBDA is too large
        ENDIF
      ENDIF

C//  accept new values and continue
110   DO J=1,NP
         A(J) = B(J)
      END DO
      FLAMDA = FLAMDA/5.
      IT=IT+1

111   IF (CHISQR.LT.1E-20) IERR=1 ! exit on exact match
      IF (ABS(ABS(CHISQ1/CHISQR)-1.).LT.0.1*TOL) IERR=2 ! exit on < 0.1*TOL change
      IF (IT.GE.MAXIT) IERR=-2 ! exit on maximum iteration number

      if (isil<1) call XML_PROGRESS(SMES,1,10,mod(IT,10),'Fitting in progress')

      IF (IERR.EQ.0) GOTO 41

      SELECT CASE (IERR)
      CASE (1)
        if (isil.lt.1)  WRITE(SOUT,*) 'OK, exact match'
      CASE (2)
        !write(*,*) 'RESULT: ',CHISQ1,CHISQR,ABS(ABS(CHISQ1/CHISQR)-1.),FLAMDA
        if (isil.lt.1)  WRITE(SOUT,*) 'OK, exit on TOL limit'
      CASE (-1)
        if (isil.lt.2)  WRITE(SOUT,*) 'lambda>1000'
      CASE (-2)
        if (isil.lt.2)  WRITE(SOUT,*) 'iteration limit exceeded'
      END SELECT
      if (isil<1) call XML_PROGRESS(SMES,2,10,10,'Fitting in progress')

      END SUBROUTINE LMOPT


!---------------------------------------------------------------------
      integer function LMFIT(FUNCTN,X,Y,DY,N,OPT,A,DA,AMI,AMA,INCA,NP,TOL,ISIL)
!  fit FUNCTN function by Levenberg-Marquardt algorithm to X,Y data
!  A(i)        ...  parameters
!  AMI,AMA(i)  ...  min/max parameter values
!  INCA(i)     ...  minimum increments for numerical derivatives of A
!  TOL         ...  tolerance indicator
!  ISIL        ...  silence level (no output for ISIL>1)
!  result is indicated through the  COMMON /ERROR/IERR
!  IERR=-2 ... max. iteration reached
!  IERR=-1 ... lambda diverges
!  IERR=1  ... O.K. exact match (linear fit)
!  IERR=2  ... O.K. TOL reached
!----------------------------------------------------------------------
      use SIMATH
      use IO
      INTEGER,intent(in) :: N,NP,ISIL
      REAL,intent (inout) :: A(:),DA(:)
      REAL,intent (in) :: X(:),Y(:),DY(:),OPT(:)
      REAL,intent (in) :: AMI(:),AMA(:),INCA(:),TOL
      integer,PARAMETER :: M=192

      INTEGER :: MAXIT,I,J,K,INS(M),IT,IERR
      REAL :: FLAMDA
      REAL(kind(1.D0)) :: ARRAY(M,M),ALPHA(M,M),BETA(M)
      REAL :: B(M),C(M),DELTAA(M),DER(M),F1,F2
      REAL :: CHISQR,CHISQ1,DB,KSI

      INTERFACE           ! for EXTERNAL procedures
        REAL FUNCTION FUNCTN(PAR,OP,NP,V)
          REAL,intent(inout) :: PAR(:)
          integer,intent(in) :: NP
          REAL,intent(in) :: OP(:),V
        END FUNCTION FUNCTN
      END INTERFACE

1     format(/,'CHISQR1: ',G12.5,2x,'CHISQR: ',G12.5,2x,'PAR: ',6(G12.5,2x))
      if (isil<2) call XML_PROGRESS(SMES,0,10,1,'Fitting in progress')
      !write(*,*) 'LMFIT'

! initializations
      IERR=0
      IT=0      ! iter. counter
      MAXIT=50 ! max. 100 iterations
      DO J=1,NP
        DELTAA(J)=ABS(A(J)*TOL) ! derive increments from TOL
        IF(DELTAA(J).LT.INCA(J)) DELTAA(J)=INCA(J)
        IF (A(J).GT.AMA(J).OR.A(J).LT.AMI(J)) THEN
          write(*,*) 'initial values outside limits ! ',J,A(J),AMA(J),AMI(J)
          if (isil<2) call XML_PROGRESS(SMES,2,10,10,'Fitting in progress')
          goto 112 ! initial values outside limits !!
        ENDIF
      ENDDO
      FLAMDA=0.1
! Get initial Chi^2
      CHISQ1 = 0.0

    !  do I=1,NP
    !    if (A(I)>1.e6) write(*,1) 'LMFIT A: ',I,A(I)
    !  enddo

      CHISQR = GETCHI2(FUNCTN,X,Y,DY,N,OPT,A,NP)

! Start main cycle
!=================================================
41    if (isil.lt.1)  write(SOUT,1) CHISQ1,CHISQR,(A(J),J=1,MIN(NP,3))
      CHISQ1 = CHISQR

! calculate derivatives:
      BETA=0.0
      ALPHA=0.0
      C(1:NP)=A(1:NP)
      DO I=1,N
      if (DY(I)>0.0) then
        DO J=1,NP
          IF (DELTAA(J)>0.0) then
            C(J)=A(J)-DELTAA(J)/2.0
            F1=FUNCTN(C,OPT,NP,X(I))
            C(J)=A(J)+DELTAA(J)/2.0
            F2=FUNCTN(C,OPT,NP,X(I))
            DER(J)=(F2-F1)/DELTAA(J)
            C(J)=A(J)
          endif
        END DO
        do J=1,NP
          DO K=1,NP
            ALPHA(J,K)=ALPHA(J,K)+DER(J)*DER(K)/DY(I)**2
          enddo
          F1=FUNCTN(C,OPT,NP,X(I))
          BETA(J)=BETA(J)+DER(J)*(F1-Y(I))/DY(I)**2
          !if ((J==1).and.(Y(I)>1.0)) write(*,*) 'DER,BETA: ',DER(J),BETA(J),F1,Y(I)
        enddo
      endif
      ENDDO

! check zeros on digonal
      DO  J=1,NP
        INS(J) = 0
        IF(ABS(ALPHA(J,J)).LT.1.E-19) INS(J)=1
      END DO

! create Hessian matrix
71    DO J=1,NP
        DO K=1,NP
          IF((INS(J)+INS(K)).EQ.0) THEN
            ARRAY(J,K) = ALPHA(J,K)
          ELSE
            ARRAY(J,K) = 0.D0
          ENDIF
        END DO
        ARRAY(J,J) = ARRAY(J,J)*(1.D0+FLAMDA)
      END DO

! invert Hessian matrix
      J=INVERT2(NP,ARRAY,M)
      if (J.ne.0) then
        IERR=-3
        goto 112
      endif
      !write(*,*) 'LMOPT invert=',J

! increment with limits checking
      KSI=1.0
      DO J=1,NP
10      DB = 0.0
        DO K=1,NP
          IF((INS(J)+INS(K)).EQ.0) THEN
            DB=DB-KSI*BETA(K)*ARRAY(J,K)
          ENDIF
        END DO
        IF ((A(J)+DB > AMA(J)).OR.(A(J)+DB < AMI(J))) THEN
          KSI=KSI/2.0
          GOTO 10
        ENDIF
!2       format('STEP: ',6(G13.5,2x))
!        if (J==1) write(*,2) J,A(J),DB,DER(J),BETA(J),ARRAY(J,J)
        B(J)=A(J)+DB
      END DO

! get new CHISQR
      CHISQR = GETCHI2(FUNCTN,X,Y,DY,N,OPT,B,NP)
      !write(*,*) 'ITER: ',CHISQ1,CHISQR,FLAMDA

C// check stop conditions
      IF(CHISQR.LT.1E-20) GOTO 110 ! exact match ... speed-up linear fit
      IF ((CHISQ1-CHISQR).LE.0.D0) THEN
        FLAMDA = 5.*FLAMDA
        IF(FLAMDA.LE.1000.) THEN ! try again with larger FLAMBDA
           GOTO 71
        ELSE
           IERR=-1
           GOTO 111  ! exit if FLAMBDA is too large
        ENDIF
      ENDIF

C//  accept new values and continue
110   A(1:NP)=B(1:NP)
      FLAMDA = FLAMDA/5.
      IT=IT+1

111   IF (CHISQR.LT.1E-20) IERR=1 ! exit on exact match
      IF (ABS(ABS(CHISQ1/CHISQR)-1.).LT.0.01*TOL) IERR=2 ! exit on < 0.1*TOL change
      IF (IT.GE.MAXIT) IERR=-2 ! exit on maximum iteration number

      if (isil<2) call XML_PROGRESS(SMES,1,10,mod(IT,10),'Fitting in progress')

      IF (IERR.EQ.0) GOTO 41

      DA(1:NP)=0.0
      do J=1,NP
        if (abs(ALPHA(J,J))>0.D0) DA(J)=1.D0/sqrt(abs(ALPHA(J,J)))
      enddo

112   SELECT CASE (IERR)
      CASE (1)
        if (isil.lt.1)  WRITE(SOUT,*) 'OK, exact match'
      CASE (2)
        !write(*,*) 'RESULT: ',CHISQ1,CHISQR,ABS(ABS(CHISQ1/CHISQR)-1.),FLAMDA
        if (isil.lt.1)  WRITE(SOUT,*) 'OK, exit on TOL limit'
      CASE (-1)
        if (isil.lt.2)  WRITE(SOUT,*) 'END, lambda>1000'
      CASE (-2)
        if (isil.lt.2)  WRITE(SOUT,*) 'END, iteration limit exceeded'
      CASE (-3)
        if (isil.lt.2)  WRITE(SOUT,*) 'ERROR: can''t invert matrix'
      END SELECT
      if (isil<2) call XML_PROGRESS(SMES,2,10,10,'Fitting in progress')
      LMFIT=IERR
      END function  LMFIT

!-------------------------------------------------------------------
      REAL FUNCTION GETCHI2(FUNCTN,X,Y,DY,N,OPT,PAR,NP)
! CHI2 for FUNCTN and given data X,Y,DY
! NP .. number of parameters
! PAR .. array with dimension 3xN
! Initial values should already be set in PAR
!-------------------------------------------------------------------
      REAL,intent(in) :: X(:),Y(:),DY(:),OPT(:)
      REAL :: PAR(:)
      integer,intent(in) :: N,NP
      REAL(KIND(1.D0)) :: Z
      REAL :: CHI
      INTEGER :: I,nd
      INTERFACE           ! for EXTERNAL procedures
        REAL FUNCTION FUNCTN(PAR,OP,NP,V)
          REAL,intent(inout) :: PAR(:)
          integer,intent(in) :: NP
          REAL,intent(in) :: OP(:),V
        END FUNCTION FUNCTN
      END INTERFACE
      CHI=0.0
      nd=0
      if (np<=0) then
        GETCHI2=0.0
        return
      endif
      DO I=1,N
        IF(DY(I).gt.0) THEN
          nd=nd+1
          Z=FUNCTN(PAR,OPT,NP,X(I))
          CHI=CHI+((Y(I)-Z)/DY(I))**2
        ENDIF
      ENDDO
      if (nd>np) then
        CHI=CHI/(nd-np)
      else
        CHI=0.0
      endif

      GETCHI2=CHI
      END FUNCTION GETCHI2

!-------------------------------------------------------------------
      REAL FUNCTION MGAUSS(PAR,NP,X)
! Get value of multiple Gauss function at X
! NP .. number of parameters
! PAR .. array with dimension 3xNP
!-------------------------------------------------------------------
      use CONSTANTS
      REAL,intent(inout) :: PAR(:)
      REAL,intent(in) :: X
      integer,intent(in) :: NP
      REAL(KIND(1.D0)) :: Z
      INTEGER :: K
      Z=0.D0
      DO K=0,NP-3,3
        PAR(K+3)=ABS(PAR(K+3))
        IF (ABS(PAR(K+3)).LT.1E-10) PAR(K+3)=1E-10
        Z=Z+PAR(K+1)*EXP(-0.5*(X-PAR(K+2))**2/(PAR(K+3)/R8LN2)**2)
      enddo
      MGAUSS=Z
      END FUNCTION MGAUSS



      end module OPTIMIZATION
