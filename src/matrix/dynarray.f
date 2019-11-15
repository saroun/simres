!//////////////////////////////////////////////////////////////////////
!////  $Id: dynarray.f,v 1.3 2010/12/11 20:44:02 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2010, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.3 $
!////     $Date: 2010/12/11 20:44:02 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Dynamically allocated arrays and basic operations
!////  The arrays are defined as types incl. dimensions
!////  Not optimized for speed, but simplifies the syntax of matrix operations.
!////
!//////////////////////////////////////////////////////////////////////
      MODULE DYNARRAY
      implicit none
      PRIVATE

      type TDARRAYI
        integer :: DX  ! allocated dimension
        integer :: NX  ! used dimension
        integer, allocatable :: X(:)
      end type TDARRAYI

      type TDARRAY
        integer :: DX  ! allocated dimension
        integer :: NX  ! used dimension
        real(kind(1.D0)), allocatable :: X(:)
      end type TDARRAY

      type TDMATRIX
        integer :: DX,DY  ! allocated dimensions
        integer :: NX,NY  ! used dimensions
        real(kind(1.D0)), allocatable :: M(:,:)
      end type TDMATRIX

      interface DYN_NEWARRAY
        module procedure DYN_NEWARRAY_I
        module procedure DYN_NEWARRAY_D
        module procedure DYN_NEWARRAY_FROM_I
        module procedure DYN_NEWARRAY_FROM_D
      end interface  DYN_NEWARRAY

      interface DYN_NEWMATRIX
        module procedure DYN_NEWMATRIX_D
        module procedure DYN_NEWMATRIX_FROM_D
      end interface  DYN_NEWMATRIX

      interface DYN_DISPOSE
        module procedure DYN_DISPOSE_I
        module procedure DYN_DISPOSE_X
        module procedure DYN_DISPOSE_M
      end interface  DYN_DISPOSE

      public TDARRAY,TDMATRIX,TDARRAYI
      public DYN_NEWARRAY,DYN_NEWMATRIX,DYN_UNIMAT,DYN_DISPOSE,DYN_DETJAC
      public DYN_ADDBTAV,DYN_ADDBTAB,DYN_MXM,DYN_JACOBI,DYN_MltpBySubmatrix

      CONTAINS


!---------------------------------------------------------
      subroutine DYN_NEWARRAY_D(X,NX)
! allocate new array of dim=NX if necessary and set X%X=0
!----------------------------------------------------------
      integer,intent(in) :: NX
      type(TDARRAY) :: X
      integer :: ierr
      ierr=0
    ! allocate when necessary
      if (.not.allocated(X%X).or.(X%NX.ne.NX)) then
        X%NX=0
        X%DX=0
        if (allocated(X%X)) deallocate(X%X,STAT=ierr)
        allocate(X%X(1:NX),STAT=ierr)
        if (ierr.eq.0) X%DX=NX
      endif
    ! set X%X=0
      if (ierr.eq.0) then
        X%X(1:NX)=0.D0
        X%NX=NX
      endif
      end subroutine DYN_NEWARRAY_D

!---------------------------------------------------------
      subroutine DYN_NEWARRAY_FROM_D(X,A,NA)
! allocate new array of dim=NX if necessary and set X%X=A
!----------------------------------------------------------
      real(kind(1.D0)), intent(in) :: A(:)
      integer,intent(in) :: NA
      type(TDARRAY) :: X
      call DYN_NEWARRAY_D(X,NA)
      X%X(1:X%NX)=A(1:X%NX)
      end subroutine DYN_NEWARRAY_FROM_D

!---------------------------------------------------------
      subroutine DYN_NEWARRAY_I(X,NX)
! allocate new array of dim=NX if necessary and set X%X=0
!----------------------------------------------------------
      integer,intent(in) :: NX
      type(TDARRAYI) :: X
      integer :: ierr
      ierr=0
    ! allocate when necessary
      if (.not.allocated(X%X).or.(X%NX.ne.NX)) then
        X%NX=0
        X%DX=0
        if (allocated(X%X)) deallocate(X%X,STAT=ierr)
        allocate(X%X(1:NX),STAT=ierr)
        if (ierr.eq.0) X%DX=NX
      endif
    ! set X%X=0
      if (ierr.eq.0) then
        X%X(1:NX)=0.D0
        X%NX=NX
      endif
      end subroutine DYN_NEWARRAY_I

!---------------------------------------------------------
      subroutine DYN_NEWARRAY_FROM_I(X,A,NA)
! allocate new array of dim=NX if necessary and set X%X=A
!----------------------------------------------------------
      integer,intent(in) :: A(:)
      integer,intent(in) :: NA
      type(TDARRAYI) :: X
      call DYN_NEWARRAY_I(X,NA)
      X%X(1:X%NX)=A(1:X%NX)
      end subroutine DYN_NEWARRAY_FROM_I


!---------------------------------------------------------
      subroutine DYN_NEWMATRIX_D(M,NX,NY)
! allocate new 2D matrix of dim=NXxNY if necessary and set M%M=0
!----------------------------------------------------------
      integer,intent(in) :: NX,NY
      type(TDMATRIX) :: M
      integer :: ierr
      ierr=0
    ! allocate when necessary
      if (.not.allocated(M%M).or.(M%NX.ne.NX).or.(M%NY.ne.NY)) then
        M%NX=0
        M%NY=0
        M%DX=0
        M%DY=0
        if (allocated(M%M)) deallocate(M%M,STAT=ierr)
        allocate(M%M(1:NX,1:NY),STAT=ierr)
        if (ierr.eq.0) then
          M%DX=NX
          M%DY=NY
        endif
      endif
    ! set A%X=0
      if (ierr.eq.0) then
        M%M(1:NX,1:NY)=0.D0
        M%NX=NX
        M%NY=NY
      endif
      end subroutine DYN_NEWMATRIX_D

!---------------------------------------------------------
      subroutine DYN_NEWMATRIX_FROM_D(M,A,DX,NX,NY)
! allocate new array of dim=NX if necessary and set X%X=A
!----------------------------------------------------------
      integer,intent(in) :: DX,NX,NY
      integer,intent(in) :: A(DX,NY)
      type(TDMATRIX) :: M
      integer :: j
      call DYN_NEWMATRIX_D(M,NX,NY)
      do j=1,NY
        M%M(1:M%NX,j)=A(1:M%NX,j)
      enddo
      end subroutine DYN_NEWMATRIX_FROM_D


!---------------------------------------------------------
      subroutine DYN_DISPOSE_M(M)
! dispose 2D matrix
!----------------------------------------------------------
      type(TDMATRIX) :: M
      integer :: ierr
      if (allocated(M%M)) deallocate(M%M,STAT=ierr)
      M%NX=0
      M%NY=0
      M%DX=0
      M%DY=0
      end subroutine DYN_DISPOSE_M

!---------------------------------------------------------
      subroutine DYN_DISPOSE_X(X)
! dispose 2D matrix
!----------------------------------------------------------
      type(TDARRAY) :: X
      integer :: ierr
      if (allocated(X%X)) deallocate(X%X,STAT=ierr)
      X%NX=0
      X%DX=0
      end subroutine DYN_DISPOSE_X

!---------------------------------------------------------
      subroutine DYN_DISPOSE_I(X)
! dispose 2D matrix
!----------------------------------------------------------
      type(TDARRAYI) :: X
      integer :: ierr
      if (allocated(X%X)) deallocate(X%X,STAT=ierr)
      X%NX=0
      X%DX=0
      end subroutine DYN_DISPOSE_I

!---------------------------------------
      subroutine DYN_MXM(IT,A,B,C)
! C = A.B
! requires A%NY=B%NX, not checked
!---------------------------------------
      integer,intent(in) :: IT
      type(TDMATRIX) :: A,B,C
      integer :: i,j,k
      call DYN_NEWMATRIX(C,A%NX,B%NY)
      do i=1,A%NX
        do j=1,A%NX
          if (IT.gt.0) then
            do k=1,A%NY
              C%M(i,j)=C%M(i,j)+A%M(i,k)*B%M(k,j)
            enddo
          else
            do k=1,A%NY
              C%M(i,j)=C%M(i,j)+A%M(k,i)*B%M(k,j)
            enddo
          endif
        enddo
      enddo
      end subroutine DYN_MXM


!---------------------------------------
      subroutine DYN_MXV(M,V,C)
! C = M.V
! requires A%NY=B%NX, not checked
!---------------------------------------
      type(TDMATRIX) :: M
      type(TDARRAY) :: V,C
      integer :: i,j,k
      call DYN_NEWARRAY(C,M%NX)
      do i=1,M%NX
        do k=1,M%NY
          C%X(i)=C%X(i)+M%M(i,k)*V%X(k)
        enddo
      enddo
      end subroutine DYN_MXV

!------------------------------------------------------------
      subroutine DYN_AddBTAB(C,B,A,INDX,NA)
! Add a submatrix A to the matrix C, using conversion matrix B
! C = C + B^T.A.B
! NA ... dimension of the variables subspace for A,V
! A is defined on a subspace of B,C space using INDX = array of indices
! C,B must have appropriate dimensions already defined - not checked here
!------------------------------------------------------------
      type(TDMATRIX) :: B,C
      real(kind(1.D0)),intent(in) :: A(:,:)
      integer,intent(in) :: INDX(:),NA
      integer :: i,j,k,m,xk,xm
      real(kind(1.D0)) :: suma
      do i=1,C%NX
        do j=1,C%NY
          do k=1,NA
            xk=INDX(k)
            suma=0.D0
            do m=1,NA
              xm=INDX(m)
              suma=suma+A(k,m)*B%M(xm,j)
            enddo
            C%M(i,j)=C%M(i,j)+B%M(xk,i)*suma
          enddo
        enddo
      enddo
      end subroutine DYN_AddBTAB

!------------------------------------------------------------
      subroutine DYN_AddBTAV(C,B,A,V,INDX,NA)
! Add a submatrix A to the matrix C, using conversion matrix B
! C = C + B^T.A.V
! NA ... dimension of the variables subspace for A,V
! A,V are defined on a subspace of B,C space using INDX = array of indices
! C,B must have appropriate dimensions already defined - not checked here
!------------------------------------------------------------
      type(TDMATRIX) :: B
      type(TDARRAY) :: C
      real(kind(1.D0)),intent(in) :: A(:,:),V(:)
      integer,intent(in) :: INDX(:),NA
      integer :: i,j,k,m,xk,xm
      real(kind(1.D0)) :: suma
      do i=1,C%NX
          do k=1,NA
            xk=INDX(k)
            suma=0.D0
            do m=1,NA
              xm=INDX(m)
              suma=suma+A(k,m)*V(m)
            enddo
            C%X(i)=C%X(i)+B%M(xk,i)*suma
          enddo
      enddo
      end subroutine DYN_AddBTAV

!------------------------------------------------------------
      subroutine DYN_MltpBySubmatrix(ISIDE,A,B,INDX,NX)
! Multiply matrix A by a submatrix B
! from left (A=B.A, ISIDE=0) or from right (A=A.B, ISIDE=1)
! B is defined on a subspace of A-space using INDX = array of indices
!------------------------------------------------------------
      integer, intent(in) :: ISIDE,NX
      type(TDMATRIX) :: A
      real(kind(1.D0)),intent(in) :: B(:,:)
      integer,intent(in) :: INDX(:)
      real(kind(1.D0)) :: C(NX)
      integer :: i,j,k,xk,xj
      select case(ISIDE)
      case(0)
  ! from left: A = B.A
        do i=1,A%NY ! scan columns
          do j=1,NX
            C(j)=0.D0
            do k=1,NX
              xk=INDX(k)
              C(j)=C(j)+B(j,k)*A%M(xk,i)
            enddo
          enddo
          do j=1,NX
            xj=INDX(j)
            A%M(xj,i)=C(j)
          enddo
        enddo
      case(1)
  ! from right: A = A.B
        do i=1,A%NX ! scan rows
          do j=1,NX
            C(j)=0.D0
            do k=1,NX
              xk=INDX(k)
              C(j)=C(j)+A%M(i,xk)*B(k,j)
            enddo
          enddo
          do j=1,NX
            xj=INDX(j)
            A%M(i,xj)=C(j)
          enddo
        enddo
      end select
      end subroutine DYN_MltpBySubmatrix


!--------------------------------------------------------
      SUBROUTINE DYN_UNIMAT(U,N)
! create NxN identity matrix
!--------------------------------------------------------
      type(TDMATRIX) :: U
      integer,intent(in) :: N
      INTEGER :: I
      call DYN_NEWMATRIX(U,N,N)
      DO I=1,N
        U%M(I,I)=1.D0
      ENDDO
      END SUBROUTINE DYN_UNIMAT


!------------------------------------------
      SUBROUTINE DYN_JACOBI(AA,D,V,NROT)
! Diagonalize matrix AA,
! returns the diagonal, D
! and the transformation matrix, V
! D=V^T*AA*V
! (Numerical Recipes)
!------------------------------------------
      integer,intent(out) :: NROT
      type(TDMATRIX) :: AA,V
      type(TDARRAY) :: D
      type(TDMATRIX) :: A
      type(TDARRAY) :: B,Z
      integer :: IQ,IP,I
      real(KIND(1.D0)) :: G,SM,TRESH

      call DYN_NEWMATRIX(A,AA%NX,AA%NY)
      call DYN_NEWARRAY(B,AA%NX)
      call DYN_NEWARRAY(Z,AA%NX)
      Z%X=0.D0
      call DYN_UNIMAT(V,A%NX)
      DO IP=1,A%NX
        DO IQ=1,A%NX
          A%M(IP,IQ)=AA%M(IP,IQ)
        enddo
      enddo
      DO IP=1,A%NX
        B%X(IP)=A%M(IP,IP)
        D%X(IP)=B%X(IP)
      enddo
      NROT=0
      DO I=1,50
        SM=0.D0
        DO IP=1,A%NX-1
          DO IQ=IP+1,A%NX
            SM=SM+ABS(A%M(IP,IQ))
          enddo
        enddo
      !  write(*,*) 'JACOBI iter=',I,SM,A%NX,A%DX
        IF(SM.LE.1.D-30) GOTO 99
        IF(I.LT.4) THEN
          TRESH=0.2D0*SM/A%NX**2
        ELSE
          TRESH=0.D0
        ENDIF
        DO IP=1,A%NX-1
          DO IQ=IP+1,A%NX
            G=100.D0*ABS(A%M(IP,IQ))
            IF((I.GT.4).AND.(ABS(D%X(IP))+G.EQ.ABS(D%X(IP))).AND.(ABS(D%X(IQ))+G.EQ.ABS(D%X(IQ)))) THEN
              A%M(IP,IQ)=0.D0
            ELSE IF (ABS(A%M(IP,IQ)).GT.TRESH) THEN
              call DYN_JACOBI_ROT(IP,IQ,A,V,D,Z)
              NROT=NROT+1
            ENDIF
          enddo
        enddo
        DO IP=1,A%NX
          B%X(IP)=B%X(IP)+Z%X(IP)
          D%X(IP)=B%X(IP)
          Z%X(IP)=0.D0
        enddo
      enddo
      WRITE(*,*) 'JACOBI : 50 ITERATION SHOULD NEVER HAPPEN'
99    CONTINUE
      call DYN_DISPOSE(A)
      call DYN_DISPOSE(B)
      call DYN_DISPOSE(Z)
      END SUBROUTINE DYN_JACOBI


!------------------------------------------------------------
      SUBROUTINE DYN_JACOBI_ROT(IP,IQ,A,V,D,Z)
!  Make one Jacobi rotation for IP,IQ element
!------------------------------------------------------------
      integer,intent(in) :: IP,IQ
      type(TDMATRIX) :: A,V
      type(TDARRAY) :: D,Z
      real(KIND(1.D0)) :: G,H,S,C,T,THETA,TAU
      integer :: J
        G=100.D0*ABS(A%M(IP,IQ))
        H=D%X(IQ)-D%X(IP)
        IF (ABS(H)+G.EQ.ABS(H)) THEN
          T=A%M(IP,IQ)/H
        ELSE
          THETA=0.5D0*H/A%M(IP,IQ)
          T=1.D0/(ABS(THETA)+SQRT(1.+THETA**2))
          IF(THETA.LT.0.D0) T=-T
        ENDIF
        S=SQRT(1.D0+T**2)
        C=1.D0/S
        S=C*T
        TAU=S/(1.D0+C)
        H=T*A%M(IP,IQ)
        Z%X(IP)=Z%X(IP)-H
        Z%X(IQ)=Z%X(IQ)+H
        D%X(IP)=D%X(IP)-H
        D%X(IQ)=D%X(IQ)+H
        A%M(IP,IQ)=0.D0
        DO J=1,IP-1
          G=A%M(J,IP)
          H=A%M(J,IQ)
          A%M(J,IP)=G-S*(H+G*TAU)
          A%M(J,IQ)=H+S*(G-H*TAU)
        enddo
        DO J=IP+1,IQ-1
          G=A%M(IP,J)
          H=A%M(J,IQ)
          A%M(IP,J)=G-S*(H+G*TAU)
          A%M(J,IQ)=H+S*(G-H*TAU)
        enddo
        DO J=IQ+1,A%NY
          G=A%M(IP,J)
          H=A%M(IQ,J)
          A%M(IP,J)=G-S*(H+G*TAU)
          A%M(IQ,J)=H+S*(G-H*TAU)
        enddo
        DO J=1,V%NX
          G=V%M(J,IP)
          H=V%M(J,IQ)
          V%M(J,IP)=G-S*(H+G*TAU)
          V%M(J,IQ)=H+S*(G-H*TAU)
        enddo
      end SUBROUTINE DYN_JACOBI_ROT

!------------------------------------------
      REAL(KIND(1.D0)) FUNCTION DYN_DETJAC(A)
! Calculate matrix determinant using Jacobi diagonalization
!------------------------------------------
      type(TDMATRIX) :: A
      REAL(KIND(1.D0)) :: D
      type(TDMATRIX) :: VT
      type(TDARRAY) :: DIA
      integer :: NROT,I
      call DYN_NEWMATRIX_D(VT,A%NX,A%NY)
      call DYN_NEWARRAY(DIA,A%NX)
      call DYN_JACOBI(A,DIA,VT,NROT)
      D=1.D0
      do i=1,DIA%NX
        D=D*DIA%X(i)
      enddo
      DYN_DETJAC=D
      call DYN_DISPOSE(VT)
      call DYN_DISPOSE(DIA)
      end FUNCTION DYN_DETJAC


      end module DYNARRAY