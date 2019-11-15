!//////////////////////////////////////////////////////////////////////
!////  $Id: matrix.f,v 1.3 2010/12/11 20:44:02 saroun Exp $
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
!////  Module for matrix description of neutron transport
!////  Theoretical background:
!////  see I. Ionita, A.D. Stoica, J. Appl. Cryst. (2000). 33, 1067-1074
!////
!//////////////////////////////////////////////////////////////////////
!
!  Usage guide:
!  1) get actual number of random variables, N
!  2) initialize before parsing the instrument: call MX_INIT(N)
!  3) add individual components in sequence as they appear on the beamline:
!     in this order:
!     CALL MX_ADD_VARIABLES(CAP,IX,NV)
!     CALL MX_ADD_TRMAT(C,T,D,INDX,NT)
!     CALL MX_ADD_CONSTRAINT(C,B,D,DIMB,NB,INDX)
!     CALL MX_ADD_TRANSPORT(C,INDX,NC)
!  4) evaluate arrays to obtain final transmision and covariance matrices
!     as well as misfit array
!  5) dispose dynamic arrays at the end
!
      MODULE MATRIX
      use DYNARRAY
      use SIMATH
      use XMLINFO
      implicit none
      PRIVATE
      save

! shape constants = fwhm/sigma for various distribution profiles
      real(KIND(1.D0)), parameter :: s_step=0.288675134595
      real(KIND(1.D0)), parameter :: s_tria=0.408248290464
      real(KIND(1.D0)), parameter :: s_circ=0.25D0
      real(KIND(1.D0)), parameter :: s_gaus=0.424660900144

      integer :: NTRM ! actual number of random variables
      integer :: NTRM_PRIM ! dtto for primary section only
! arrays describing transport + recurrent formulae:
! C,T,D are arrays corresponding to n-th component in its LOCAL coordinates
! conversion of X from (i-1)-th exit coordinates to i-th entry coordinates is
!  X_i = C_i.X_(i-1) + D_i
! Arrays describing the whole instrument:
!  CMAT  ... coordinate conversion matrix: CMAT_i = C_i.CMAT_i-1
!  TMAT  ... transmission matrix, TMAT = TMAT + CMAT^T.Tn.CMAT
!            where Tn = C^T.T.C for n-th component
!  TVEC  ... misfit vector: TVEC = TVEC + CMAT_i^T.Tn.Dn
!            where Dn=C^T.T.D for n-th component
!  TCNT  ... constant term: TCNT = TCNT + D^T.T.D
! transmission function exponent in source coordinates, X
! 1/2.(X - TMAT^-1.TVEC)^T.TMAT.(X - TMAT^-1.TVEC) - 1/2.TCNT
      type(TDMATRIX) :: CMAT,TMAT
      type(TDARRAY) :: TVEC
      real(KIND(1.D0)) :: TCNT

! arrays describing constraints:
! constraint for a component in its entry coordinates is expressed as
! SUM_i {b(i).(X(i) - E(i))} = 0
! in Gaussian approximation, the exponent of partial constraint function of a component is
! 1/(2 lambda).(X - E_k)^T.B_k.(X - E_k)
! where B(i,j) = b(i)b(j)
! Arrays describing the whole instrument:
!  BMAT  ... constraint matrix, BMAT = BMAT + CMAT^T.Bn.CMAT
!            where Bn = C^T.B.C for n-th component
!  BVEC  ... misfit vector: BVEC = BVEC + CMAT_i^T.Bn.En
!            where En=C^T.B.E for n-th component
! constraint function exponent in source coordinates, X
! 1/(2 lambda).(X - BMAT^-1.BVEC)^T.BMAT.(X - BMAT^-1.BVEC)
! for lambda -> +INF
      type(TDMATRIX) :: BMAT
      type(TDARRAY) :: BVEC
      real(KIND(1.D0)) :: BCNT

! index map for non-integrable variables used by components
      type(TDARRAYI) :: IXRND

! TOF mode=true if there was a time limited component found on the beamline
      logical :: MX_TOF

! Matrices created in RESMAT_EVAL for further use :
!   COV    ... covariance matrix
!   UMAT   ... matrix diagonalizing COV (eigenvectors)
!   DIAG   ... vector with diagonalized COV (eigenvalues)
!   DETCOV ... determinant COV
      type(TDMATRIX) :: COV,UMAT
      type(TDARRAY) :: DIAG
      real(KIND(1.D0)) :: DETCOV
  ! legend for matrix elements
      character*128 MXLEG


  ! debug options
      logical :: DBG=.false.     ! flag
      integer :: IUDBG=6   ! Log-file unit

      interface DeleteIndex
        module procedure DeleteIndex_M
        module procedure DeleteIndex_V
      end interface DeleteIndex


      public s_tria,s_gaus,s_step,s_circ
      public MX_INIT,MX_DISPOSE,MX_ADD_TRMAT,MX_ADD_CONSTRAINT,MX_ADD_VARIABLES
      public MX_ADD_TRANSPORT,MX_EIGEN,MX_EVAL,MX_REVERSE
      public IXRND
      CONTAINS

!-----------------------------------------------------
! Initialize arrays for diension = N
! Basic variable space is (x,y,z,kx,ky,kz,t)
! All other variables (i>7) are specific to components
!-----------------------------------------------------
      subroutine MX_INIT(N)
      integer,intent(in) :: N
        call DYN_NEWMATRIX(TMAT,N,N)
        call DYN_NEWMATRIX(BMAT,N,N)
        call DYN_UNIMAT(CMAT,N)
        call DYN_NEWARRAY(TVEC,N)
        call DYN_NEWARRAY(BVEC,N)
        call DYN_NEWARRAY(IXRND,N)
        TCNT=0.D0
        BCNT=0.D0
        NTRM=0
        MXLEG=''
        MX_TOF=.false.
        call MX_ADD_VARIABLES('x|y|z|kx|ky|kz|t',(/1,2,0,4,5,6,0/),7)
      end subroutine MX_INIT

!-----------------------------------------------------
! Dispose dynamically allocated arrays
!-----------------------------------------------------
      subroutine MX_DISPOSE
        call DYN_DISPOSE(CMAT)
        call DYN_DISPOSE(TMAT)
        call DYN_DISPOSE(BMAT)
        call DYN_DISPOSE(TVEC)
        call DYN_DISPOSE(BVEC)
        call DYN_DISPOSE(IXRND)
      end subroutine MX_DISPOSE

!---------------------------------------------------------
      subroutine MX_ADD_TRMAT(C,T,D,INDX,NT)
! Add next transmission matrix
! Component arrays:
!  T  ... transmission matrix
!  C  ... coordinate conversion matrix:
!  D  ... misfit vector
!  DIMT ... leading dimension
!  NT  ... true dimension
!  INDX ... index array mapping T subspace on the global space of variables
! C is defined so that C.X+D  converts X from preceding component's exit coordinates to
! this component's local coordinates
!--------------------------------------------------------
      integer,intent(in) :: INDX(:),NT
      real(kind(1.D0)), intent(in) :: C(:,:),T(:,:),D(:)
      real(kind(1.D0)) :: Tn(NT,NT),Dn(NT)
! Tn = C^T.T.C
      call BTAB(T,C,Tn,NT)
! Update transmission arrays
!  TMAT  ... transmission matrix, TMAT = TMAT + CMAT^T.Tn.CMAT
      call DYN_AddBTAB(TMAT,CMAT,Tn,INDX,NT)
! Dn=C^T.T.D
      call BTAC(C,T,D,Dn,NT)
! update misfit term
!  TVEC  ... misfit vector: TVEC = TVEC + CMAT_i^T.Tn.Dn
      call DYN_AddBTAV(TVEC,CMAT,Tn,Dn,INDX,NT)
!  TCNT  ... constant term: TCNT = TCNT + D^T.T.D
      TCNT=TCNT+XTMX(D,T,NT)
      end subroutine MX_ADD_TRMAT

!---------------------------------------------------------
      subroutine MX_ADD_TRANSPORT(C,INDX,NC)
! Updates transport matrix CMAT for the next component exit coordinates
! Component arrays:
!  C  ... coordinate conversion matrix:
!  NC   ... variables subspace dimension for C and INDX
!  INDX ... index array mapping C subspace on the global space of variables
! C is defined so that C.X  converts X from preceding component's exit coordinates to
! this component's exit coordinates
!--------------------------------------------------------
      integer,intent(in) :: INDX(:),NC
      real(kind(1.D0)), intent(in) :: C(:,:)
! CMAT  ... transport matrix, CMAT = C.CMAT
      call DYN_MltpBySubmatrix(0,CMAT,C,INDX,NC)
      end subroutine MX_ADD_TRANSPORT

!---------------------------------------------------------
      subroutine MX_ADD_CONSTRAINT(C,B,D,INDX,NB)
! Add next transmission matrix
! Component arrays:
!  B  ... constraint matrix
!  C  ... coordinate conversion matrix:
!  D  ... misfit vector
!  DIMB ... leading dimension
!  NB  ... true dimension
!  INDX ... index array mapping B subspace on the global space of variables
! C is defined so that C.X+D  converts X from preceding component's exit coordinates to
! this component's local coordinates
!--------------------------------------------------------
      integer,intent(in) :: INDX(:),NB
      real(kind(1.D0)), intent(in) :: C(:,:),B(:,:),D(:)
      real(kind(1.D0)) :: Bn(NB,NB),Dn(NB)
! Bn = C^B.B.C
      call BTAB(B,C,Bn,NB)
! Update transmission arrays
!  TMAT  ... transmission matrix, TMAT = TMAT + CMAT^B.Bn.CMAT
      call DYN_AddBTAB(BMAT,CMAT,Bn,INDX,NB)
! Dn=C^B.B.D
      call BTAC(C,B,D,Dn,NB)
! update misfit term
!  TVEC  ... misfit vector: TVEC = TVEC + CMAT_i^B.Bn.Dn
      call DYN_AddBTAV(BVEC,CMAT,Bn,Dn,INDX,NB)
!  BCNT  ... constant term: BCNT = BCNT + D^T.B.D
      BCNT=BCNT+XTMX(D,B,NB)
      end subroutine MX_ADD_CONSTRAINT


!---------------------------------------------------------
      subroutine MX_ADD_VARIABLES(CAP,IX,NV)
! Add new variables
! new variables added by the component:
! CAP ... caption string (| separated list of variable names)
! IX  ... index vector, <>0 for non-integable variables
! NV  ... number of new variables
!--------------------------------------------------------
      integer,intent(in) :: IX(:),NV
      character(*) :: CAP
      integer :: i
      if (NV.gt.0) then
        do i=1,NV
          NTRM=NTRM+1
          if (IX(i).gt.0) IXRND%X(NTRM)=NTRM
        enddo
        MXLEG=trim(MXLEG)//'|'//trim(CAP)
      endif
      end subroutine MX_ADD_VARIABLES


!-----------------------------------------------------------------------
!C Evaluate matrices:
!C Integrate over variables other than x,y,z,kx,ky,kz
!C Calculate matrices required for further use (see file header)
!-----------------------------------------------------------------------
      subroutine MX_EVAL
      INTEGER :: i,j
      type(TDMATRIX) :: AUX,AUX1
      real(kind(1.D0)) :: DUM

! check consistency of NTRM and matrix dimension
      if (TMAT%NX.ne.NTRM) then
        write(*,*) 'ERROR in MX_EVAL: unequal number of rand. variables: NTRM=',NTRM,' TMAT%NX=',TMAT%NX
        return
      endif

      if (DBG) call XML_RSXDUMP(IUDBG,' ',1)

! Set TOF mode on if there is a time-limited component
      MX_TOF=(TMAT%M(7,7).gt.1.D-10)

! Constraints on z=0 at starting surface
! (either source or detector, depending on tracing direction).
      BMAT%M(3,3)=BMAT%M(3,3)+1.D0
! time is set to t=0 at the source
      if (.NOT.MX_TOF) BMAT%M(7,7)=BMAT%M(7,7)+1.D0

! only secondary section: integrate over x,y,.. kz,t
! to be moved to another procedure
!      if (TROPT%MODE.eq.tr_secondary) IXRND%X(1:NTRM_PRIM)=0

! first integrate over z
      call IntegrateVariable(3,IXRND)
! integrate all other constrained variables
      i=NTRM
      do while(i.gt.0)
        call IntegrateVariable(i,IXRND)
        i=i-1
      enddo

! check it !
      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CMAT integrated ',MXLEG,MXLEG,CMAT%M,CMAT%DX)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'TRM integrated ',MXLEG,MXLEG,TMAT%M,TMAT%DX)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'CON integrated ',MXLEG,MXLEG,BMAT%M,BMAT%DX)
      endif

! COV
      if (DBG) write(*,*) 'RESMAT_EVAL  INVERT start'
      call DYN_NEWMATRIX(COV,TMAT%DX,TMAT%DX)
      J=INVERT(NTRM,TMAT%M,TMAT%DX,COV%M,COV%DX)

      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV ',MXLEG,MXLEG,COV%M,COV%DX)
      endif
      if (DBG) write(*,*) 'RESMAT_EVAL  INVERT OK'

! determinants
      DETCOV=DYN_DETJAC(COV)
      if (DBG) write(*,*) 'RESMAT_EVAL  DET OK'

! DIAG,UMAT
      CALL DYN_JACOBI(COV,DIAG,UMAT,I)
      if (DBG) write(*,*) 'RESMAT_EVAL  JACOBI OK ',I

  ! check it !
      if (DBG) then
        call DYN_NEWMATRIX(AUX,COV%NX,COV%NX)
        call DYN_NEWMATRIX(AUX1,COV%NX,COV%NX)

        call DYN_MxM(1,COV,UMAT,AUX)
        call DYN_MxM(-1,UMAT,AUX,AUX1)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV ',MXLEG,MXLEG,COV%M,COV%DX)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'UT.COV.U ',MXLEG,MXLEG,AUX1%M,COV%DX)
        call DYN_DISPOSE(AUX)
        call DYN_DISPOSE(AUX1)
      endif

      if (DBG) then
        call XML_MATRIX(IUDBG,NTRM,NTRM,'COV',MXLEG,MXLEG,COV%M,COV%DX)
        call XML_MATRIX(IUDBG,NTRM,NTRM,'UMAT',MXLEG,MXLEG,UMAT%M,COV%DX)
        call XML_ARRAY(IUDBG,NTRM,'DIAG',MXLEG,DIAG%X)
        DUM=1.D0
        do i=1,NTRM
          DUM=DUM*DIAG%X(i)
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
      end subroutine MX_EVAL


!-----------------------------------------------------------------------
!C Reverse tracing direction to up-stream, i.e. invert transport matrix etc.
!-----------------------------------------------------------------------
      subroutine MX_REVERSE
      type(TDMATRIX) :: AUX
      integer :: I
    ! INVERT works with a local copy of CMAT,
    ! so we can use the same variable for output and input
      I=INVERT(NTRM,CMAT%M,CMAT%DX,CMAT%M,CMAT%DX)
    ! transform transmission matrix
      call DYN_NEWMATRIX(AUX,CMAT%DX,CMAT%DX)
      call BTAB(TMAT%M,CMAT%M,AUX%M,NTRM)
      TMAT%M=AUX%M
    ! transform constraint matrix
      call DYN_NEWMATRIX(AUX,CMAT%DX,CMAT%DX)
      call BTAB(BMAT%M,CMAT%M,AUX%M,NTRM)
      BMAT%M=AUX%M
      call DYN_UNIMAT(CMAT,CMAT%DX)
      call DYN_DISPOSE(AUX)
      end subroutine MX_REVERSE

!--------------------------------------------------------------
!C Return eigenvalues and eigenvectors of the covariance matrix
!C does not evaluate! RESMAT_EVAL must be called before
!C   EIVALS  ... vector diagonalized covariance matrix (eigenvalues)
!C   EIVECT  ... matrix diagonalizing the covariance matrix (eigenvectors)
!C   IDIR     ... transport direction
!C   NL,N ... leading dimension and rank of the matrices
!C NOTE: the matrix rank, N is for output !
!---------------------------------------------------------------
      subroutine MX_EIGEN(EIVECT,EIVALS,NL,N)
      integer, intent(in) :: NL
      REAL(KIND(1.D0)),intent(out) :: EIVALS(NL),EIVECT(NL,NL)
      integer, intent(out) :: N
      integer :: i,j
      do j=1,NTRM
        do i=1,NTRM
          EIVECT(i,j)=UMAT%M(i,j)
        enddo
        EIVALS(j)=DIAG%X(j)
      enddo
      N=NTRM
      end subroutine MX_EIGEN

!------------------------------------------------------------------------
! Integrate over the k-th variable, reducing matrices TRM,CMAT and CON
!----------------------------------------------------------------------
      subroutine IntegrateVariable(k,IDX)
      integer,intent(in) :: k
      type(TDARRAYI) :: IDX
      integer :: i,j,IS,IL,LL
      real(KIND(1.D0)),parameter :: EPS=1.D-12
      real(KIND(1.D0)) :: Z
      logical :: is_integrable,isConstrained
      character*256 TMPLEG
      character*16 :: VN
      ! get variable ID
        VN=' '
        call FINDSTRPAR(MXLEG,'|',k,IS,IL)
        if (IL.GT.0) VN=MXLEG(IS:IS+IL-1)
        is_integrable=.TRUE.
      ! check if the variable is integrable
        do j=1,IDX%NX
          is_integrable=(is_integrable.and.(IDX%X(j).ne.k))
        enddo
        if (is_integrable.and.(ABS(TMAT%M(k,k)).le.EPS)) then
          is_integrable=.false.
          write(*,*) 'ERROR:  ',trim(VN),' - can''t integrate variable, diag=',TMAT%M(k,k)
        endif
        if (is_integrable) then
          if (DBG) write(IUDBG,*) '   ',trim(VN),' is integrable '
        else
          if (DBG) write(IUDBG,*) '   ',trim(VN),' is not integrable '
          return
        endif

! remove legend item from the list
        LL=len_trim(MXLEG)
        if (IL.GT.0) then
          TMPLEG=MXLEG
          if (LL.gt.IS+IL) then
            if (IS.eq.1) then
              MXLEG=TMPLEG(IS+IL+1:LL)
            else
              MXLEG=TMPLEG(1:IS-2)//TMPLEG(IS+IL:LL)
            endif
          else
            MXLEG=TMPLEG(1:IS-2)
          endif
        endif

! check if the variable is constrained (BMAT diagobal element <> 0)
        isConstrained=(ABS(BMAT%M(k,k)).gt.0.D0)
        if (isConstrained) then
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) then
                  isConstrained = (abs(BMAT%M(i,k)*BMAT%M(k,j)).lt.BMAT%M(k,k)*1.D10)
                endif
                if (.not.isConstrained) exit
              enddo
            endif
            if (.not.isConstrained) exit
          enddo
        endif
! there is a non-zero constraint k-element
        if (isConstrained) then
          if (DBG) write(IUDBG,*) '   ',trim(VN),' is constrained: ',BMAT%M(k,k)
        ! new transmission matrix
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) then
                  Z=BMAT%M(i,j)*TMAT%M(k,k)-BMAT%M(k,j)*TMAT%M(i,k)-BMAT%M(i,k)*TMAT%M(k,j)
                  TMAT%M(i,j)=TMAT%M(i,j)+Z/BMAT%M(k,k)
                endif
              enddo
            endif
            TVEC%X(i)=TVEC%X(i) + (TMAT%M(k,k)*BVEC%X(i)-BMAT%M(i,k)*TVEC%X(k)-TMAT%M(i,k)*BVEC%X(k))/BMAT%M(k,k)
          enddo
          TCNT=TCNT+(TMAT%M(k,k)*BCNT-2*TVEC%X(k)*BVEC%X(k))/BMAT%M(k,k)
        ! new constraint matrix
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) BMAT%M(i,j)=BMAT%M(i,j)-BMAT%M(i,k)*BMAT%M(k,j)/BMAT%M(k,k)
              enddo
            endif
            BVEC%X(i)=BVEC%X(i) - BMAT%M(i,k)*BVEC%X(k)/BMAT%M(k,k)
          enddo
          BCNT=BCNT-BVEC%X(k)**2/BMAT%M(k,k)
! there is a no constraint for the k-element
        else
          if (DBG) write(IUDBG,*) '   ',trim(VN),' unconstrained, diag=',TMAT%M(k,k)
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) TMAT%M(i,j)=TMAT%M(i,j)-TMAT%M(i,k)*TMAT%M(k,j)/TMAT%M(k,k)
              enddo
            endif
            TVEC%X(i)=TVEC%X(i) - TMAT%M(i,k)*TVEC%X(k)/TMAT%M(k,k)
          enddo
          TCNT=TCNT-TVEC%X(k)**2/TMAT%M(k,k)
          do i=1,NTRM
            if (i.ne.k) then
              do j=1,NTRM
                if (j.ne.k) BMAT%M(i,j)=BMAT%M(i,j)-(BMAT%M(k,j)*TMAT%M(i,k)+BMAT%M(i,k)*TMAT%M(k,j))/TMAT%M(k,k)
              enddo
            endif
            BVEC%X(i)=BVEC%X(i) - (TMAT%M(i,k)*BVEC%X(k)+BMAT%M(i,k)*TVEC%X(k))/TMAT%M(k,k)
          enddo
          BCNT=BCNT-2*TVEC%X(k)*BVEC%X(k)/TMAT%M(k,k)
        endif
  ! reduce dimension
        call DeleteIndex(TMAT,k,k)
        call DeleteIndex(BMAT,k,k)
        call DeleteIndex(CMAT,k,k)
        call DeleteIndex(TVEC,k)
        call DeleteIndex(BVEC,k)
  ! update pointers to XRND array (used by components)
        do j=1,IDX%NX
          if (IDX%X(j).gt.k) IDX%X(j)=IDX%X(j)-1
        enddo
        NTRM=NTRM-1
        if (DBG) write(IUDBG,*) '   ','new dimension: ',NTRM
      end subroutine IntegrateVariable

!-----------------------------------------------------
      subroutine DeleteIndex_M(A,IROW,ICOL)
! Delete one row and column (i.e. compact and lower matrix rank by 1)
! IROW,ICOL are the indices of row and column to be deleted
!-----------------------------------------------------
      integer, intent(in) :: IROW,ICOL
      type(TDMATRIX) :: A
      integer :: i,j
  ! remove row
      do i=1,A%NY
        do j=IROW,A%NX-1
          A%M(j,i)=A%M(j+1,i)
        enddo
        A%M(A%NX,i)=0.D0
      enddo
  ! remove column
      do i=1,A%NX-1
        do j=ICOL,A%NY-1
          A%M(i,j)=A%M(i,j+1)
        enddo
        A%M(i,A%NY)=0.D0
      enddo
      A%NX=A%NX-1
      A%NY=A%NY-1
      end subroutine DeleteIndex_M

!-----------------------------------------------------
      subroutine DeleteIndex_V(A,IROW)
! Delete one row of a vector
!-----------------------------------------------------
      integer, intent(in) :: IROW
      type(TDARRAY) :: A
      integer :: i
  ! remove row
      do i=IROW,A%NX-1
        A%X(i)=A%X(i+1)
      enddo
      A%X(A%NX)=0.D0
      A%NX=A%NX-1
      end subroutine DeleteIndex_V

      end module MATRIX
