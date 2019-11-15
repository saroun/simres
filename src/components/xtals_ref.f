!//////////////////////////////////////////////////////////////////////
!////  $Id: xtals_ref.f,v 1.27 2017/04/11 10:52:41 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.27 $
!////     $Date: 2017/04/11 10:52:41 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Stores data for Bragg reflections used by XTAL components
!////
!////////////////////////////////////////////////////////////////////////
      MODULE XTALS_REF
      use CONSTANTS
      use LATTICE
      use TABLE_CRYSTALS
      use XMLINFO
      use SIMATH
      use FIELDDEF
      use FIELDDATA
      implicit none

      integer,parameter :: MAX_CR = 32
      integer,parameter :: MAX_NODES=256
      integer,parameter :: MAX_XREF = 1024
! reciprocal lattice  node
      type RLNODE; SEQUENCE
        integer :: NP             ! number of subsequent nodes (NEXT and PATH elements)
        integer :: LEVEL          ! level on the reflections sequence from (0,0,0), =1 for primary reflections
        integer :: HKL(3)         ! hkl of this node
        real(kind(1.D0)) :: LAMBDA ! nominal wavelengths needed to reach this node
        integer :: ORD(MAX_NODES) ! order index (for ordering by DHKL) of secondary nodes
        integer :: NEXT(MAX_NODES)  ! pointers to nodes in RLNODE_TAB achievable by secondary reflections
        integer :: PATH(MAX_NODES)  ! pointers to reflections in XREF_TAB leading to the NEXT nodes
        real(kind(1.D0)) :: LAM(MAX_NODES) ! wavelengths for the reflections
        real(kind(1.D0)) :: DHKL(MAX_NODES) ! dhkl for the reflections
      end type RLNODE

! encapsulated pointer to RLNODE
      type PRLNODE
        TYPE(RLNODE),pointer :: N
      end TYPE PRLNODE
      TYPE(PRLNODE) :: RLNODE_TAB(MAX_NODES,MAX_CR)
      integer :: RLNODE_N(MAX_CR) =0
      integer :: RLNODE_ORD(MAX_NODES,MAX_CR) =0

! data type describing Bragg reflection
      TYPE XREF; SEQUENCE
        INTEGER :: HKL(3)  ! Miller indices
        real(kind(1.D0)) :: DHKL    ! dhkl [A]
        real(kind(1.D0)) :: G(3)    ! G-vector local coordinates
        real(kind(1.D0)) :: RHO(2)  ! curvature of the neutral surface
        real(kind(1.D0)) :: GAMMA(3)    ! vector perpendicular to G in (K,G) plane, normalized
        real(kind(1.D0)) :: DG_DR(3,3)    ! G-vector gradient in local coordinates
        real(kind(1.D0)) :: QML     ! Maier-Leibnitz reflectivity, 4*PI*(F*dhkl/V0)**2 [ A^-1 cm^-1]
        real(kind(1.D0)) :: GRADKK  ! gradient of delta_k/k [mm-1]
      end TYPE XREF

! encapsulated pointer to XREF
      type PXREF
        TYPE(XREF),pointer :: X
      end TYPE PXREF
      TYPE(PXREF) :: XREF_TAB(MAX_XREF,MAX_CR)
      integer :: XREF_N(MAX_CR) =0
! counter of crystals with associaed table items.
! ensures that each crystal instance uses unique table fields
      integer :: XREF_NCR=0

      contains

!---------------------------------------------------------
      subroutine XTALS_REF_DISPOSE
!---------------------------------------------------------
      call XREF_CLEAR_TABLE
      call RLNODE_CLEAR_TABLE
      XREF_NCR=0
      end subroutine XTALS_REF_DISPOSE



!---------------------------------------------------------
      subroutine XREF_REPORT(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      integer :: INODE,I
1     format('ipath   reflection --> node  dhkl  lambda QML theta -grad_dkk')
      if ((ICR.le.0).or.(ICR.gt.MAX_CR)) return
      write(*,1)
      do INODE=1,RLNODE_N(ICR)
        I=RLNODE_ORD(INODE,ICR)
        call XREF_REPORT_NODE(ICR,I)
      enddo
      end subroutine XREF_REPORT


!---------------------------------------------------------
      subroutine RLNODE_REPORT(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      integer :: INODE,I
      character(128) :: CNUM,CNUM1
1     format('inode  ord  level (h k l) lambda')
2     format(I3,2x,I3,2x,I3,2x,'(',a,') ',a)
      if ((ICR.le.0).or.(ICR.gt.MAX_CR)) return
      write(*,1)
      do INODE=1,RLNODE_N(ICR)
        I=RLNODE_ORD(INODE,ICR)
        if ((I.gt.0).and.(I.le.RLNODE_N(ICR))) then
          call IARRAY2STR(RLNODE_TAB(I,ICR)%N%HKL,3,CNUM)
          call FLOAT2STR(RLNODE_TAB(I,ICR)%N%LAMBDA,CNUM1)
          write(*,2) INODE,I,RLNODE_TAB(I,ICR)%N%LEVEL,trim(CNUM),trim(CNUM1)
        endif
      enddo
      end subroutine RLNODE_REPORT


!---------------------------------------------------------
      subroutine XREF_REPORT_NODE(ICR,INODE)
!---------------------------------------------------------
      integer,intent(in) :: ICR,INODE
      integer :: REF(3),I,J,IPATH
      real(kind(1.D0)) :: QML,LAM,DHKL,GRDKK
      TYPE(RLNODE) :: NODE

1     format(I3,') ',3(1x,I3),' --> ',3(1x,I3),5(2x,G10.4))
      if ((ICR.le.0).or.(ICR.gt.MAX_CR)) return
      if ((INODE.le.0).or.(INODE.gt.RLNODE_N(ICR))) return
      NODE=RLNODE_TAB(INODE,ICR)%N
      if (NODE%NP.gt.0) write(*,*) '(',NODE%HKL,') np=',NODE%NP
      do I=1,NODE%NP
        J=NODE%ORD(I)
        LAM=NODE%LAM(J)
        IPATH=NODE%PATH(J)
        if ((IPATH.gt.0).and.(IPATH.le.XREF_N(ICR))) then
          REF=XREF_TAB(IPATH,ICR)%X%HKL
          QML=XREF_TAB(IPATH,ICR)%X%QML
          DHKL=XREF_TAB(IPATH,ICR)%X%DHKL
          GRDKK=XREF_TAB(IPATH,ICR)%X%GRADKK
          write(*,1) IPATH,REF,NODE%HKL+REF,DHKL,LAM,QML,ASIN(LAM/2/DHKL)*180/PI,GRDKK
        else
          write(*,*) IPATH,') out of range'
        endif
      enddo
      end subroutine XREF_REPORT_NODE

!---------------------------------------------------------
      subroutine XREF_LIST_NODES(ICR,SELHKL,ARG,IERR)
! convert list of nodes to SELECT data type
! SELHKL(3) is the selected node hkl vector
!---------------------------------------------------------
      integer,intent(in) :: ICR
      integer,intent(in) :: SELHKL(3)
      TYPE(TFIELD) :: ARG
      integer,intent(out) :: ierr
      character(128) :: CNUM
      integer :: i,j,isel
        ierr=1
        if (FIELD_IS_TYPE(ARG,ID_FIELD_SELECT)) then
          if (ALLOC_FIELD(ARG,RLNODE_N(ICR))) then
            isel=RLNODE_GET(ICR,SELHKL)
            ARG%EORD=-1
            do i=1,RLNODE_N(ICR)
              j=RLNODE_ORD(I,ICR)
              if ((j.gt.0).and.(j.le.RLNODE_N(ICR))) then
                call IARRAY2STR(RLNODE_TAB(j,ICR)%N%HKL,3,CNUM)
                call ITEM2FIELD(ARG,i,trim(CNUM))
                ARG%NP=ARG%NP+1
                if (isel.eq.j) ARG%EORD=i-1
              else
                write(*,*) 'WARNING in XREF_LIST_NODES: wrong sort index,icr,j= ',icr,j,', max=',RLNODE_N(ICR)
              endif
            enddo
            write(*,*) 'XREF_LIST_NODES sel=',SELHKL,' isel=',isel,' ord=',ARG%EORD
          !  call RLNODE_REPORT(ICR)
            ierr=0
          endif
        endif
      end subroutine XREF_LIST_NODES


!---------------------------------------------------------
      integer function XREF_GET(ICR,HKL)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      integer,intent(in) :: HKL(3)
      integer :: i,j
      integer :: Z
      XREF_GET=0
      if ((ICR.gt.0).and.(ICR.le.MAX_CR)) then
        do i=1,XREF_N(ICR)
          Z=0
          do j=1,3
            Z=Z+ABS(XREF_TAB(i,ICR)%X%HKL(j)-HKL(j))
          enddo
          if (Z.eq.0) then
            XREF_GET=i
            exit
          endif
        enddo
      endif
      end function XREF_GET

!---------------------------------------------------------
      integer function XREF_GET_LAST(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      XREF_GET_LAST=0
      if ((ICR.gt.0).and.(ICR.le.MAX_CR)) XREF_GET_LAST=XREF_N(ICR)
      end function XREF_GET_LAST

!---------------------------------------------------------
      subroutine XREF_ADD(X,ICR,IERR)
!---------------------------------------------------------
      TYPE(XREF),intent(in) :: X
      integer,intent(in) :: ICR
      integer,intent(out) :: IERR
      character(32) :: S
      IERR=1
      if ((ICR.le.0).or.(ICR.gt.MAX_CR)) then
        call INT2STR(ICR,S)
        call MSG_ERROR('XREF_ADD','Crystal index out of range '//trim(S),0,1)
      else
        if (XREF_N(ICR).lt.MAX_XREF) then
          if ( XREF_ALLOC(XREF_TAB(XREF_N(ICR)+1,ICR)) ) then
            XREF_N(ICR)=XREF_N(ICR)+1
            XREF_TAB(XREF_N(ICR),ICR)%X=X
            IERR=0
          else
            IERR=4
          endif
        else
          call INT2STR(MAX_XREF+1,S)
          call MSG_ERROR('XREF_ADD','Too may reflections ('//trim(S)//')',0,1)
          IERR=3
        endif
      endif
      end subroutine XREF_ADD


!---------------------------------------------------------
      subroutine XREF_CLEAR_TABLE
!---------------------------------------------------------
      integer :: i
      do i=1,MAX_CR
        call XREF_CLEAR_CR(i)
      enddo
      end subroutine XREF_CLEAR_TABLE

!---------------------------------------------------------
      integer function XREF_ADD_CR()
! increment counter of crystal items and return its value
!---------------------------------------------------------
      if (XREF_NCR.LT.MAX_CR) then
        XREF_NCR=XREF_NCR+1
      endif
      XREF_ADD_CR=XREF_NCR
      end function XREF_ADD_CR

!---------------------------------------------------------
      subroutine XREF_CLEAR_CR(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      if ((ICR.gt.0).and.(ICR.LE.MAX_CR)) then
        do while (XREF_N(ICR).gt.0)
          call XREF_DISPOSE(XREF_TAB(XREF_N(ICR),ICR))
          XREF_N(ICR)=XREF_N(ICR)-1
        enddo
      endIF
      end subroutine XREF_CLEAR_CR

!---------------------------------------------------------
      subroutine XREF_CLEAR(X)
!---------------------------------------------------------
      TYPE(XREF),intent(out) :: X
        X%HKL=0
        X%G=0.D0
        X%RHO=0.D0
        X%DHKL=0.D0
        X%DG_DR=0.D0
        X%QML=0.D0
        X%GRADKK=0.D0
      end subroutine XREF_CLEAR

!---------------------------------------------------------
      logical function XREF_ALLOC(X)
!---------------------------------------------------------
      TYPE(PXREF) :: X
      integer :: ierr
        allocate(X%X,STAT=ierr)
        if(ierr.eq.0) then
          call XREF_CLEAR(X%X)
        endif
        XREF_ALLOC=(ierr.eq.0)
      end function XREF_ALLOC

!---------------------------------------------------------
      subroutine XREF_DISPOSE(X)
!---------------------------------------------------------
      TYPE(PXREF) :: X
        if (associated(X%X)) then
          DEallocate(X%X)
          nullify(X%X)
        endif
      end subroutine XREF_DISPOSE

!---------------------------------------------------------
      subroutine XREF_GET_ITEM(ICR,IREF,XR)
! return reflection data as PXREF object for given crystal and reflection index
!---------------------------------------------------------
      integer,intent(in) :: ICR,IREF
      TYPE(PXREF),intent(out) :: XR
        if ((IREF.gt.0).and.(IREF.le.XREF_GET_LAST(ICR))) then
          XR=XREF_TAB(IREF,ICR)
        else
          nullify(XR%X)
        endif
      end subroutine XREF_GET_ITEM

!=================================================================================================

!---------------------------------------------------------
      subroutine RLNODE_CLEAR_TABLE
!---------------------------------------------------------
      integer :: i
      do i=1,MAX_CR
        call RLNODE_CLEAR_CR(i)
      enddo
      end subroutine RLNODE_CLEAR_TABLE


!---------------------------------------------------------
      subroutine RLNODE_CLEAR_CR(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      if ((ICR.gt.0).and.(ICR.LE.MAX_CR)) then
        do while (RLNODE_N(ICR).gt.0)
          call RLNODE_DISPOSE(RLNODE_TAB(RLNODE_N(ICR),ICR))
          RLNODE_N(ICR)=RLNODE_N(ICR)-1
        enddo
      endIF
      end subroutine RLNODE_CLEAR_CR

!---------------------------------------------------------
      logical function RLNODE_ALLOC(N)
!---------------------------------------------------------
      TYPE(PRLNODE) :: N
      integer :: ierr
        allocate(N%N,STAT=ierr)
        if(ierr.eq.0) then
          call RLNODE_CLEAR(N%N)
        endif
        RLNODE_ALLOC=(ierr.eq.0)
      end function RLNODE_ALLOC

!---------------------------------------------------------
      subroutine RLNODE_DISPOSE(N)
!---------------------------------------------------------
      TYPE(PRLNODE) :: N
        if (associated(N%N)) then
          DEallocate(N%N)
          nullify(N%N)
        endif
      end subroutine RLNODE_DISPOSE

!---------------------------------------------------------
      subroutine RLNODE_CLEAR(N)
!---------------------------------------------------------
      TYPE(RLNODE),intent(out) :: N
        N%LEVEL=0
        N%HKL=0
        N%NEXT=0
        N%PATH=0
        N%LAM=0.D0
        N%DHKL=0.D0
        N%NP=0
        N%LAMBDA=0.D0
      end subroutine RLNODE_CLEAR

!-------------------------------------------------------------
      integer function RLNODE_GET(ICR,HKL)
! get index of the node with given HKL = pointer to CR%NODES
! return 1 if the next node is (0,0,0)
! return 0 if not found
!-------------------------------------------------------------
      integer,intent(in) :: ICR
      integer,intent(in) :: HKL(3)
      integer :: i,j
      integer :: Z
      RLNODE_GET=0
      if ((ICR.gt.0).and.(ICR.le.MAX_CR)) then
        do i=1,RLNODE_N(ICR)
          Z=0
          do j=1,3
            Z=Z+ABS(RLNODE_TAB(i,ICR)%N%HKL(j)-HKL(j))
          enddo
          if (Z.eq.0) then
            RLNODE_GET=i
            exit
          endif
        endDO
      endif
      end function RLNODE_GET

!---------------------------------------------------------
      integer function RLNODE_GET_LAST(ICR)
!---------------------------------------------------------
      integer,intent(in) :: ICR
      RLNODE_GET_LAST=0
      if ((ICR.gt.0).and.(ICR.le.MAX_CR)) RLNODE_GET_LAST=RLNODE_N(ICR)
      end function RLNODE_GET_LAST

!---------------------------------------------------------
      subroutine RLNODE_GET_ITEM(ICR,INODE,PNODE)
! return reflection data as PXREF object for given crystal and reflection index
!---------------------------------------------------------
      integer,intent(in) :: INODE,ICR
      TYPE(PRLNODE),intent(out) :: PNODE
        if ((INODE.gt.0).and.(INODE.le.RLNODE_GET_LAST(ICR))) then
          PNODE=RLNODE_TAB(INODE,ICR)
        else
          nullify(PNODE%N)
        endif
      end subroutine RLNODE_GET_ITEM

!---------------------------------------------------------
      subroutine RLNODE_ADD(N,ICR,IERR)
!---------------------------------------------------------
      TYPE(RLNODE),intent(in) :: N
      integer,intent(in) :: ICR
      integer,intent(out) :: IERR
      character(32) :: S
      IERR=1
      if ((ICR.le.0).or.(ICR.gt.MAX_CR)) then
        call INT2STR(ICR,S)
        call MSG_ERROR('RLNODE_ADD','Crystal index out of range '//trim(S),0,1)
      else
        if (RLNODE_N(ICR).lt.MAX_NODES) then
          if ( RLNODE_ALLOC(RLNODE_TAB(RLNODE_N(ICR)+1,ICR)) ) then
            RLNODE_N(ICR)=RLNODE_N(ICR)+1
            RLNODE_TAB(RLNODE_N(ICR),ICR)%N=N
            IERR=0
          else
            IERR=3
            call INT2STR(RLNODE_N(ICR),S)
            call MSG_ERROR('RLNODE_ADD','Allocation failed for node ('//trim(S)//')',0,1)
          endif
        else
          IERR=2
          call INT2STR(MAX_NODES+1,S)
          call MSG_ERROR('RLNODE_ADD','Too may reflections ('//trim(S)//')',0,1)
        endif
      endif
      end subroutine RLNODE_ADD

!---------------------------------------------------------
      integer function RLNODE_ASSIGN_XREF(NODE,ICR,XR,K0,IERR)
! assign new child reflection to a node
!---------------------------------------------------------
      TYPE(RLNODE),intent(inout) :: NODE
      integer,intent(in) :: ICR
      type(XREF),intent(in) :: XR
      real(kind(1.D0)),intent(in) :: K0
      integer,intent(out) :: IERR
      integer :: HKL(3),INEXT,IREF
      TYPE(RLNODE) :: NEWNODE
      character(32) :: S
      IERR=1
      RLNODE_ASSIGN_XREF=0
      if (NODE%NP.ge.MAX_NODES) then
        call INT2STR(MAX_NODES+1,S)
      !  call MSG_ERROR('RLNODE_ASSIGN_XREF','Too may reflections ('//trim(S)//')',0,1)
        return
      endif
! coordinates of the next node
      HKL=NODE%HKL+XR%HKL
! reflection leading to the next node
      IREF=XREF_GET(ICR,XR%HKL)
    ! add a new reflection if not yet defined
      if (IREF.le.0) then
        call XREF_ADD(XR,ICR,IERR)
        if (IERR.ne.0) return
        IREF=XREF_GET_LAST(ICR)
      !  write(*,*) '   RLNODE_ASSIGN_XREF new ref: hkl=',XR%HKL,' iref=',IREF
      endif
! identify the next node
      INEXT=RLNODE_GET(ICR,HKL)
    ! define a new node if not found
      if (INEXT.le.0) then
        call RLNODE_CLEAR(NEWNODE)
        NEWNODE%HKL=HKL
        NEWNODE%LEVEL=NODE%LEVEL+1
        NEWNODE%LAMBDA=2.D0*PI/K0
        call RLNODE_ADD(NEWNODE,ICR,IERR)
        if (IERR.ne.0) return
        INEXT=RLNODE_GET_LAST(ICR)
      !  write(*,*) '   RLNODE_ASSIGN_XREF new node: hkl=',HKL,' inode=',INEXT
      ! else
      !  write(*,*) '   added  old: hkl=',HKL,' idx=',INEXT
      endif
      NODE%NP=NODE%NP+1
      NODE%PATH(NODE%NP)=IREF
      NODE%NEXT(NODE%NP)=INEXT
      NODE%LAM(NODE%NP)=2.D0*PI/K0
      NODE%DHKL(NODE%NP)=XR%DHKL
      IERR=0
      RLNODE_ASSIGN_XREF=IREF
      end function RLNODE_ASSIGN_XREF

!---------------------------------------------------------
      subroutine RLNODE_SORT(NODE)
! sort by DHKL (descending) using Heapsort algorithm
!---------------------------------------------------------
      TYPE(RLNODE),intent(inout) :: NODE
        call HEAPSORT(NODE%DHKL(1),NODE%NP,-1,NODE%ORD)
      end subroutine RLNODE_SORT


      end MODULE XTALS_REF
