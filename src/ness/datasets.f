!C//////////////////////////////////////////////////////////////////////
!C////  $Id: datasets.f,v 1.15 2019/08/16 17:16:26 saroun Exp $
!C////
!C////  R E S T R A X - Simulation of neutron three-axis spectrometers
!C////
!C////     Copyright (C) 1995-2006, All rights reserved
!C////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!C////     Institut Laue Langevin, Grenoble, France
!C////
!C////     Written by:  Jan Saroun
!C////     $Revision: 1.15 $
!C////     $Date: 2019/08/16 17:16:26 $
!C//////////////////////////////////////////////////////////////////////
!C////
!C////  Hanles parameters and manipulation with data sets
!C////  Data set may represent a single channel in multichannel device or
!C////  a completely different configuration, but with the same set of components.
!C////  Subsets are defined for each dataset. They represent different scan positions.
!C//////////////////////////////////////////////////////////////////////
      module DATASETS
      use XMLINFO
      use NSTORE,ONLY: NSTORE_UNREGISTER,NSTORE_REGISTER
      use CLASSES,ONLY:BEAMLINE_DIM
      use CLASSDEF,ONLY:MAX_FIELDS
      IMPLICIT NONE

      private

      integer,parameter :: MDATS=128   ! maximum number of datasets
      integer,parameter :: MSETS=32   ! maximum number of subsets
      integer :: DSET_N            ! numer of defined datasets
      integer :: DSET_NS(MDATS)    ! numer of defined subsets in each dataset
      integer :: DSET_NP    ! numer of parameters for a single dataset
      integer :: DSET_ISEL  ! selected dataset
      ! maximum number of parameters per one beamline
      ! add 16000 to account for free components with many parameters (SGUIDES)
      integer,parameter :: DSET_DIM = 16000+2*MAX_FIELDS*BEAMLINE_DIM+MAX_FIELDS
      ! max. number of differences per one beamline
      integer,parameter :: DSET_MAXDIF = MAX_FIELDS*BEAMLINE_DIM
      ! reference dataset parameters
      real(KIND(1.D0)) :: DSET_REF(DSET_DIM)
      ! differences w.r.t. DSET_REF for all datasets
      real(KIND(1.D0)) :: DSET_DIF(MDATS*DSET_MAXDIF)
      ! pointers of values in DSET_DIF to corresponding fields in DSET_REF
      integer :: DSET_IDIF(MDATS*DSET_MAXDIF)
      ! partitioning of DSET_DIF on individual datastes
      integer :: DSET_NDIF(0:MDATS)

      TYPE TDATASET  ; SEQUENCE
        real(KIND(1.D0)) ::  KI,KF(3)     ! nominal setting in Lab coordinates
        real(KIND(1.D0)) ::  MCL(3,3)    ! trans. matrix from Lab to C&N
        real(KIND(1.D0)) ::  MFL(3,3)    ! trans. matrix from Lab to KF coordinates
        integer :: REG1,REG2  ! registry references for primary and secondary paths
      end type TDATASET

      TYPE PDATASET
        TYPE(TDATASET),POINTER :: P
      END TYPE PDATASET

! store nominal position of the instrument for each subset and dataset
      TYPE(PDATASET) :: PDSET(1:MDATS)  ! dataset position
      TYPE(PDATASET) :: PSSET(1:MSETS,1:MDATS)  ! subset positions

      public DSET_WRITE
      public DSET_READ
      public DSET_ADD
      public DSET_CLEAR
      public DSET_N
      public DSET_ISEL
      public MDATS,MSETS
      public TDATASET,PDATASET
      public DSET_SETPOS
      public DSET_DATA
      public DSET_MKREF
      public DSET_GETK0
      contains


!--------------------------------------------------------
      subroutine DSET_CLEAR
! clear all datasets
!--------------------------------------------------------
      integer :: i,j
      !write(*,*) 'DSET_CLEAR ',DSET_N,DSET_NS(DSET_N)
      do i=1,DSET_N
        do j=1,DSET_NS(i)
          if (associated(PSSET(i,j)%P)) DEALLOCATE(PSSET(i,j)%P)
        enddo
        if (associated(PDSET(i)%P)) then
          call NSTORE_UNREGISTER(PDSET(i)%P%REG1)
          call NSTORE_UNREGISTER(PDSET(i)%P%REG2)
          DEALLOCATE(PDSET(i)%P)
        endif
      enddo
      DSET_N=0
      DSET_NS=0
      DSET_NDIF=0
      DSET_NP=0
      DSET_ISEL=0

      end subroutine DSET_CLEAR

!--------------------------------------------------------
      TYPE(TDATASET) function DSET_GETPOS(IDAT)
! get dataset nominal position
!--------------------------------------------------------
      integer,intent(in) :: IDAT
      type(TDATASET) :: DT
  !    write(*,*) 'DSET_GETPOS ',IDAT,DSET_N
      if (check_index(IDAT)) then
        if (associated(PDSET(IDAT)%P)) then
          DT=PDSET(IDAT)%P
        else
          call MSG_ERROR('DSET_GETPOS','Required dataset is not allocated',0,1)
        endif
      endif
      DSET_GETPOS=DT
      end function DSET_GETPOS


!--------------------------------------------------------
      TYPE(PDATASET) function DSET_DATA(IDAT)
! get dataset as pointer
!--------------------------------------------------------
      integer,intent(in) :: IDAT
      type(PDATASET) :: DT
      if (check_index(IDAT)) then
        if (associated(PDSET(IDAT)%P)) then
          DT=PDSET(IDAT)
        else
          NULLIFY(DT%P)
          call MSG_ERROR('DSET_GETPOS','Required dataset is not allocated',0,1)
        endif
      endif
      DSET_DATA=DT
      end function DSET_DATA


!--------------------------------------------------------
      REAL(kind(1.D0)) function DSET_GETK0(IDAT,ISPEC)
! get nominal K0 for given dataset
! ISPEC=1 ... KI
! ISPEC=2 ... KF
!--------------------------------------------------------
      integer,intent(in) :: IDAT,ISPEC
      type(TDATASET) :: DT
      REAL(kind(1.D0)) :: K0
      K0=0.D0
  !    write(*,*) 'DSET_GETPOS ',IDAT,DSET_N
      if (check_index(IDAT)) then
        if (associated(PDSET(IDAT)%P)) then
          DT=PDSET(IDAT)%P
          select case (ISPEC)
          case(1)
            K0=DT%KI
          case default
            K0=SQRT(DT%KF(1)**2+DT%KF(2)**2+DT%KF(3)**2)
          end select
        else
          call MSG_ERROR('DSET_GETPOS','Required dataset is not allocated',0,1)
        endif
      endif
      DSET_GETK0=K0
      end function DSET_GETK0


!--------------------------------------------------------
      subroutine DSET_SETPOS(IDAT,DT)
! set dataset nominal position
!--------------------------------------------------------
      integer,intent(in) :: IDAT
      type(TDATASET),intent(in) :: DT
  !    write(*,*) 'DSET_SETPOS ',IDAT,DSET_N
      if (check_index(IDAT)) then
        if (associated(PDSET(IDAT)%P)) then
      ! copy fileds individually => preserve REG1 and REG2 references !!
          PDSET(IDAT)%P%KI=DT%KI
          PDSET(IDAT)%P%KF=DT%KF
          PDSET(IDAT)%P%MCL=DT%MCL
          PDSET(IDAT)%P%MFL=DT%MFL
        else
          call MSG_ERROR('DSET_SETPOS','Can''t access nonallocated dataset',0,1)
        endif
      endif
      end subroutine DSET_SETPOS

!--------------------------------------------------------
      subroutine DSET_ADD
! add a new dataset
!--------------------------------------------------------
      if (check_dmax(DSET_N+1)) then
        DSET_N=DSET_N+1
        if (associated(PDSET(DSET_N)%P)) then
          DEallocate(PDSET(DSET_N)%P)
          ! write(*,*) 'DSET_ADD deallocated ',DSET_N
        endif
        allocate(PDSET(DSET_N)%P)
        PDSET(DSET_N)%P%REG1=0
        PDSET(DSET_N)%P%REG2=0
        ! write(*,*) 'DSET_ADD allocated ',DSET_N
      endif
      end subroutine DSET_ADD

!--------------------------------------------------------
      subroutine DSET_DEL(IX)
! delete given dataset
!--------------------------------------------------------
      integer,intent(in) :: IX
      integer :: NP,I
  !    write(*,*) 'DSET_DEL ',IX
      if (check_index(IX)) then
        NP=DSET_NDIF(IX)-DSET_NDIF(IX-1)
        call DSET_SHRINK(IX,NP)
        do i=IX,MDATS-1
          DSET_NDIF(i)=DSET_NDIF(i+1)
        enddo
        DSET_NDIF(MDATS)=DSET_NDIF(MDATS-1)
        DEallocate(PDSET(DSET_N)%P)
        DSET_N=DSET_N-1
      endif
      end subroutine DSET_DEL

!--------------------------------------------------------
      subroutine DSET_MKREF(PAR,NPAR)
! define datasets reference
! clear all datasets
!--------------------------------------------------------
      real(KIND(1.D0)),intent(in) :: PAR(*)
      integer,intent(in) :: NPAR
      character*32 :: CNUM
      integer :: i
      if (NPAR>DSET_DIM) then
        call INT2STR(DSET_DIM,CNUM)
        ! write(*,*) 'DSET_MKREF','number of instrument parameters exceeds maximum '//trim(CNUM)
        call MSG_ERROR('DSET_MKREF','number of instrument parameters exceeds maximum '//trim(CNUM),1,1)
        !return
      endif
      call DSET_CLEAR
      DSET_NP=NPAR
      do i=1,DSET_NP
        DSET_REF(i)=PAR(i)
      enddo
      ! write(*,*) 'DSET_MKREF ',DSET_NP
      end subroutine DSET_MKREF


!--------------------------------------------------------
      integer function DSET_COUNTDIF(PAR,NPAR)
! count the number of differences w.r.t. DSET_REF
!--------------------------------------------------------
      real(KIND(1.D0)),intent(in) :: PAR(*)
      integer,intent(in) :: NPAR
      integer :: N,i
      N=0
      if (check_npar(npar)) then
        do i=1,DSET_NP
          if (DSET_REF(i).ne.PAR(i)) then
            N=N+1
          endif
        enddo
      endif
      DSET_COUNTDIF=N
      end function DSET_COUNTDIF

!--------------------------------------------------------
      subroutine DSET_READ(IX,PAR,NPAR)
! read dataset parameters from PAR to DSET_DIF
!--------------------------------------------------------
      real(KIND(1.D0)),intent(in) :: PAR(*)
      integer,intent(in) :: IX,NPAR
      integer :: i,ND

  !    write(*,*) 'DSET_READ ',IX,DSET_N,NPAR,DSET_NP
      if (.not.check_npar(NPAR)) return
  ! new number of fields in DSET_DIF array for IX-th dataset
      ND=DSET_COUNTDIF(PAR,NPAR)
  ! expand array
      if (ND.gt.(DSET_NDIF(IX)-DSET_NDIF(IX-1))) then
        call DSET_EXPAND(IX,ND-DSET_NDIF(IX)+DSET_NDIF(IX-1))
  ! shrink array
      else if (ND.lt.(DSET_NDIF(IX)-DSET_NDIF(IX-1))) then
        call DSET_SHRINK(IX,DSET_NDIF(IX)-DSET_NDIF(IX-1)-ND)
      endif
  ! store differences in DSET_DIF
      DSET_NDIF(IX)=DSET_NDIF(IX-1)
      do i=1,DSET_NP
        if (DSET_REF(i).ne.PAR(i)) then
          DSET_NDIF(IX)=DSET_NDIF(IX)+1
          DSET_DIF(DSET_NDIF(IX))=PAR(i)
          DSET_IDIF(DSET_NDIF(IX))=i
        endif
      endDO
      end subroutine DSET_READ

!--------------------------------------------------------
      subroutine DSET_EXPAND(IX,NS)
! expand IX-th field in the DIF array by NS items
!--------------------------------------------------------
      integer,intent(in) :: IX,NS
      integer :: i
  !    write(*,*) 'DSET_EXPAND ',IX,NS
      if (check_index(IX)) then
        if (check_maxdif(DSET_NDIF(MDATS)+NS)) then
          i=DSET_NDIF(MDATS)
          do while(i.gt.DSET_NDIF(IX))
            DSET_DIF(i+NS)=DSET_DIF(i)
            DSET_IDIF(i+NS)=DSET_IDIF(i)
          enddo
          do i=IX,MDATS
            DSET_NDIF(i)=DSET_NDIF(i)+NS
          enddo
        endif
      endif
      end subroutine DSET_EXPAND

!--------------------------------------------------------
      subroutine DSET_SHRINK(IX,NS)
! shrink IX-th field in the DIF array by NS items
!--------------------------------------------------------
      integer,intent(in) :: IX,NS
      integer :: i
  !    write(*,*) 'DSET_SHRINK ',IX,NS
      if (check_index(IX)) then
        do i=DSET_NDIF(IX-1)+1,DSET_NDIF(MDATS)
          DSET_DIF(i)=DSET_DIF(i+NS)
          DSET_IDIF(i)=DSET_IDIF(i+NS)
        enddo
        do i=IX,MDATS
          DSET_NDIF(i)=DSET_NDIF(i)-NS
        enddo
      endif
      end subroutine DSET_SHRINK

!--------------------------------------------------------
      logical function DSET_WRITE(IX,PAR,NPAR)
! write dataset parameters to PAR
!--------------------------------------------------------
      real(KIND(1.D0)),intent(out) :: PAR(*)
      integer,intent(in) :: IX,NPAR
      integer :: j,i,ND
  !    write(*,*) 'DSET_WRITE ',IX,NPAR
      if (check_npar(NPAR).and.check_index(IX)) then
        j=0
        ND=DSET_NDIF(IX)-DSET_NDIF(IX-1) ! nubr of differences
        do i=1,DSET_NP
          if ((j.lt.ND).and.(DSET_IDIF(j+1).eq.i)) then
            j=j+1
            PAR(i)=DSET_DIF(DSET_NDIF(IX-1)+j)
          else
            PAR(i)=DSET_REF(i)
          endif
        endDO
        DSET_WRITE=.true.
      else
        DSET_WRITE=.false.
      endif
      end function DSET_WRITE

! ****  CHECK FUNCTIONS WITH ERROR MESSAGES ****

!--------------------------------------------------
      logical function check_maxdif(n)
! check dimension of DSET_DIF array
!--------------------------------------------------
      integer,intent(in) :: N
      if (n.lt.MDATS*DSET_MAXDIF) then
        check_maxdif=.true.
      else
        check_maxdif=.false.
        call MSG_ERROR('check_maxdif','too many diferences w.r.t. the reference dataset',1,1)
      endif
      end function check_maxdif

!--------------------------------------------------
      logical function check_npar(n)
! check number of parameters
!--------------------------------------------------
      integer,intent(in) :: N
      character(32) :: CNUM1,CNUM2
      if (n.eq.DSET_NP) then
        check_npar=.true.
      else
        check_npar=.false.
        call INT2STR(n,CNUM1)
        call INT2STR(DSET_NP,CNUM2)
        call MSG_ERROR('check_npar','unequal number of parameters, REQ='//trim(CNUM1)//' DEF='//trim(CNUM2),1,1)
      endif
      end function check_npar

!--------------------------------------------------
      logical function check_dmax(n)
! check max. number of datasets
!--------------------------------------------------
      integer,intent(in) :: N
      if (n.le.MDATS) then
        check_dmax=.true.
      else
        check_dmax=.false.
        call MSG_ERROR('check_dmax','reached maximum number of datasets',1,1)
      endif
      end function check_dmax

!--------------------------------------------------
      logical function check_index(n)
! check max. number of datasets
!--------------------------------------------------
      integer,intent(in) :: N
      if ((n.le.DSET_N).and.(n.gt.0)) then
        check_index=.true.
      else
        check_index=.false.
        call MSG_ERROR('check_index','required dataset is not defined',1,1)
      endif
      end function check_index

      end module DATASETS