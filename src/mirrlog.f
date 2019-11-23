!//////////////////////////////////////////////////////////////////////
!////  $Id: mirrlog.f,v 1.4 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.4 $
!////     $Date: 2019/08/15 15:02:06 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Logger for mirror reflections
!////
!////////////////////////////////////////////////////////////////////////
      MODULE MIRRLOG
      use CONSTANTS
      use IO
      implicit none
      private

      integer, parameter :: mlog_none=0
      integer, parameter :: mlog_left=1
      integer, parameter :: mlog_right=2
      integer, parameter :: mlog_top=3
      integer, parameter :: mlog_bottom=4
      integer, parameter :: mlog_leftright=5
      integer, parameter :: mlog_topbottom=6
      integer, parameter :: mlog_all=7

! parameters of the bounce recorder (records distribution of m-values along z)
      integer, parameter :: MLOGS_DIM=128
      integer, parameter :: MLOGS_MAXZ=100
      integer, parameter :: MLOGS_MAXM=60
      real(kind(1.D0)) :: MLOGS_MAXVAL=6.D0
! pointer type for mirror logger
      type PMLOGS
        character(16) :: ID
        integer :: CCLS,INST  ! instance index
        ! record bounces for walls (left,right,top,bootom)
        integer :: NBNC,LOGBNC
        real(kind(1.D0)) :: ZRANGE
        real(kind(1.D0)), allocatable :: BNC(:,:)
        real(kind(1.D0)), allocatable :: TMPBNC(:,:)
        real(kind(1.D0)), allocatable :: TMPCOORD(:,:)
      end type PMLOGS

      integer :: MLOGS_NC
      type(PMLOGS) :: MLOGS(1:MLOGS_DIM)


      integer,private :: idbg
      logical,private :: dbg

      public MLOGS_FREEALL,MLOGS_REGISTER,MLOGS_REC_CLR,MLOGS_NC
      public MLOGS_REC_TMP,MLOGS_REC_UPD,MLOGS_REC_OUT,MLOGS_RESET,MLOGS_GET_REC
      public PMLOGS,MLOGS,MLOGS_PRINT_TMP


      contains



!---------------------------------------------------------
      integer function MLOGS_GET_REC(CCLS,INST)
! get record index
!---------------------------------------------------------
      integer,intent(in) :: CCLS,INST
      integer :: i,ires
      ires=0
      write(*,*) 'MLOGS_GET_REC ',MLOGS_NC
      i=0
      do while ((i<MLOGS_NC).and.(ires==0))
        i=i+1
         write(*,*) 'MLOGS_GET_REC ',i,MLOGS_NC,MLOGS(i)%CCLS,MLOGS(i)%INST
        if ((MLOGS(i)%CCLS==CCLS).and.(MLOGS(i)%INST==INST)) then
          if (MLOGS_isValid(i)) ires=i
        endif
      enddo
      MLOGS_GET_REC=ires
      end function MLOGS_GET_REC

!---------------------------------------------------------
      subroutine MLOGS_RESET(IREC)
! reset temporary logger for a next history
!---------------------------------------------------------
      integer,intent(in) :: IREC
      if ((IREC>0).and.(IREC<MLOGS_DIM)) then
        MLOGS(IREC)%NBNC=0
      endif
      end subroutine MLOGS_RESET

!---------------------------------------------------------
      subroutine MLOGS_FREEALL
! free all loggers
!---------------------------------------------------------
      integer :: i
      do i=1,MLOGS_NC
        call MLOGS_REC_FREE(i)
      enddo
      MLOGS_NC=0
      idbg=0
      dbg=.false.
      end subroutine MLOGS_FREEALL

!---------------------------------------------------------
      integer function MLOGS_REGISTER(ID,CCLS,INST,LOGBNC,ZRANGE)
! register new logger
!---------------------------------------------------------
      integer,intent(in) :: CCLS,INST,LOGBNC
      real(kind(1.D0)),intent(in) :: ZRANGE
      character(*) :: ID
      integer :: i,ires
      ires=0
      if (MLOGS_NC<MLOGS_DIM) then
        MLOGS_NC=MLOGS_NC+1
        i=MLOGS_NC
        MLOGS(i)%CCLS=CCLS
        MLOGS(i)%INST=INST
        MLOGS(i)%ID=trim(ID)
        MLOGS(i)%NBNC=0
        MLOGS(i)%ZRANGE=ZRANGE
        MLOGS(i)%LOGBNC=LOGBNC
        call MLOGS_REC_ALLOC(i)
        if (MLOGS_isValid(i))  ires=i
      endif
      MLOGS_REGISTER=ires
      end function MLOGS_REGISTER


!---------------------------------------------------------
      logical function MLOGS_isValid(IREC)
! clear bounce logger
!---------------------------------------------------------
      integer,intent(in) :: IREC
      logical :: res
      res = ((IREC>0).and.(IREC<=MLOGS_NC))
      if (res) res = allocated(MLOGS(IREC)%BNC)
      MLOGS_isValid=res
      end function MLOGS_isValid

!---------------------------------------------------------
      subroutine MLOGS_REC_CLR(IREC)
! clear bounce logger
!---------------------------------------------------------
      integer,intent(in) :: IREC
      integer :: i,j
      !return
      if (.not.MLOGS_isValid(IREC)) return
        !write(*,*) 'MLOGS_REC_CLR '//trim(MLOGS(IREC)%X%FRAME%ID),' ',IREC
      do j=1,MLOGS_MAXZ
        DO i=1,4*MLOGS_MAXM
          MLOGS(IREC)%BNC(i,j)=0.D0
        enddo
      enddo
      MLOGS(IREC)%NBNC=0
      end subroutine MLOGS_REC_CLR

!---------------------------------------------------------
      subroutine MLOGS_REC_ALLOC(IREC)
! allocate bounce logger
!---------------------------------------------------------
      integer,intent(in) :: IREC
      integer :: ierr
      if (.NOT.allocated (MLOGS(IREC)%BNC)) then
        allocate(MLOGS(IREC)%BNC(4*MLOGS_MAXM,MLOGS_MAXZ),STAT=ierr)
        if (IERR==0) then
          write(*,*) 'MLOGS_REC_ALLOC '//trim(MLOGS(IREC)%ID),' ',IREC
          call MLOGS_REC_CLR(IREC)
        else
          write(*,*) 'failed MLOGS_REC_ALLOC '//trim(MLOGS(IREC)%ID),' ',IREC
        endif
      endif
      if (.NOT.allocated (MLOGS(IREC)%TMPBNC)) then
        allocate(MLOGS(IREC)%TMPBNC(3,MLOGS_MAXZ),STAT=ierr)
        if (IERR==0) then
          write(*,*) '   OK'
        else
          write(*,*) '   failed'
        endif
      endif
      if (.NOT.allocated (MLOGS(IREC)%TMPCOORD)) then
      ! record x,y,z,kx,ky,kz
        allocate(MLOGS(IREC)%TMPCOORD(6,MLOGS_MAXZ),STAT=ierr)
      endif

      end subroutine MLOGS_REC_ALLOC

!---------------------------------------------------------
      subroutine MLOGS_REC_FREE(IREC)
! free bounce logger
!---------------------------------------------------------
      integer,intent(in) :: IREC
      if (allocated (MLOGS(IREC)%BNC)) then
        write(*,*) 'MLOGS_REC_FREE '//trim(MLOGS(IREC)%ID),' ',IREC
        deallocate(MLOGS(IREC)%BNC)
      endif
      if (allocated (MLOGS(IREC)%TMPBNC)) then
        deallocate(MLOGS(IREC)%TMPBNC)
      endif
      if (allocated (MLOGS(IREC)%TMPCOORD)) then
        deallocate(MLOGS(IREC)%TMPCOORD)
      endif
      MLOGS(IREC)%LOGBNC=0
      MLOGS(IREC)%NBNC=0
      end subroutine MLOGS_REC_FREE

!---------------------------------------------------------
      subroutine MLOGS_REC_TMP(IREC,IC,R,K,QM)
! update temporary bounce logger
! IC = left (1), right (2), top (3), bottom (4)
! Z = position along Z [mm], local coord.
! QM = m-value for the reflection
!---------------------------------------------------------
      integer,intent(in) :: IREC,IC
      real(kind(1.D0)),intent(in) :: R(3),K(3),QM
      integer :: np
1     format(a,6(1x,G12.5))

      !dbg=(idbg<50)
      !return

      if (.not.MLOGS_isValid(IREC)) then
        if (dbg) write(*,1) 'MLOGS_REC_TMP invalid IREC ',IREC
        return
      !else
      !  if (dbg) write(*,1) 'MLOGS_REC_TMP IREC= ',IREC
      endif


        np=MLOGS(IREC)%NBNC
        !idbg=idbg+1
        !if (dbg) write(*,1) 'MLOGS_REC_TMP '//trim(MLOGS(IREC)%X%FRAME%ID),IREC,idbg
        if (np<MLOGS_MAXZ) then
          np=np+1
          MLOGS(IREC)%NBNC=np
          !if (dbg) write(*,1) '    ',MLOGS(IREC)%NBNC,IC,Z,QM
          MLOGS(IREC)%TMPBNC(1:3,np)=(/IC*1.D0,R(3),QM/)
          MLOGS(IREC)%TMPCOORD(1:6,np)=(/R(1),R(2),R(3),K(1),K(2),K(3)/)
        endif

      end subroutine MLOGS_REC_TMP


!---------------------------------------------------------
      subroutine MLOGS_REC_UPD(IREC,P)
! update bounce logger with the prehiousle stored bouncees history in TMPBNC
! P = event weight
!---------------------------------------------------------
      integer,intent(in) :: IREC
      real(kind(1.D0)),intent(in) :: P
      integer :: iz,im,k,ioff,ic,np
      real(kind(1.D0)) :: qm,z,dz,dq
1     format(a,6(1x,G12.5))

      if (.not.MLOGS_isValid(IREC)) then
        write(*,1) 'MLOGS_REC_UPD invalid IREC ',IREC
        return
      endif
      !return
      if (MLOGS(IREC)%NBNC>0) then
        dz=MLOGS(IREC)%ZRANGE
        dq=MLOGS_MAXVAL
        np=MLOGS(IREC)%NBNC
        idbg=idbg+1
        !dbg=(idbg<50)
        if (dbg) write(*,1) 'MLOGS_REC_UPD ',IREC,np,dz,dq
        do k=1,np
          ic=NINT(MLOGS(IREC)%TMPBNC(1,k))
          z=MLOGS(IREC)%TMPBNC(2,k)
          qm=MLOGS(IREC)%TMPBNC(3,k)
          if (dbg) write(*,1) '    k,ic,z,q ',k,ic,z,qm
          ioff=(ic-1)*MLOGS_MAXM
          iz=nint(MLOGS_MAXZ*z/dz+0.5)
          im=ioff+nint(MLOGS_MAXM*qm/dq+0.5)
          if (dbg) write(*,1) '    ioff, iz,im ',ioff,iz,im
          if ((iz>0).and.(iz<=MLOGS_MAXZ)) then
            if ((im>0).and.(im<=4*MLOGS_MAXM)) then
              MLOGS(IREC)%BNC(im,iz)=MLOGS(IREC)%BNC(im,iz)+P
            endif
          endif
        enddo
      endif
      end subroutine MLOGS_REC_UPD


!---------------------------------------------------------
      subroutine MLOGS_PRINT_TMP(IREC,IC,IU)
! write guide bounce logger map in the given file unit
! IC ... which coordinate: x(1) or y(2)
!---------------------------------------------------------
      integer,intent(in) :: IREC
      integer,intent(in) :: IU,IC
      integer :: k,is
      real(kind(1.D0)) :: RC(4)
1     format(a,6(1x,G12.5))
      write(IU,*) '# MLOGS history, '//trim(MLOGS(IREC)%ID)
      write(IU,*) '# i, side, z, x, kx, qm'
      if (MLOGS(IREC)%NBNC>0) then
        do k=1,MLOGS(IREC)%NBNC
          is=NINT(MLOGS(IREC)%TMPBNC(1,k))
          RC(1)=MLOGS(IREC)%TMPCOORD(3,k)
          RC(2)=MLOGS(IREC)%TMPCOORD(IC,k)
          RC(3)=MLOGS(IREC)%TMPCOORD(3+IC,k)
          RC(4)=MLOGS(IREC)%TMPBNC(3,k)
          write(IU,1) '   ',k,is,RC
        enddo
      endif
      end subroutine MLOGS_PRINT_TMP

!---------------------------------------------------------
      subroutine MLOGS_REC_OUT(IREC,IU)
! write guide bounce logger map in the given file unit
!---------------------------------------------------------
      integer,intent(in) :: IREC
      integer,intent(in) :: IU
      character(LEN_LINE) :: LINE1,LINE2
      real(kind(1.D0)) :: X(2),MV(MLOGS_MAXM)
      integer :: i,j,k,k1,k2,IX(2)
1     format(a,6(1x,G12.5))
      if (.not.MLOGS_isValid(IREC)) then
        write(*,1) 'MLOGS_REC_OUT invalid IREC ',IREC
        return
      endif
      write(IU,*) '# guide ID: '//trim(MLOGS(IREC)%ID)
      X=(/0.D0,MLOGS(IREC)%ZRANGE/)
      call ARRAY2STR(X(1:2),2,LINE2)
      write(IU,*) '# z-range [mm]: '//trim(LINE2)
      X=(/0.D0,MLOGS_MAXVAL/)
      call ARRAY2STR(X(1:2),2,LINE2)
      write(IU,*) '# m-range [Ni]: '//trim(LINE2)
      select case (MLOGS(IREC)%LOGBNC)
      case(mlog_left,mlog_right,mlog_top,mlog_bottom)
        IX=(/MLOGS(IREC)%LOGBNC,MLOGS_MAXM/)
        call IARRAY2STR(IX,2,LINE2)
        write(IU,*) '# wall no., columns: ['//trim(LINE2)//']'
        k1=(MLOGS(IREC)%LOGBNC-1)*MLOGS_MAXM+1
        k2=k1+MLOGS_MAXM-1
        do i=1,MLOGS_MAXZ
          call ARRAY2STR(MLOGS(IREC)%BNC(k1:k2,i),MLOGS_MAXM,LINE2)
          write(IU,*) trim(LINE2)
        enddo
      case(mlog_leftright)
      ! left,right
        write(IU,*) '# sum of left & right walls'
        k1=0*MLOGS_MAXM
        k2=1*MLOGS_MAXM
        do i=1,MLOGS_MAXZ
          do j=1,MLOGS_MAXM
            MV(j)=MLOGS(IREC)%BNC(k1+j,i)+MLOGS(IREC)%BNC(k2+j,i)
          enddo
          call ARRAY2STR(MV,MLOGS_MAXM,LINE1)
          write(IU,*) trim(LINE1)
        enddo
      case(mlog_topbottom)
      ! top, bottom
        write(IU,*) '# sum of top & bottom walls'
        k1=2*MLOGS_MAXM
        k2=3*MLOGS_MAXM
        do i=1,MLOGS_MAXZ
          do j=1,MLOGS_MAXM
            MV(j)=MLOGS(IREC)%BNC(k1+j,i)+MLOGS(IREC)%BNC(k2+j,i)
          enddo
          call ARRAY2STR(MV,MLOGS_MAXM,LINE1)
          write(IU,*) trim(LINE1)
        enddo
      case(mlog_all)
      ! all
        write(IU,*) '# sum of all walls'
        do i=1,MLOGS_MAXZ
          do j=1,MLOGS_MAXM
            MV(j)=0.D0
            do k=0,3
              MV(j)=MV(j)+MLOGS(IREC)%BNC(k*MLOGS_MAXM+j,i)
            enddo
          enddo
          call ARRAY2STR(MV,MLOGS_MAXM,LINE1)
          write(IU,*) trim(LINE1)
        enddo
      end select
      end subroutine MLOGS_REC_OUT


      end module MIRRLOG