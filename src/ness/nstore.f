!//////////////////////////////////////////////////////////////////////
!////  $Id: nstore.f,v 1.34 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2012, All rights Reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.34 $
!////     $Date: 2019/08/15 15:02:08 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Storage of neutron events, new version
!////
!////////////////////////////////////////////////////////////////////////
      MODULE NSTORE
      use TRACINGDATA
      use XMLINFO
      use FILETOOLS
      implicit none
      save
      private
! message ID's
      integer, parameter :: ERR_UNALLOC_ITEM=101
      integer, parameter :: ERR_UNALLOC_EVENT=102
      integer, parameter :: ERR_EXCEED_MAXDATA=103
      integer, parameter :: ERR_EXCEED_MAXSETS=104
      integer, parameter :: ERR_CANNOT_ALLOCATE=105
      integer, parameter :: ERR_CANNOT_DEALLOCATE=106
      integer, parameter :: ERR_CANNOT_REGISTER=107
      integer, parameter :: WARN_EXCEED_EVENT=201
      integer, parameter :: INFO_ALLOCATED=301
      integer, parameter :: INFO_REGISTERED=302
      integer, parameter :: INFO_UNREGISTERED=303

! EVENTS stores following data from NEUTRON type at given indices
! 1..3   <-- R(3)
! 4..6   <-- K(3)
! 7..9   <-- S(3) .. spin
!    10  <-- TOF
!    11  <-- phase
!    12  <-- weight
!    13  <-- T0 (time zero for ToF)

! REFS stores references
! 1 <-- another NRECORD reference
! 2 <-- another NRECORD's index reference
! 3 <-- User's flag 1
! 4 <-- User's flag 2

      integer,parameter :: MAX_MON=2 ! maximum number of monitors
      integer,parameter :: MAX_SUBSETS=8 ! maximum number of subsets
      integer,parameter :: MAX_DATASETS=128 ! maximum number of datasets (channels)
      integer,parameter :: LREC=13 ! record length for events
      integer,parameter :: LREF=3 ! record length for references
      integer,parameter :: MREC=10000 ! segment size in records for memory allocation

      ! reference between a component and event storage
      ! if subset=0, then the reference is made to the monitors
      TYPE NSTOREREF
        character(LEN_ID) :: ID  ! components ID
        integer :: ISET,IDAT,ISTACK  ! subset and dataset indices
      end type NSTOREREF

      TYPE  NRECORD
        integer :: NEV       ! events stored
        integer :: NALLOC    ! events allocated
        integer :: IREF      ! reference to the REGSTORE registry
        real(kind(1.D0)) :: K0(3) ! reference k-vector
        real(kind(1.D0)) :: MCNORM ! normalizing factor
        integer , allocatable :: REFS(:,:)
        real(kind(1.D0)), allocatable :: EVENTS(:,:)
      END TYPE NRECORD

      ! encapsulated pointer to NRECORD
      type PNRECORD
        logical :: OK  ! true if allocated
        character(LEN_ID) :: ID  ! components ID
        TYPE(NRECORD),pointer :: R
      end TYPE PNRECORD

! list of registered storage objects
      integer :: REGSTORE_N
      integer,parameter :: REGSTORE_MAX=MAX_DATASETS*MAX_SUBSETS+MAX_MON
      TYPE(PNRECORD) :: REGSTORE(REGSTORE_MAX)


      public NSTORE_REGISTER,NSTORE_UNREGISTER,NSTORE_UNREGISTER_ALL,NSTORE_GETK0
      public NREC_ALLOCATE,NSTORE_MAXN,NSTORE_GETEVENT,NSTORE_SETEVENT,NSTORE_SETNORM
      public NSTORE_GETNORM,NSTORE_XML_LIST,NSTORE_CLEAR,NSTORE_GETID
      public nstore_dump

      contains

!--------------------------------------------------------------------------
      integer function FIND_FREE_REGISTER()
! return index of the 1st free register = unallocated element of REGSTORE
!--------------------------------------------------------------------------
      integer :: i
      logical :: LOG1
      i=1
      LOG1=.true.
      do while((i<=REGSTORE_MAX).and.LOG1)
        LOG1=associated(REGSTORE(i)%R)
        if (LOG1) i=i+1
      enddo
      if (i>REGSTORE_MAX) i=0
      FIND_FREE_REGISTER=i
      end function FIND_FREE_REGISTER

!--------------------------------------------------------------------------
      integer function NSTORE_REGISTER(ID)
! register a new storage for given subset, dataset and beamline
! return the new reference index
!--------------------------------------------------------------------------
      character(*),intent(in) :: ID
      integer :: ierr,N
1     format(a,6(1x,G12.5))
      N=FIND_FREE_REGISTER()
      if (N>0) then
        ierr=0
        !write(*,1) 'NSTORE_REGISTER '//trim(ID),N
        if (.not.associated(REGSTORE(N)%R)) then
          allocate(REGSTORE(N)%R,STAT=ierr)
        endif
        if (ierr==0) then
          REGSTORE(N)%ID=trim(ID)
          REGSTORE(N)%OK=.true.
          REGSTORE(N)%R%NEV=0
          REGSTORE(N)%R%NALLOC=0
          REGSTORE(N)%R%IREF=N
          REGSTORE(N)%R%K0=0.D0
          REGSTORE(N)%R%MCNORM=0.D0
          call info_msg(INFO_REGISTERED,(/N/))
        else
          N=0
          call err_msg(ERR_CANNOT_REGISTER,(/0/))
        endif
      endif
      REGSTORE_N=MAX(N,REGSTORE_N)
      NSTORE_REGISTER=N
      end function NSTORE_REGISTER

!--------------------------------------------------------------------------
      subroutine NSTORE_UNREGISTER(IREF)
! unregister and destroy registry
!--------------------------------------------------------------------------
      integer,intent(in) :: IREF
      integer :: i
      logical :: LOG1
      !write(*,*) 'NSTORE_UNREGISTER ',IREF
      if ((IREF>0).and.(IREF<=REGSTORE_MAX)) then
        if (REGSTORE(IREF)%OK) then
          call NREC_DESTROY(REGSTORE(IREF)%R)
          deallocate(REGSTORE(IREF)%R)
          NULLIFY(REGSTORE(IREF)%R)
          REGSTORE(IREF)%OK=.false.
          REGSTORE(IREF)%ID=' '
          ! find the last associated registry
          if (REGSTORE_N==IREF) then
            i=REGSTORE_N-1
            LOG1=.false.
            !write(*,*) 'NSTORE_UNREGISTER loop ',i
            do while((i>0).and.(.not.LOG1))
              LOG1=(REGSTORE(i)%OK)
              if (.not.LOG1) i=i-1
            enddo
            !write(*,*) 'NSTORE_UNREGISTER loop end ',i
            REGSTORE_N=i
          endif
          call info_msg(INFO_UNREGISTERED,(/IREF/))
        endif
      endif
      end subroutine NSTORE_UNREGISTER


!--------------------------------------------------------------------------
      subroutine NSTORE_UNREGISTER_ALL
! unregister all registers
!--------------------------------------------------------------------------
      integer :: i
      i=REGSTORE_N
      do while(i>0)
        call NSTORE_UNREGISTER(i)
        i=i-1
      enddo
      end subroutine NSTORE_UNREGISTER_ALL

!--------------------------------------------------------------------------
      subroutine NREC_ALLOCATE(IREF,KREF,NR,PRESERVE)
! Allocates NR records in IREFs storage.
! Allocates more space in memory and preserves data when required.
!--------------------------------------------------------------------------
      integer, intent(in) :: IREF,NR
      real(kind(1.D0)), intent(in) :: KREF(3)
      logical, intent(in) :: PRESERVE
      type(NRECORD),pointer :: RC
      integer :: i,ierr,NAL,NREC,NOLD
      logical :: REPLACE
      real(kind(1.D0)), allocatable :: AUX(:,:)
      integer, allocatable :: IAUX(:,:)
1     format(a,6(1x,G12.5))
      ! not registered

      if (.NOT.((IREF>0).and.REGSTORE(IREF)%OK)) return
      !write(*,1) 'NREC_ALLOCATE ',IREF,NR,REGSTORE(IREF)%R%NALLOC
      RC=>REGSTORE(IREF)%R
      NAL=NR/MREC+1
      NREC=NAL*MREC

      if (.not.PRESERVE) RC%NEV=0
      RC%K0=KREF
      ! already allocated
      !if (NR.le.RC%NALLOC) return
      ierr=0
      ! don't re-allocate memory if not necessary
      if (RC%NALLOC==NREC) return

      REPLACE=(allocated(RC%EVENTS).and.(RC%NEV>0))
      if (REPLACE) then
        NOLD=RC%NEV
        allocate(AUX(1:LREC,1:NOLD),STAT=ierr)
        allocate(IAUX(1:LREF,1:NOLD),STAT=ierr)
        if (IERR==0) then
          do i=1,NOLD
            AUX(1:LREC,i)= RC%EVENTS(1:LREC,i)
            IAUX(1:LREF,i)= RC%REFS(1:LREF,i)
          enddo
        endif
      endif
      if (IERR.ne.0) return
      if (allocated(RC%EVENTS)) deallocate(RC%EVENTS,STAT=ierr)
      if (allocated(RC%REFS)) deallocate(RC%REFS,STAT=ierr)
      !write(*,1) 'NREC_ALLOCATE ierr= ',ierr
      ! number of segments to be allocated
      if (IERR==0) then
        RC%NALLOC=0
        RC%NEV=0
        allocate(RC%EVENTS(1:LREC,1:NREC),STAT=ierr)
        !write(*,1) 'NREC_ALLOCATE EVENTS ierr= ',ierr
        if (ierr==0) allocate(RC%REFS(1:LREF,1:NAL*MREC),STAT=ierr)
        !write(*,1) 'NREC_ALLOCATE REFS ierr= ',ierr
        if (ierr==0) then
          RC%NALLOC=NREC
          if (REPLACE) then
            do i=1,NOLD
              RC%EVENTS(1:LREC,i)= AUX(1:LREC,i)
              RC%REFS(1:LREF,i)= IAUX(1:LREF,i)
            enddo
            RC%NEV=NOLD
            deallocate(AUX,STAT=ierr)
            deallocate(IAUX,STAT=ierr)
          endif
          call info_msg(INFO_ALLOCATED,(/RC%NALLOC,IREF/))
          !write(*,1) 'NREC_ALLOCATE OK ',IREF,REGSTORE(IREF)%R%NALLOC
        else
          call err_msg(ERR_CANNOT_ALLOCATE,(/NREC,IREF/))
        endif
      else
        RC%NEV=0
        RC%NALLOC=0
        call err_msg(ERR_CANNOT_DEALLOCATE,(/IREF/))
      endif
      !write(*,1) 'NREC_ALLOCATE OK '//trim(REGSTORE(IREF)%ID),IREF,REGSTORE(IREF)%R%NALLOC
      end subroutine NREC_ALLOCATE

!--------------------------------------------------------------------------
      subroutine NREC_DESTROY(RC)
! Allocates NR records. Allocates more space in memory when necessary.
! Preserves axisting data.
!--------------------------------------------------------------------------
      type(NRECORD) :: RC
      if (allocated(RC%EVENTS)) deallocate(RC%EVENTS)
      if (allocated(RC%REFS)) deallocate(RC%REFS)
      RC%NALLOC=0
      RC%NEV=0
      end subroutine NREC_DESTROY

!--------------------------------------------------------------------------
      INTEGER FUNCTION NSTORE_MAXN(IREF)
! return number of valid records for given store
!--------------------------------------------------------------------------
        integer , intent(in):: IREF
        integer :: N
        N=0
        !write(*,*) 'NSTORE_MAXN ',IREF,REGSTORE(IREF)%OK,REGSTORE(IREF)%R%NALLOC
        if ((IREF>0).and.REGSTORE(IREF)%OK) then
          if (REGSTORE(IREF)%R%NALLOC>0) then
            N=REGSTORE(IREF)%R%NEV
            !write(*,*) 'NSTORE_MAXN ',N
          endif
        endif
        NSTORE_MAXN=N
      end FUNCTION NSTORE_MAXN

!-----------------------------------------------------------------------
      subroutine NSTORE_GETEVENT(IREF,NEV,NEU)
! store event in given record
!-----------------------------------------------------------------------
      integer, intent(in) :: NEV,IREF
      type(NEUTRON) :: NEU
      REAL(KIND(1.D0)) :: VK(3)
      TYPE(NRECORD),pointer :: RC
      if (CHECK_VALID(IREF,NEV)) then
        RC=>REGSTORE(IREF)%R
        NEU%R(1:3)=RC%EVENTS(1:3,NEV)
        VK(1:3)=RC%EVENTS(4:6,NEV)
        NEU%K=VK+RC%K0
        NEU%S(1:3)=RC%EVENTS(7:9,NEV)
        NEU%T=RC%EVENTS(10,NEV)
        NEU%PHI=RC%EVENTS(11,NEV)
        NEU%P=RC%EVENTS(12,NEV)
        NEU%T0=RC%EVENTS(13,NEV)
        NEU%K0=SQRT(NEU%K(1)**2+NEU%K(2)**2+NEU%K(3)**2)
        NEU%REF(1:2)=RC%REFS(1:2,NEV)
        NEU%LABEL=RC%REFS(3,NEV)
        NEU%CNT=NEV
      endif
      end subroutine NSTORE_GETEVENT


!-----------------------------------------------------------------------
      subroutine NSTORE_SETEVENT(IREF,NEU)
!Write incident neutron to NREC_KI
!-----------------------------------------------------------------------
      integer, intent(in) :: IREF
      type(NEUTRON), intent(in) :: NEU
      REAL(KIND(1.D0)) :: VK(3)
      TYPE(NRECORD),pointer :: RC
      logical :: dbg=.false.
      integer :: NEV
1     format(a,': ',6(G10.3,1x),I5)
      NEV=NEU%CNT
      !dbg=(NEV<100)
      if (NSTORE_CHECK_ACCESS(IREF,NEV)) then
          if (dbg) write(*,1) 'NSTORE_SETEVENT ',IREF,NEV,REGSTORE(IREF)%R%NALLOC
          RC => REGSTORE(IREF)%R
          VK=NEU%K-RC%K0
          RC%EVENTS(1:3,NEV)=NEU%R(1:3)
          RC%EVENTS(4:6,NEV)=VK(1:3)
          RC%EVENTS(7:9,NEV)=NEU%S(1:3)
          RC%EVENTS(10,NEV)=NEU%T
          RC%EVENTS(11,NEV)=NEU%PHI
          RC%EVENTS(12,NEV)=NEU%P
          RC%EVENTS(13,NEV)=NEU%T0
          RC%NEV=MAX(NEV,RC%NEV)
          RC%REFS(1:2,NEV)=NEU%REF(1:2)
          RC%REFS(3,NEV)=NEU%LABEL
          if (dbg) write(*,1) 'NSTORE_SETEVENT OK ',RC%NEV
      endif
      end subroutine NSTORE_SETEVENT

!-----------------------------------------------------------------------
      subroutine NSTORE_SETNORM(IREF,NORM)
! set norm for calculation of absolute intensities
!-----------------------------------------------------------------------
      integer, intent(in) :: IREF
      REAL(KIND(1.D0)), intent(in) :: NORM
      if ((IREF>0).and.REGSTORE(IREF)%OK) then
        REGSTORE(IREF)%R%MCNORM=NORM
        !if (NORM>0.D0) write(*,*) 'NSTORE_SETNORM ' ,IREF,NORM
      endif
      end subroutine NSTORE_SETNORM


!-----------------------------------------------------------------------
      subroutine NSTORE_CLEAR(IREF)
! clear counter of events for given registry
!-----------------------------------------------------------------------
      integer, intent(in) :: IREF
      if ((IREF.GT.0).and.NSTORE_CHECK_ACCESS(IREF,1)) then
        REGSTORE(IREF)%R%NEV=0
      endif
      end subroutine NSTORE_CLEAR


!--------------------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION NSTORE_GETK0(IREF)
! return modulus of the reference k-vector
!--------------------------------------------------------------------------
        integer , intent(in):: IREF
        REAL(KIND(1.D0)) :: K0,K3(3)
        K0=0.D0
        if ((IREF>0).and.REGSTORE(IREF)%OK) then
          K3=REGSTORE(IREF)%R%K0
          K0=SQRT(K3(1)**2+K3(2)**2+K3(3)**2)
        endif
        NSTORE_GETK0=K0
      end FUNCTION NSTORE_GETK0


!--------------------------------------------------------------------------
      SUBROUTINE NSTORE_GETID(IREF,ID)
! return ID of the referenced storage
!--------------------------------------------------------------------------
        integer , intent(in):: IREF
        character(*) :: ID
        integer :: LL,LS,L
        LL=LEN(ID)
        ID=''
        if ((IREF>0).and.REGSTORE(IREF)%OK) then
          LS=LEN(REGSTORE(IREF)%ID)
          L=min(LS,LL)
          ID=REGSTORE(IREF)%ID(1:L)
        endif
      end SUBROUTINE NSTORE_GETID

!--------------------------------------------------------------------------
      REAL(KIND(1.D0)) FUNCTION NSTORE_GETNORM(IREF)
! get norm for calculation of absolute intensities
!--------------------------------------------------------------------------
        integer , intent(in):: IREF
        REAL(KIND(1.D0)) :: Z
        Z=0.D0
        if ((IREF>0).and.REGSTORE(IREF)%OK) then
          Z=REGSTORE(IREF)%R%MCNORM
        endif
        NSTORE_GETNORM=Z
      end FUNCTION NSTORE_GETNORM


!--------------------------------------------------------------------------
      LOGICAL FUNCTION CHECK_VALID(IREF,NEV)
! Check validity for given event = storage must exist and event written
!--------------------------------------------------------------------------
        integer , intent(in):: IREF,NEV
        logical :: LOG1
        LOG1=(IREF>0)
        if (LOG1) then
        ! valid storage
            LOG1=((IREF>0).and.REGSTORE(IREF)%OK)
        ! valid event
            LOG1=(LOG1.and.(REGSTORE(IREF)%R%NEV>=NEV))
        endif
        CHECK_VALID=LOG1
      end FUNCTION CHECK_VALID

!--------------------------------------------------------------------------
      LOGICAL FUNCTION NSTORE_CHECK_ACCESS(IREF,NEV)
! Check access for given event = storage must exist and memory allocated
!--------------------------------------------------------------------------
        integer , intent(in):: IREF,NEV
        logical :: LOG1
1       format(a,': ',6(G12.5,1x))
        !write(*,1) 'NSTORE_CHECK_ACCESS ',IREF,NEV
    ! valid storage
        LOG1=(IREF>0)
        if (LOG1) then
            LOG1=(LOG1.and.(NEV>0))
            LOG1=(LOG1.and.REGSTORE(IREF)%OK)
    ! valid event
            LOG1=(LOG1.and.(REGSTORE(IREF)%R%NALLOC>=NEV))
        endif
        NSTORE_CHECK_ACCESS=LOG1
      end FUNCTION NSTORE_CHECK_ACCESS


!---------------------------------------------------------
      subroutine NSTORE_XML_LIST(IU,SEL)
! fill ARG with SELECT data = list of available storages
!---------------------------------------------------------
      integer,intent(inout) :: IU,SEL
      TYPE(TFIELD) :: ARG
      character(16) :: CNUM1,CNUM2
      integer :: i,isel,n
1     format('<NSTORES length="',a,'" selected="',a,'">')
2     format('</NSTORES>')
        isel=-1
        n=0
        do i=1,REGSTORE_N
          !write(*,*) 'NSTORE_XML_LIST ',REGSTORE(i)%OK,' ',trim(REGSTORE(i)%ID)
          if (CHECK_VALID(i,0)) then
            if (i==SEL) isel=i-1
            n=n+1
          endif
        enddo
        if (isel==-1) isel=n
      ! SEL is indexed from 1, isel from 0 !!
        SEL=isel+1
        call INT2STR(n,CNUM1)
        call INT2STR(isel,CNUM2)
        call XML_RSXDUMP(IU,' ',1)
        write(IU,1) trim(CNUM1),trim(CNUM2)
        do i=1,REGSTORE_N
          if (CHECK_VALID(i,0)) then
            call XML_TAG(IU,'ITEM',trim(REGSTORE(i)%ID),2)
          endif
        enddo
        write(IU,2)
        call XML_RSXDUMP(IU,' ',0)
      end subroutine NSTORE_XML_LIST

!--------------------------------------------------------------------------
      subroutine err_msg(imsg,ipar)
! Fatal error message, STOP
!--------------------------------------------------------------------------
        integer, intent(in) :: imsg
        integer, intent(in) :: ipar(:)
        character(256) :: MSGSTR
1       format('Record is not allocated: REC=',I8,' SUBSET=',I4,' DATASET=',I4,' STACK=',I2)
2       format('Dataset is not allocated: SUBSET=',I4,' DATASET=',I4,' STACK=',I2)
3       format('Exceeded max. number of datasets: ',I4,' MAX=',I4)
4       format('Can''t allocate additional memory for ',I8,' records in stack ',I2)
5       format('Can''t deallocate memory')
6       format('Can''t register new event storage - memory problem?')
        select case (imsg)
        case(ERR_UNALLOC_EVENT)
          WRITE(MSGSTR,1) ipar(1),ipar(2),ipar(3),ipar(4)
        case(ERR_UNALLOC_ITEM)
          WRITE(MSGSTR,2) ipar(1),ipar(2),ipar(3)
        case(ERR_EXCEED_MAXDATA)
          WRITE(MSGSTR,3) ipar(1),ipar(2)
        case(ERR_CANNOT_ALLOCATE)
          WRITE(MSGSTR,4) ipar(1),ipar(2)
        case(ERR_CANNOT_DEALLOCATE)
          WRITE(MSGSTR,5)
        case(ERR_CANNOT_REGISTER)
          WRITE(MSGSTR,6)
        end select
        call MSG_ERROR('NSTACK',MSGSTR,1,1)
        read(*,*)
        stop
      end subroutine err_msg

!--------------------------------------------------------------------------
      subroutine warn_msg(imsg,ipar)
! warning message
!--------------------------------------------------------------------------
        integer, intent(in) :: imsg
        integer, intent(in) :: ipar(:)
        character(256) :: MSGSTR
1       format('Maximum allowed number of events: ',I8,' / required: ',I8)
        select case (imsg)
        case(WARN_EXCEED_EVENT)
          WRITE(MSGSTR,1) ipar(2),ipar(1)
        end select
        call MSG_WARN(MSGSTR,1)
      end subroutine warn_msg

!--------------------------------------------------------------------------
      subroutine info_msg(imsg,ipar)
! info message
!--------------------------------------------------------------------------
        integer, intent(in) :: imsg
        integer, intent(in) :: ipar(:)
        character(256) :: MSGSTR
1       format('Allocated memory for ',I8,' records in storage ',I8)
2       format('Registered storage ',I8,' for ',a)
3       format('Unregistered storage ',I8)
        select case (imsg)
        case(INFO_ALLOCATED)
          WRITE(MSGSTR,1) ipar(1),ipar(2)
        case(INFO_REGISTERED)
          WRITE(MSGSTR,2) ipar(1),trim(REGSTORE(ipar(1))%ID)
        case(INFO_UNREGISTERED)
          WRITE(MSGSTR,3) ipar(1)
        end select
        write(*,*) trim(MSGSTR) ! info messages only to console
      end subroutine info_msg


!--------------------------------------------------------------------------
      subroutine dump_register(iu,ireg)
! dump given register to the given output unit
!--------------------------------------------------------------------------
      integer, intent(in) :: iu,ireg
      character(128) :: CNUM
      TYPE(NRECORD),pointer :: RC
      integer :: i
1     format('# ',a)
2     format(13(G13.5,1x),3(G13.5,1x))
      RC => REGSTORE(ireg)%R
      write(iu,1) 'registry='//trim(REGSTORE(ireg)%ID)
      call INT2STR(RC%NEV,CNUM)
      write(iu,1) 'nrec='//trim(CNUM)
      call ARRAY2STR(RC%K0,3,CNUM)
      write(iu,1) 'k_ref='//trim(CNUM)
      call FLOAT2STR(RC%MCNORM,CNUM)
      write(iu,1) 'mc_norm='//trim(CNUM)
      write(iu,1) 'columns: x,y,z,kx,ky,kz,sx, sy, sz, t,phi,p,t0,ref,ref_cnt,label'
      write(iu,1) 'DATA:'
      do i=1,RC%NEV
        write(iu,2) RC%EVENTS(1:13,i),RC%REFS(1:3,i)
      enddo
      end subroutine dump_register


!--------------------------------------------------------------------------
      subroutine nstore_dump(sarg)
! dump nstore content in a text file
!--------------------------------------------------------------------------
      character(*) :: sarg
      integer :: ireg
      integer :: imin,imax,i,iu,ires
      character(MAX_FNAME_LENGTH) :: FRES,fname
      character(16) :: CNUM

      call READ_I4('REG','REG='//trim(SARG),ireg,ires)
      if (ires.ne.0) then
        ireg=0
      endif
      if (ireg<=0) then
        imin=1
        imax=REGSTORE_N
      else
        imin=ireg
        imax=ireg
      endif
      checkOverwrite=.false.
      do i=imin,imax
        if (NSTORE_CHECK_ACCESS(i,1)) then
          call INT2STR(i,CNUM)
          fname='ndump_'//trim(CNUM)//'_'//trim(REGSTORE(i)%ID)
          call DLG_FILEOPEN(fname,OUTPATH,'dat',0,1,IRES,FRES)
          if (IRES<=0) return
          iu=OPENFILEUNIT(trim(FRES),.false.)
          if (iu>0) then
            call dump_register(iu,i)
            close(iu)
            write(*,*) 'Saved dump file: '//trim(fres)
          endif
        endif
      enddo


      end subroutine nstore_dump


      end MODULE NSTORE