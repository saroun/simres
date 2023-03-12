!C//////////////////////////////////////////////////////////////////////
!C////  $Id: strain_table.f,v 1.9 2019/06/22 22:57:16 saroun Exp $
!C////
!C////  R E S T R A X - Simulation of neutron three-axis spectrometers
!C////
!C////     Copyright (C) 1995-2006, All rights reserved
!C////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!C////     Institut Laue Langevin, Grenoble, France
!C////
!C////     Written by:  Jan Saroun
!C////     $Revision: 1.9 $
!C////     $Date: 2019/06/22 22:57:16 $
!C//////////////////////////////////////////////////////////////////////
!C////
!C////  N E S S - Ray-tracing package for RESTRAX
!C////
!C////  Table for macroscopic strain (depth distribution)
!C////
!C//////////////////////////////////////////////////////////////////////
      module STRAIN_TABLE
      use CONSTANTS
      use XMLINFO
      USE SIMATH
      IMPLICIT NONE
      private
      save

      integer, parameter :: STRAIN_DIM=5   ! max. number of tables
      INTEGER,parameter :: s_nrow=256 ! max. number of table rows

      integer :: s_num=0               ! number of actually used tables
      INTEGER :: s_n(STRAIN_DIM)       ! number of rows for each table
      integer :: s_used(6,STRAIN_DIM)  ! 1 if s_eps(i)<>0
      integer :: p_used(3,STRAIN_DIM)  ! 1 if p_depth(i)<>0
      character(512) :: s_files(STRAIN_DIM) ! filenames
      real(kind(1.D0)) :: s_depth(s_nrow,STRAIN_DIM)  ! depth [mm]
      real(kind(1.D0)) :: s_eps(s_nrow,6,STRAIN_DIM)  ! strain tensor elements
! s_eps(6) = e_11,e_22,e_33,e_12,e_23,e_13
      real(kind(1.D0)) :: p_depth(s_nrow,3,STRAIN_DIM)  ! scattering probability along x,y,z minus 1

      public STRAIN_CLEARALL,STRAIN_GET,SET_STABLE,STRAIN_ISDEF
      contains


!-----------------------------------------------------------
      LOGICAL FUNCTION STRAIN_ISDEF(itab)
! returns true if the table is defined
!-----------------------------------------------------------
      integer,intent(in) :: itab
      STRAIN_ISDEF=((itab>0).and.(itab.le.s_num))
      end FUNCTION STRAIN_ISDEF

!---------------------------------------------------------------
      subroutine STRAIN_CLEARALL
! clear all tables
!---------------------------------------------------------------
      integer :: i
      s_num=0
      s_n=0
      do i=1,STRAIN_DIM
        s_files(i)=""
        s_used(1:6,i)=0
        p_used(1:3,i)=0
      enddo
      end subroutine STRAIN_CLEARALL

!---------------------------------------------------------------
      subroutine STRAIN_WRITEXML(IU)
! write list of tables in XML
!---------------------------------------------------------------
      integer, intent(in) :: IU
      integer :: i
      character(32) :: CNUM1
1     FORMAT('<STRAINLIST length="',a,'">')
2     FORMAT('<ITEM key="',a,'" >',a,'</ITEM>')
3     FORMAT('</STRAINLIST>')

      call XML_RSXDUMP(IU,'',1)
      call INT2STR(s_num,CNUM1)
      write(IU,1) trim(CNUM1)
      do i=1,s_num
        call INT2STR(i,CNUM1)
        write(IU,2) trim(CNUM1),trim(s_files(i))
      enddo
      write(IU,3)
      call XML_RSXDUMP(IU,'',0)
      end subroutine STRAIN_WRITEXML


!---------------------------------------------------------------
      INTEGER FUNCTION STRAIN_READ(FNAME,ID)
! read table with depth profiles of strain tensor components
! use FILETOOLS wrappers - getting file data via GUI
! file structure: 1 line header + 7 columns: depth, eps(6)
! returns:
! if OK, index to lookup table
! -1, if table is full
! -2, if file format error occures (less than 3 valid lines)
! 0, if ID outside valid range (1,STRAIN_DIM)
!---------------------------------------------------------------
      use FILETOOLS
      use IO
      IMPLICIT NONE
      INTEGER,intent(in) :: ID
      character*(*) :: FNAME
      character(128) :: MSG,LINE
      character(32) :: CNUM1,CNUM2,CNUM3
      INTEGER :: ierr,i,j,itab,res,IU,ird,nrow,iline
      real(kind(1.D0)) :: a(10), sumeps(6), sump(3)
1     format(a,10(1x,G10.4))
3     FORMAT('Strain table (',a,'/',a') read, ',a,' lines')
4     FORMAT('Error while reading strain table ',a,', line ',a, ', ierr=',a)
10    FORMAT(a)
! JS5/19      IU=OPENFILE_READ(FNAME,11)
      IU = OPENFILEUNIT(FNAME,.true.)
      nrow=0
      iline=0
      res=0
      MSG=''
      if ((ID<1).or.(ID>STRAIN_DIM)) then
        call INT2STR(ID,CNUM1)
        call INT2STR(STRAIN_DIM,CNUM2)
        MSG='Table index outside range: '//trim(CNUM1)//', max='//trim(CNUM2)
        goto 99
      endif
      itab=ID
      ! JS5/19 ird=FILE_READ(IU,LINE) ! header
      Read(IU,10,iostat=IERR,end=88,err=99) LINE
      iline=iline+1
      ird = len_trim(LINE)
      !write(*,*) ird, trim(LINE)
      ierr=0
      s_used(1:6,itab)=0
      p_used(1:3,itab)=0
      do while((ird>3).and.(nrow<s_nrow).and.(ierr==0))
        ! ird=FILE_READARRAY(IU,A,10,1,10)
        ! JS5/19 ird=FILE_READ(IU,LINE)
        Read(IU,10,iostat=IERR,end=88,err=99) LINE
        iline=iline+1
        ! skip comments
        if (LINE(1:1).ne.'#') then
        ! read table row
          call READ_ARRAY8(LINE,10,a,ird)
          ! write(*,1) 'STRAIN_READ: ',a(1:6),ird
          if (ird>3) then
            nrow=nrow+1
            s_depth(nrow,itab)=a(1)
            ! e_11, e_22, e_33, e_23, e_13, e_12
            do i=2,min(ird,7)
              s_eps(nrow,i-1,itab)=a(i)*1.D-3 ! table is in units [0.001]
            enddo
            ! clear e_23, e_13, e_12 if not read
            do i=ird+1,7
              s_eps(nrow,i-1,itab)=0.D0
            enddo
            ! read px,    py,    pz if present
            if (ird>9) then
              do i=1,3
                p_depth(nrow,i,itab)=A(i+7)
              enddo
            ! clear  px,    py,    pz  if not read
            else
              do i=1,3
                p_depth(nrow,i,itab)=0.D0
              enddo
            endIF
          !else if (ird.ne.0) then
            !ierr=1
          endif
        endif
      enddo
      ! indicate which columnts are used (must have at least 1 row <> 0)
88    sumeps=0.D0
      sump=0.D0
      ierr=0
      do i=1,nrow
        do j=1,6
          sumeps(j) = sumeps(j) + abs(s_eps(i,j,itab))
        enddo
        do j=1,3
          sump(j) = sump(j) + abs(p_depth(i,j,itab))
        enddo
      enddo
      do j=1,6
        s_used(j,itab)=LOG2INT(sumeps(j).gt.0.D0)
      enddo
      do j=1,3
        p_used(j,itab)=LOG2INT(sump(j).gt.0.D0)
      enddo

      !write(*,1) ' last row ',nrow,ird,A
      s_n(itab)=nrow
      !write(*,1) ' STRAIN_READ ',itab,s_r(itab),nrow
      ! JS5/19 call CLOSEFILE_READ(IU)
      close(IU)
      if ((ierr==0).and.(nrow>2)) then
    ! write message
        s_num=max(itab,s_num)
        call INT2STR(itab,CNUM1)
        call INT2STR(s_num,CNUM2)
        call INT2STR(nrow,CNUM3)
        write(*,3) trim(CNUM1),trim(CNUM2),trim(CNUM3)
        write(*,1) 'file '//trim(FNAME)
        !write(*,1) '-------------------'
        res=itab
      else
        res=-2
        s_n(itab)=0
        goto 99
      endif
      STRAIN_READ=res
      return
! error:
! JS5/19
!99    call CLOSEFILE_READ(IU)
99    close(IU)
      call INT2STR(iline,CNUM1)
      call INT2STR(ierr,CNUM3)
      write(*,4) trim(FNAME),trim(CNUM1),CNUM3
      if (len_trim(MSG)>0) write(*,1) trim(MSG)
      STRAIN_READ=res

      end FUNCTION STRAIN_READ

!-----------------------------------------------------------
      subroutine STRAIN_GET(db, itab,depth,E, P)
! returns 6 components of the strain tensor for given depth 
!-----------------------------------------------------------
      logical,intent(in) :: db
      integer,intent(in) :: itab
      REAL(kind(1.D0)),intent(in) :: depth
      REAL(kind(1.D0)),intent(out) :: E(6), P(3)
      REAL(kind(1.D0)) :: dd,z,d(3),delta
      integer :: iz,NR,N,i
1     format(a,': ',6(G13.6,1x))
      E(1:6)=0.D0
      P(1:3)=1.D0
      if ((itab<1).or.(itab>s_num)) return
      NR=itab
      N=s_n(NR)
      dd=(s_depth(N,NR)-s_depth(1,NR))/(N-1)
      z=(depth-s_depth(1,NR))/dd
      if (db) write(*,1) 'STRAIN_GET ', N, z, dd
	  if (db) write(*,1) 's_used', s_used(1:6,NR)
      if (db) write(*,1) 'p_used', p_used(1:3,NR)
      if (z<0.D0) then
        do i=1,6
          E(i)=s_eps(1,i,NR)
        enddo
        do i=1,3
		   if (p_used(i,NR)>0) then
              P(i)=p_depth(1,i,NR)
		   endif
        enddo
        return
      else if (z>=1.D0*(N-1)) then
        do i=1,6
          E(i)=s_eps(N,i,NR)
        enddo
        do i=1,3
		   if (p_used(i,NR)>0) then
              P(i)=p_depth(1,i,NR)
		   endif
        enddo
        return
      else
        iz=INT(z)+1
        if (iz>N-2) iz=N-2
        delta=depth-s_depth(iz,NR)
        if (db) write(*,1) '   delta', delta, iz

        do i=1,6
          if (s_used(i,NR)>0) then
            call quadinterp3(s_depth(iz:iz+2,NR),s_eps(iz:iz+2,i,NR),d)
            if (db.and.(i==1)) then
              write(*,1) 'interp', s_depth(iz:iz+2,NR), s_eps(iz:iz+2,i,NR)
              write(*,1) '     d', d
            endif
            E(i)=d(1)*delta**2+d(2)*delta+d(3)
          endif
        enddo
        do i=1,3
          if (p_used(i,NR)>0) then
            call quadinterp3(s_depth(iz:iz+2,NR),p_depth(iz:iz+2,i,NR),d)
            P(i)=d(1)*delta**2+d(2)*delta+d(3)
          else
            P(i)=1.D0
          endif
        enddo
      ENDIF
	  if (db) write(*,1)  'STRAIN_GET P', P
      END subroutine STRAIN_GET

!--------------------------------------------------------------------
      subroutine SET_STABLE(CMD)
! wrapper for STRAIN_READ function, to be called from command interpreter
! interpret 1st argument as ID and the 2nd as filename
!--------------------------------------------------------------------
      character(*),intent(in) :: CMD
      integer :: ID=0.D0
      character*(256) :: FNAME
      integer :: IS,IL,IRES,ierr
      IRES=-3
      IS=1
      CALL FINDPAR(CMD,1,IS,IL)
      !write(*,*) CMD(IS:IS+IL-1),' : ',trim(CMD(IS+IL:))
      if (CMD(IS:IS+IL-1).eq.'CLEAR') then
        call STRAIN_CLEARALL
        return
      else if (CMD(IS:IS+IL-1).eq.'DUMP') then
        ! write(*,*) 'SET_MTABLE ',CMD(IS:IS+IL-1)
        call STRAIN_WRITEXML(smes)
        return
      endif

      read(CMD(IS:IS+IL-1),*,ERR=99,END=10,IOSTAT=ierr) ID
10    CALL FINDPAR(CMD,2,IS,IL)
      FNAME=CMD(IS:IS+IL-1)

      IRES=STRAIN_READ(FNAME,ID)
      !write(*,*) 'SET_MTABLE ires=',ires
      if (IRES>0) then
        s_files(IRES)=trim(FNAME)
        return
      endif

99    select case(IRES)
      case(-3)
        call MSG_ERROR('STRAIN_READ','Wrong format of table key, must be a number',0,1)
      case(-2)
        call MSG_ERROR('STRAIN_READ','Wrong table format: '//trim(FNAME),0,1)
      case(-1)
        call MSG_ERROR('STRAIN_READ','No more space for tables.',0,1)
      case(0)
        call MSG_ERROR('STRAIN_READ','Table ID is out of range',0,1)
      !case default

      end select
      end subroutine SET_STABLE

      end  module STRAIN_TABLE

