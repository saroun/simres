!//////////////////////////////////////////////////////////////////////
!////  $Id: tables.f,v 1.6 2013/07/12 08:16:51 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.6 $
!////     $Date: 2013/07/12 08:16:51 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Tables management
!//////////////////////////////////////////////////////////////////////
      MODULE TABLES
      use FIELDDEF
      use ENUMDEF
      USE FILETOOLS
      USE CLASSES
      USE TABLE_ATOMS
      implicit none

      TYPE TCR_PARAM; SEQUENCE
        character*8  ID               !   ID string[8] of the crystal/reflection, e.g. " Si 111 "
        real(kind(1.D0)) ::  dhkl     !   d-spacing [A]
        real(kind(1.D0)) ::  QML      !   4*PI*(F*dhkl/V0)**2 [ A^-1 cm^-1] (Maier-Leibnitz reflectivity)
        real(kind(1.D0)) ::  sigmab   !   bound-atom scattering cross-section [barn]
        real(kind(1.D0)) ::  sigmaa   !   absorption for 1A neutrons [barn*A^-1]
        real(kind(1.D0)) ::  sigmai   !   incoherent scattering cross-section [barn]
        real(kind(1.D0)) ::  VOL      !   volume [A^3]/atom
        real(kind(1.D0)) ::  A        !   atomic number
        real(kind(1.D0)) ::  thetaD   !   Debye temperature (K)
        real(kind(1.D0)) ::  C2       !   constant from the Freund's paper  [A^-2 meV^-1]
        real(kind(1.D0)) ::  poi      !   Poisson elastic constant
      end type TCR_PARAM

      TYPE(TCR_PARAM) :: CR_TABLE(100)
      integer,private :: N_TABLE=0

      contains

!---------------------------------------------------------------------
      SUBROUTINE READ_ATOMS
! read atoms.dat and fill the table in TABLE_ATOMS
!---------------------------------------------------------------------
      character(MAX_FNAME_LENGTH) :: FRES
      CHARACTER*128 :: LINE
      integer :: ie,ILINE,IU,IRES
      logical :: error
      character*(32) :: S
1     FORMAT(a)
      error=.false.
      ! CALL DLG_FILEOPEN('atoms',trim(CFGPATH)//'/tables','lib',0,0,IRES,FRES)
      CALL DLG_FILEOPEN('atoms',trim(CFGPATH)//'/tables','dat',0,0,IRES,FRES)
      if (IRES.gt.0) then
        call XML_RSXDUMP(SMES,' ',1)
        IU=OPENFILEUNIT(FRES,.true.)
        if (IU.le.0) then
          call MSG_ERROR('READ_ATOMS','Cannot open atoms.dat',0,0)
        ELSE
          ILINE=0
          ie=0
          DO WHILE ((ie.EQ.0).and.(.not.error))
            ILINE=ILINE+1
            LINE=' '
            READ(IU,1,iostat=ie) LINE
            call TABLE_ATOMS_READLN(LINE,error)
          enddo
          close(IU)
          if (error) then
            call INT2STR(ILINE,S)
            call MSG_WARN('Format error in atoms.dat, line'//trim(S),0)
          endif
          call MSG_INFO('Atoms library read from '//trim(FRES),0)
        endif
        call XML_RSXDUMP(SMES,' ',0)
      endif
      end SUBROUTINE READ_ATOMS


!---------------------------------------------------------------------
      SUBROUTINE READCRYST
! read crystal.dat and fill the table CR_TABLE
!---------------------------------------------------------------------
      character(MAX_FNAME_LENGTH) :: FRES
      CHARACTER*128 LINE
      character(LEN_ID) :: CRNAME
      INTEGER*4 ie,ILINE,IU,ires,is,il,NC,i,j
      TYPE(TCR_PARAM) :: R
      integer :: ITYPE,IORD
1     FORMAT(a)
2     FORMAT('Format error in the crystal library, line ',I5 )
4     FORMAT('Cannot open crystal.tab' )
5     FORMAT('WARNING: Incompatible crystal library, ',a,' is not defined in ',a )

      ! CALL DLG_FILEOPEN('crystal','|'//trim(RESPATH)//'|'//trim(CFGPATH),'lib',0,0,IRES,FRES)
      CALL DLG_FILEOPEN('crystal',trim(RESPATH)//'/tables|'//trim(CFGPATH)//'/tables','tab',0,0,IRES,FRES)
      if (IRES.gt.0) then
        IU=OPENFILEUNIT(FRES,.true.)
        if (IU.le.0) then
          write(*,4)
          return
        endif
        call XML_RSXDUMP(SMES,' ',1)
        ILINE=0
        ie=0
        N_TABLE=0
        do I=1,100
          CR_TABLE(I)%ID=""
        enddo
        DO WHILE (ie.EQ.0)
          ILINE=ILINE+1
          LINE=' '
          READ(IU,1,iostat=ie) LINE
        ! skip comments and empty lines
          if ((LINE(1:1).ne.'#').and.(len_trim(line).gt.0)) then
          if (LEN_TRIM(LINE).le.25) ie=1
          if (ie.eq.0) then
            CALL BOUNDS(LINE(1:8),is,il)
            R%ID=LINE(is:is+il-1)
          ! IORD is indexed from 0 !!
            call GetEnumOrd('CRYSTAL_REFNAME',R%ID,ITYPE,IORD)
          !write(*,*) 'READCRYST [',trim(R%ID),'][',LINE(is:is+il-1),']',ITYPE,IORD
            if (IORD.lt.0) then
              write(LINE,5) trim(R%ID),trim(XCONFIG_CLASSES)
              call MSG_WARN(trim(LINE),0)
            ELSE
              if (IORD.lt.100) then
                READ(LINE(9:),*,iostat=ie) R%DHKL,R%QML,R%sigmab,R%sigmaa,R%sigmai,R%VOL,R%A,R%thetaD,R%C2,R%poi
                R%C2=R%C2/1000  ! from eV^-1 to meV^-1
                if (ie.eq.0) then
                  N_TABLE=MAX(N_TABLE,IORD+1)
                  CR_TABLE(IORD+1)=R
                else
                  write(*,2) ILINE
                  exit
                endif
              ENDIF
            ENDIF
          endif
          endif
        endDO
        close(IU)
! problem with type definition => return
        if (ITYPE.le.0) return
! pack the array
        do i=1,N_TABLE
          if (len_trim(CR_TABLE(I)%ID).eq.0) then
            do j=i,N_TABLE-1
              CR_TABLE(J)=CR_TABLE(J+1)
            enddo
            N_TABLE=N_TABLE-1
          endif
        ENDDO
! check consistence with GUI
      !  call COUNTPAR(ENUM_TYPE_LIST(ITYPE) ,':',NC)
        NC=GetEnumDim(ITYPE)
        do IORD=0,NC-1
          call GetEnumString('CRYSTAL_REFNAME',IORD,CRNAME)
          if (.not.FINDCRYSTPARAM(CRNAME,R)) then
            call MSG_WARN(trim(CRNAME)//' not found in the crystal library '//trim(FRES),0)
          endif
        enddo
        call MSG_INFO('Crystal library read from '//trim(FRES),0)
        call XML_RSXDUMP(SMES,' ',0)
      ENDif
      end SUBROUTINE READCRYST


!---------------------------------------------------------------------
      logical function FINDCRYSTPARAM(CRNAME,PARAM)
! find item woth given name and return a record with parameters
! return true if the item was found
!---------------------------------------------------------------------
      character(*),intent(in) :: CRNAME
      TYPE(TCR_PARAM),intent(out) :: PARAM
      character(8) :: CNAME,CITEM
      integer :: i
      logical:: log1,dbg
      CNAME=trim(CRNAME)
      call Mkupcase(CNAME)
      i=0
      log1=.false.
      dbg = .false.
      if (dbg) write(*,*) 'FINDCRYSTPARAM ',trim(CNAME),N_TABLE
      do while ((.not.log1).and.(i.lt.N_TABLE))
        i=i+1
        CITEM=CR_TABLE(i)%ID
        call Mkupcase(CITEM)
        if (dbg) write(*,*) i,' ',trim(CITEM)
        if (CITEM.eq.CNAME) then
          PARAM=CR_TABLE(i)
          log1=.true.
          if (dbg) write(*,*) 'found'
        endif
      enddo
      FINDCRYSTPARAM=log1
      end function FINDCRYSTPARAM

      end MODULE TABLES