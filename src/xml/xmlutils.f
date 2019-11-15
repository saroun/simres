!//////////////////////////////////////////////////////////////////////
!////  $Id: xmlutils.f,v 1.4 2009/02/25 16:30:53 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.4 $
!////     $Date: 2009/02/25 16:30:53 $
!//////////////////////////////////////////////////////////////////////
!////
!////  XML utilities without dependences on other modules
!////
!//////////////////////////////////////////////////////////////////////
      MODULE XMLUTILS
      implicit none

      interface readClassAttrib
        module procedure readClassAttrib_1
        module procedure readClassAttrib_3
      end interface

      private readClassAttrib_1
      private readClassAttrib_3

      contains


!-----------------------------------------------------------
      subroutine reportAttrib(attribs,IU)
! Print all attributes with values on UNIT=IU
!-----------------------------------------------------------
      character(*), dimension(:,:),intent(in)  :: attribs
      integer,intent(in) :: IU
      integer :: i
      do i=1,size(attribs,2)
        write(IU,"(a10,'= ',a)") trim(attribs(1,i)),trim(attribs(2,i))
      enddo
      end subroutine reportAttrib



!-----------------------------------------------------------
      subroutine readEnumAttrib(attribs,ENUM_TYPE,ENUM_INDEX)
! Search for enumerated type attributes
! Return enumerator type and index
!-----------------------------------------------------------
      character(*), dimension(:,:),intent(in)  :: attribs
      character(*),intent(out) :: ENUM_TYPE
      integer,intent(out) :: ENUM_INDEX
      integer :: i,j
      logical :: isInteger
      ENUM_TYPE=' '
  ! NOTE: index=-1 is allowed and interpreted as "undefined"
      ENUM_INDEX=-2
      do i=1,size(attribs,2)
        select case (trim(attribs(1,i)))
          case('enum')
            ENUM_TYPE=trim(attribs(2,i))
          case('value')
            if (isInteger(trim(attribs(2,i)),j)) ENUM_INDEX=j
        end select
      enddo
  !    if (ENUM_TYPE.eq.' ') ENUM_INDEX=-2
      end subroutine readEnumAttrib


!--------------------------------------------------------
      subroutine readClassAttrib_3(attribs,CLSNAME,IDNAME,CNAME,error)
! search for class attributes
! return class name, ID name and descriptive name of the class
!--------------------------------------------------------
        character(*), dimension(:,:),intent(in)  :: attribs
        character(*),intent(out) :: CLSNAME,IDNAME,CNAME
        logical,intent(out) :: error
        integer :: at,i
        at=0
        do i=1,size(attribs,2)
          select case(trim(attribs(1,i)))
          case('class')
            CLSNAME=trim(attribs(2,i))
            at=at+1
          case('id')
            IDNAME=trim(attribs(2,i))
            at=at+1
          case('name')
            CNAME=trim(attribs(2,i))
            at=at+1
          end select
        enddo
      ! check attributes
        if (at.lt.3) error=.true.
      end subroutine readClassAttrib_3

!--------------------------------------------------------
      subroutine readAttrib(attribs,ANAME,AVALUE,error)
! search for attribute of given name and get its value
!--------------------------------------------------------
        character(*), dimension(:,:),intent(in)  :: attribs
        character(*),intent(in) :: ANAME
        character(*),intent(out) :: AVALUE
        logical,intent(out) :: error
        integer :: at,i
        integer :: IDCOMPARE
        at=0
        AVALUE=' '
        do i=1,size(attribs,2)
          if (IDCOMPARE(trim(attribs(1,i)),trim(ANAME)).eq.0) then
            at=1
            AVALUE=trim(attribs(2,i))
            exit
          else
          endif
        enddo
        error=(at.ne.1)
      end subroutine readAttrib

!--------------------------------------------------------
      logical function readAttribBoolean(attribs,ANAME,DEF)
! search for attribute of given name
! interpret result as boolean variable ('yes/no','true/false','on/off')
! DEF is used as default answer
!--------------------------------------------------------
        character(*), dimension(:,:),intent(in)  :: attribs
        character(*),intent(in) :: ANAME
        logical,intent(in) :: DEF
        character(32) :: AVALUE
        logical :: error
        AVALUE=' '
        call readAttrib(attribs,trim(ANAME),AVALUE,error)
  !      write(*,*) 'readAttribBoolean ',trim(ANAME),' ',trim(AVALUE),error
        if (error) then
          readAttribBoolean=DEF
        else
          call MKUPCASE(AVALUE)
          select case(trim(AVALUE))
          case('YES','Y','TRUE','T','ON','1')
            readAttribBoolean=.true.
          case('NO','N','FALSE','F','OFF','0')
            readAttribBoolean=.false.
          case default
            readAttribBoolean=DEF
          end select
        endif
      end function readAttribBoolean

!--------------------------------------------------------
      integer function readAttribInt(attribs,ANAME,DEF)
! search for integer attribute of given name
! DEF is used as default answer in case of error
!--------------------------------------------------------
        character(*), dimension(:,:),intent(in)  :: attribs
        character(*),intent(in) :: ANAME
        integer,intent(in) :: DEF
        character(32) :: AVALUE
        logical :: error
        integer :: res,I
        logical :: IsInteger
        AVALUE=' '
        res=DEF
        call readAttrib(attribs,trim(ANAME),AVALUE,error)
        if (.not.error) then
          if (isInteger(trim(AVALUE),I)) res=I
        endif
        readAttribInt=res
      end function readAttribInt

!--------------------------------------------------------
      subroutine readClassAttrib_1(attribs,CLSNAME,error)
! search for class attributes
! return class name only
!--------------------------------------------------------
        character(*), dimension(:,:),intent(in)  :: attribs
        character(*),intent(out) :: CLSNAME
        logical,intent(out) :: error
        integer :: at,i
        at=0
        do i=1,size(attribs,2)
          select case(trim(attribs(1,i)))
          case('class')
            CLSNAME=trim(attribs(2,i))
            at=at+1
          end select
        enddo
      ! check attributes
        if (at.lt.1) error=.true.
      end subroutine readClassAttrib_1

      end module XMLUTILS

