!//////////////////////////////////////////////////////////////////////
!////  $Id: enumdef.f,v 1.9 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.9 $
!////     $Date: 2019/08/15 15:02:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Definitions of enumerated types
!////   Includes XML handlers for enum definitions in classes.xml
!////
!//////////////////////////////////////////////////////////////////////
      MODULE ENUMDEF
      use XMLUTILS
      use XMLINFO

      implicit none

      private
      save
      integer, parameter :: LEN_EID=16          ! length of enumerator ID
      integer, parameter :: LEN_ENAME=64        ! length of enumerator value
      integer, parameter :: MAX_EITEMS=256      ! max. number of items for one enum. type
      integer, parameter :: MAX_ENUMS=64        ! max. number of enumerated types

      TYPE TENUMDEF
        INTEGER :: IDX=0            ! type index
        character(LEN_EID) :: ID     ! field ID
        character(LEN_ENAME) :: ITEMS(MAX_EITEMS)
        integer :: DIM
      end TYPE TENUMDEF

! encapsulated pointer to TENUMDEF
      type PENUMDEF
        TYPE(TENUMDEF),pointer :: E
  !      TYPE(TENUMDEF),allocatable :: EA(:) ! field definitions
      end TYPE PENUMDEF

! physical storage of enumerators definitions
      integer :: ENUM_NUM=0
      TYPE(PENUMDEF) :: EDEF(MAX_ENUMS)

! private instances for XML import
      TYPE(PENUMDEF) :: myENUMDEF
      TYPE(TENUMDEF) :: ZERO_ENUM

      public LEN_EID
      public GetEnumOrd,GetEnumDim,GetEnumDef,GetEnumString,GetEnumID,ENUM_NUM
      public isEnumDefined,GET_ENUMINDEX,CLEAR_ENUMS
      public PENUMDEF,TENUMDEF
      public startENUMDEF,dataENUMDEF,endENUMDEF
      public myENUMDEF,GET_EMPTY_ENUM,DISPOSE_TMP_ENUMS,report_enum
! flag that enables XML parsing

      interface GetEnumString
        module procedure GetEnumStringS
        module procedure GetEnumStringI
      end interface GetEnumString

      interface GetEnumOrd
        module procedure GetEnumOrdS
        module procedure GetEnumOrdI
      end interface GetEnumOrd

      contains

!----------------------------------------------
      TYPE(TENUMDEF) function GET_EMPTY_ENUM()
!----------------------------------------------
      GET_EMPTY_ENUM=ZERO_ENUM
      end function GET_EMPTY_ENUM

!--------------------------------------------------------
      subroutine CLEAR_ENUMS
! clear definitions of all types
!--------------------------------------------------------
        integer :: i
        call CLEAR_ENUMDEF(ZERO_ENUM)
        do i=1,ENUM_NUM
          call DISPOSE_ENUMDEF(EDEF(i))
        enddo
        ENUM_NUM=0
        call DISPOSE_ENUMDEF(myENUMDEF)
      end subroutine CLEAR_ENUMS

!--------------------------------------------------------
      subroutine DISPOSE_TMP_ENUMS
! dispose temporary instaces of class definitions
!--------------------------------------------------------
        call DISPOSE_ENUMDEF(myENUMDEF)
      end subroutine DISPOSE_TMP_ENUMS

!-------------------------------------------------------------
      character(LEN_EID) function GetEnumID(IENU)
!-------------------------------------------------------------
      integer, intent(in) :: IENU
      if ((IENU.gt.0).and.(IENU.le.ENUM_NUM)) then
        GetEnumID=EDEF(IENU)%E%ID
      else
        GetEnumID=''
      endif
      end function GetEnumID

!-------------------------------------------------------------
      SUBROUTINE GetEnumDef(IENU,OBJ)
!-------------------------------------------------------------
      integer, intent(in) :: IENU
      type(TENUMDEF),intent(out) :: OBJ
      if ((IENU.gt.0).and.(IENU.le.ENUM_NUM)) then
        OBJ=EDEF(IENU)%E
      else
        OBJ=ZERO_ENUM
      endif
      end SUBROUTINE GetEnumDef

C-------------------------------------------------------------
      integer function GetEnumDim(ITYPE)
! get dimesion of enumerator with given type index
! return 0 if enum. not defined
C-------------------------------------------------------------
      integer, intent(in) :: ITYPE
      integer :: res
      res=0
      if ((ITYPE.gt.0).and.(ITYPE.lt.ENUM_NUM)) then
        res=EDEF(ITYPE)%E%DIM
      endif
      GetEnumDim=res
      end function GetEnumDim

C-------------------------------------------------------------
      SUBROUTINE GetEnumOrdS(ETYPE,STRVAL,ITYPE,IORD)
! convert enumerator type and value given by strings to integers
! note: values are indexed from 0 !
! IORD=-1 for valid enum tyoe, but undefined value
! IORD=-2 for invalid enum type
C-------------------------------------------------------------
      character(*) ETYPE,STRVAL
      integer, intent(out) :: ITYPE,IORD
      ITYPE=GET_ENUMINDEX(ETYPE)
      call GetEnumOrdI(ITYPE,STRVAL,IORD)
      end SUBROUTINE GetEnumOrdS

C-------------------------------------------------------------
      SUBROUTINE GetEnumOrdI(IE,STRVAL,IORD)
! convert enumerator type and value given by strings to integers
! note: values are indexed from 0 !
! IORD=-1 for valid enum tyoe, but undefined value
! IORD=-2 for invalid enum type
C-------------------------------------------------------------
      character(*) STRVAL
      integer,intent(in) :: IE
      integer, intent(out) :: IORD
      IORD=-1
      if (IE.gt.0) then
        IORD=GET_ENUMITEM(IE,trim(STRVAL))
      else
        IORD=-2
      endif
      end SUBROUTINE GetEnumOrdI


!--------------------------------------------------------
      subroutine GetEnumStringS(ETYPE,IORD,ESTR)
! return string value for given enumerated type and ord number
!--------------------------------------------------------
      character(*),intent(in) :: ETYPE
      integer,intent(in) :: IORD
      character(*),intent(out) :: ESTR
      integer :: IE
      IE=GET_ENUMINDEX(trim(ETYPE))
      call GetEnumStringI(IE,IORD,ESTR)
      END subroutine GetEnumStringS

!--------------------------------------------------------
      subroutine GetEnumStringI(IE,IORD,ESTR)
! return string value for given enumerated type and ord number
! IORD is indexed from 0 !
!--------------------------------------------------------
      integer ,intent(in):: IE
      integer,intent(in) :: IORD
      character(*),intent(out) :: ESTR
      if (IE.gt.0) then
        if ((IORD.ge.0).and.(IORD.lt.EDEF(IE)%E%DIM)) then
          ESTR=trim(EDEF(IE)%E%ITEMS(IORD+1))
        else if (EDEF(IE)%E%DIM>0) then
          ESTR=trim(EDEF(IE)%E%ITEMS(1))
        else
          ESTR='undefined'
        endif
      endif
      END subroutine GetEnumStringI

!////////////////////////////////////////////////////////////
!/////                                                  /////
!/////                 PRIVATE AREA                     /////
!/////                                                  /////
!////////////////////////////////////////////////////////////


!---------------------------------------------------------
      logical function isEnumDefined(ID)
!---------------------------------------------------------
      character(*) :: ID
      logical :: b
      integer :: i
      b=.false.
      do i=1,ENUM_NUM
        IF (ASSOCIATED(EDEF(i)%E)) then
          b=(trim(ID).eq.trim(EDEF(i)%E%ID))
          if (b) exit
        endif
      enddo
      isEnumDefined=b
      end function isEnumDefined

!---------------------------------------------------------
      subroutine ALLOC_ENUMDEF(E,ID)
!---------------------------------------------------------
        TYPE(PENUMDEF) :: E
        character(*) :: ID
        integer :: ierr
        allocate(E%E,STAT=ierr)
        if(ierr.eq.0) then
          call CLEAR_ENUMDEF(E%E)
          E%E%ID=trim(ID)
        endif
      end subroutine ALLOC_ENUMDEF

!---------------------------------------------------------
      subroutine CLEAR_ENUMDEF(E)
!---------------------------------------------------------
        TYPE(TENUMDEF) :: E
        integer :: i
        E%ID=' '
        do i=1,MAX_EITEMS
          E%ITEMS(i)=' '
        enddo
        E%DIM=0
      end subroutine CLEAR_ENUMDEF

!---------------------------------------------------------
      subroutine DISPOSE_ENUMDEF(E)
!---------------------------------------------------------
        TYPE(PENUMDEF) :: E
        if (associated(E%E)) then
          DEallocate(E%E)
          nullify(E%E)
        endif
      end subroutine DISPOSE_ENUMDEF

!---------------------------------------------------------
      integer function GET_ENUMINDEX(ID)
! get index of enuerator with given ID
!---------------------------------------------------------
      character(*) :: ID
      integer :: i,ix
      ix=0
      i=0
      do while ((i.lt.ENUM_NUM).and.(ix.eq.0))
        i=i+1
        if (trim(ID).eq.trim(EDEF(i)%E%ID)) ix=i
      enddo
      GET_ENUMINDEX = ix
      end function GET_ENUMINDEX

!---------------------------------------------------------
      integer function GET_ENUMITEM(IE,TXT)
! get ord number of given enum value for an enumerator with given index IE
! return -1 if not found
! NOTE: enum. ord numbers are indexed from 0 !
!---------------------------------------------------------
      character(*) :: TXT
      integer IE
      integer :: i,ix
      ix=0
      i=0
      !if (IE.eq.18) write(*,*) 'GET_ENUMITEM  ',IE,ENUM_NUM,' [',trim(TXT),']'
      if (IE.le.ENUM_NUM) then
        if (associated(EDEF(IE)%E)) then
          !if (IE.eq.23) write(*,*) '    DIM=',EDEF(IE)%E%DIM
          do while ((i.lt.EDEF(IE)%E%DIM).and.(ix.eq.0))
            i=i+1
            !if (IE.eq.23) write(*,*) i,' ID=[',trim(EDEF(IE)%E%ITEMS(i)),']'
            if (trim(TXT).eq.trim(EDEF(IE)%E%ITEMS(i))) ix=i
            !if (IE.eq.18) then
            !  write(*,*) '     [',trim(TXT),'] [',trim(EDEF(IE)%E%ITEMS(i)),'] ',ix
            !endif
          enddo
        endif
      endif
      GET_ENUMITEM = ix-1
      end function GET_ENUMITEM

!---------------------------------------------------------
      subroutine ADD_ENUMDEF(ENUM,error)
!---------------------------------------------------------
      type(TENUMDEF),pointer :: ENUM
      logical,intent(out) :: error
      integer :: I
      error=.true.
      if ((ENUM_NUM.lt.MAX_ENUMS).and.(.not.isEnumDefined(ENUM%ID))) then
        I=ENUM_NUM+1
        call ALLOC_ENUMDEF(EDEF(I),ENUM%ID)
        if (associated(EDEF(I)%E)) then
          ENUM_NUM=I
          EDEF(I)%E=ENUM
          EDEF(I)%E%IDX=I
          error=.false.
        endif
      endif
      end subroutine ADD_ENUMDEF

!---------------------------------------------------------
      subroutine ADD_ENUMITEM(E,TXT)
!---------------------------------------------------------
      TYPE(TENUMDEF),pointer :: E
      character(*) :: TXT
      if (E%DIM.lt.MAX_EITEMS) then
        E%DIM=E%DIM+1
        E%ITEMS(E%DIM)=trim(TXT)
      endif
      end subroutine ADD_ENUMITEM

!--------------------------------------------------------
      subroutine startENUMDEF(ename,error)
!--------------------------------------------------------
        character(len=*)                  :: ename
        logical                           :: error
        if (isEnumDefined(trim(ename))) then
          call MSG_ERROR('FIELDDEF','duplicite definition of enumerated type: '//trim(ename),0,1)
          error=.true.
          return
        else
          call ALLOC_ENUMDEF(myENUMDEF,trim(ename))
        endif
      end subroutine startENUMDEF

!--------------------------------------------------------
      subroutine dataENUMDEF(xmldata)
!--------------------------------------------------------
        character(len=*) :: xmldata
        if (associated(myENUMDEF%E)) then
          call ADD_ENUMITEM(myENUMDEF%E,adjustl(xmldata))
        endif
      end subroutine dataENUMDEF

!--------------------------------------------------------
      subroutine endENUMDEF(error)
!--------------------------------------------------------
      logical :: error
        if (associated(myENUMDEF%E)) then
          call ADD_ENUMDEF(myENUMDEF%E,error)
          call DISPOSE_ENUMDEF(myENUMDEF)
        endif
      end subroutine endENUMDEF

!--------------------------------------------------------
      subroutine report_enum(IENU,IU)
! report class definition
!--------------------------------------------------------
      integer, intent(in) :: IENU
      integer :: IU
      integer :: NP,j
      integer, parameter :: nt=3
      type(TENUMDEF),POINTER :: ED
      character*64  :: TABS
      TABS=' '
2     format(a,I5,'  ','[',a,']')

      if ((IENU.gt.0).and.(IENU.le.ENUM_NUM)) then
        if (associated(EDEF(IENU)%E)) then
          ED=>EDEF(IENU)%E
          NP=ED%DIM
          write(IU,'(I5,a,a16,a,I5)') ED%IDX,TABS(1:nt),adjustl(ED%ID),' items: ',NP
          do j=1,NP
            write(IU,2) TABS(1:nt*2),j-1,trim(ED%ITEMS(j))
          enddo
        endif
      endif
      end subroutine report_enum


      end module ENUMDEF
