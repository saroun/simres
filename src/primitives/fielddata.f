!//////////////////////////////////////////////////////////////////////
!////  $Id: fielddata.f,v 1.20 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.20 $
!////     $Date: 2019/08/15 17:24:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////   Defines abstract TFILED data type, which can encapsulate
!////   - array of integer, float, string or vrange
!////   - enumerated type
!////   - string list type
!////   Implements basic procedures for field data handling
!////
!//////////////////////////////////////////////////////////////////////
      MODULE FIELDDATA
      use XMLINFO
      use FIELDDEF
      use ENUMDEF
      use CLASSDEF

      implicit none

      private
      save

! string field dimensions: keep private
      integer,parameter :: FIELD_BUFFER_FLOAT=128
      integer,parameter :: FIELD_BUFFER_STR=LEN_LINE

! Intrinsic type VRANGE defines a range of a component parameter
      TYPE VRANGE;SEQUENCE
        CHARACTER(2*LEN_ID+1) :: ID ! parameter ID in format COMPONENT.PAR
        REAL(KIND(1.D0)) :: VMI     ! minimum value
        REAL(KIND(1.D0)) :: VMA     ! maximum value
        REAL(KIND(1.D0)) :: VINC    ! increment
      end type VRANGE

! intrinsic field type, can represent an array of integers, floats, strings or value ranges
! fields are dynamically allocated
      TYPE TFIELD
        type(TFIELDDEF),pointer :: DEF=>null()
        CHARACTER(LEN_ID) :: ID=''
        integer :: LENGTH ! allocated array length
        integer :: MAXROW ! allocated number of rows for tables
        integer :: NP     ! used array length (rows for tables)
        integer :: EORD ! ord value of enumerated type
        real(kind(1.D0)), allocatable :: RFIELD(:) ! float arrays
        INTEGER, allocatable :: IFIELD(:)
        CHARACTER(FIELD_BUFFER_STR), allocatable :: SFIELD(:)
        TYPE(VRANGE), allocatable :: VFIELD(:)
      end type TFIELD


      public TFIELD,VRANGE
      public FIELD_COPY,ALLOC_FIELD,ALLOC_TABLE,DISPOSE_FIELD,CLEAR_FIELD
      public FIELD_BUFFER_FLOAT,FIELD_BUFFER_STR
      public STR2FIELD,ARRAY2FIELD,FIELD2STR,FIELD2ARRAY,ITEM2FIELD,FIELD_IS_TYPE
      public STR2ROW,ARRAY2ROW,ROW2STR,ROW2ARRAY

      contains


C---------------------------------------------------------
      SUBROUTINE FIELD2ARRAY(SRC,IDX,ARG,MA,NARG)
! Convert TFIELD object to real array.
! Only for types with numeric value (FLOAT,INT,ENUM)
! return:
! ARG  ... output array
! NARG ... number of processed array elements
C---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IDX,MA
      real(KIND(1.D0)),intent(out) :: ARG(MA)
      integer,intent(out) :: NARG
      integer :: imin,imax,L
      logical :: dbg = .false.
      !  dbg=(trim(SRC%ID).eq.'SGN')
      NARG=0
      if (associated(SRC%DEF)) then
      !  if (SRC%DEF%ID.eq.'ENABLED') write(*,*) 'FIELD2ARRAY ENABLED ',SRC%LENGTH,IDX,SRC%DEF%TID
        L=SRC%LENGTH
        if (IDX.le.0) then
          imin=1
          imax=SRC%LENGTH
        else if (IDX.gt.SRC%LENGTH) then
          return
        else
          imin=IDX
          imax=IDX
        endif
        L=imax-imin+1
        if (L.LE.MA) then
          select case (SRC%DEF%TID)
          case(ID_FIELD_FLOAT)
            if (allocated(SRC%RFIELD)) then
              ARG(1:L)=SRC%RFIELD(imin:imax)
              NARG=L
            endif
          case(ID_FIELD_TABLE)
            if (allocated(SRC%RFIELD)) then
              ARG(1:L)=SRC%RFIELD(imin:imax)
            ! we read only LENGTH data, but report max. size, needed for compatibility
              NARG=SRC%DEF%DIM
            endif
          case(ID_FIELD_INT)
            if (allocated(SRC%IFIELD)) then
              ARG(1:L)=SRC%IFIELD(imin:imax)
              NARG=L
            endif
          case(ID_FIELD_ENUM)
            ARG(1)=SRC%EORD
            NARG=1
          case(ID_FIELD_SELECT)
        ! for selection list, return only selection index stoed as EORD
            ARG(1)=SRC%EORD
            NARG=1
          end select
        else
          write(*,*) trim(SRC%ID),' FIELD2ARRAY: buffer size too small for ',trim(SRC%DEF%ID),':',L,'>',MA
        endif
      else
        write(*,*) trim(SRC%ID),' FIELD2ARRAY: FIELD definition not associated'
      endif
      end subroutine FIELD2ARRAY

!---------------------------------------------------------
      SUBROUTINE ROW2ARRAY(SRC,IROW,IDX,ARG,MA,NARG)
! Convert table row to real array. IF IDX>0, set only given element
! Only for TABLE type
! return:
! ARG  ... output array
! NARG ... number of processed array elements
!---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IDX,MA,IROW
      real(KIND(1.D0)),intent(out) :: ARG(MA)
      integer,intent(out) :: NARG
      integer :: imin,imax,L,ncol
      NARG=0
      if (associated(SRC%DEF)) then
        ncol=SRC%DEF%NCOL
        if (IDX>ncol) return
        if (irow>SRC%MAXROW) return
        if (IDX<=0) then
          imin=(irow-1)*ncol+1
          imax=irow*ncol
        else
          imin=(irow-1)*ncol+IDX
          imax=imin
        endif
!1         format(a,6(1x,G12.5))
        !write(*,*) 'ROW2ARRAY ',irow,imin,imax
        L=imax-imin+1
        if (L<=MA) then
          select case (SRC%DEF%TID)
          case(ID_FIELD_TABLE)
            if (allocated(SRC%RFIELD)) then
              !write(*,1) 'ROW2ARRAY ',SRC%RFIELD(imin:imax)
              ARG(1:L)=SRC%RFIELD(imin:imax)
              NARG=L
            endif
          end select
        else
          write(*,*) trim(SRC%ID),' FIELD2ARRAY: buffer size too small for ',trim(SRC%DEF%ID),':',L,'>',MA
        endif
      else
        write(*,*) trim(SRC%ID),' FIELD2ARRAY: FIELD definition not associated'
      endif
      end subroutine ROW2ARRAY


C---------------------------------------------------------
      SUBROUTINE ARRAY2FIELD(TAR,IDX,ARG,MA,NARG)
! Convert real array to TFIELD object
! Only for types with numeric value (FLOAT,INT,ENUM)
! return:
! TAR  ... output TFIELD
! NARG ... number of read processed elements
C---------------------------------------------------------
      TYPE(TFIELD) :: TAR
      integer,intent(in) :: IDX,MA
      real(KIND(1.D0)),intent(in) :: ARG(MA)
      integer,intent(out) :: NARG
      integer :: imin,imax,L
      logical :: dbg = .false.
       ! dbg=(trim(TAR%ID).eq.'SGN')
      NARG=0
      if (associated(TAR%DEF)) then
        L=TAR%LENGTH
        if (IDX.le.0) then
          imin=1
          imax=TAR%LENGTH
        else if (IDX.gt.TAR%LENGTH) then
          return
        else
          imin=IDX
          imax=IDX
        endif
        L=imax-imin+1
        if (L.LE.MA) then
          select case (TAR%DEF%TID)
          case(ID_FIELD_FLOAT)
            if (allocated(TAR%RFIELD)) then
              TAR%RFIELD(imin:imax)=ARG(1:L)
              NARG=L
            endif
          case(ID_FIELD_TABLE)
            if (allocated(TAR%RFIELD)) then
              TAR%RFIELD(imin:imax)=ARG(1:L)
              NARG=TAR%DEF%DIM
            endif
          case(ID_FIELD_INT)
            if (allocated(TAR%IFIELD)) then
              TAR%IFIELD(imin:imax)=NINT(ARG(1:L))
              NARG=L
            endif
          case(ID_FIELD_ENUM)
            TAR%EORD=NINT(ARG(1))
            NARG=1
          case(ID_FIELD_SELECT)
        ! for selection list, update only selection index
            TAR%EORD=NINT(ARG(1))
            if (TAR%EORD.ge.TAR%NP) TAR%EORD=-1
            NARG=1
          end select
        else
          write(*,*) trim(TAR%ID),' ARRAY2FIELD: buffer size too small for ',trim(TAR%DEF%ID),':',L,'>',MA
        endif
        if (TAR%DEF%TID.ne.ID_FIELD_TABLE) then
          TAR%NP=max(TAR%NP,NARG)
        endif
      else
        write(*,*) trim(TAR%ID),' ARRAY2FIELD: FIELD definition not associated'
      endif

      end subroutine ARRAY2FIELD


C---------------------------------------------------------
      SUBROUTINE ARRAY2ROW(TAR,IROW,IDX,ARG,MA,NARG)
! Convert real array to  a row of TABLE filed
! return:
! TAR  ... output TFIELD
! NARG ... number of read processed elements
C---------------------------------------------------------
      TYPE(TFIELD) :: TAR
      integer,intent(in) :: IROW,IDX,MA
      real(KIND(1.D0)),intent(in) :: ARG(MA)
      integer,intent(out) :: NARG
      integer :: imin,imax,L,ncol
      logical :: dbg = .false.
       ! dbg=(trim(TAR%ID).eq.'SGN')
      NARG=0
      if (associated(TAR%DEF)) then
        ncol=TAR%DEF%NCOL
        if (IDX>ncol) return
        if (irow>TAR%MAXROW) return
        if (IDX<=0) then
          imin=(irow-1)*ncol+1
          imax=irow*ncol
        else
          imin=(irow-1)*ncol+IDX
          imax=imin
        endif
        L=imax-imin+1
        if (L.LE.MA) then
          select case (TAR%DEF%TID)
          case(ID_FIELD_TABLE)
            if (allocated(TAR%RFIELD)) then
              TAR%RFIELD(imin:imax)=ARG(1:L)
              NARG=L
            endif
          end select
        else
          write(*,*) trim(TAR%ID),' ARRAY2FIELD: buffer size too small for ',trim(TAR%DEF%ID),':',L,'>',MA
        endif
      else
        write(*,*) trim(TAR%ID),' ARRAY2FIELD: FIELD definition not associated'
      endif

      end subroutine ARRAY2ROW


C---------------------------------------------------------
      SUBROUTINE FIELD2STR(SRC,IDX,SARG,NARG)
C---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IDX
      character(*) :: SARG
      integer,intent(out) :: NARG
      character(256) :: CNUM
      integer :: imin,imax,LL,i
      character(LEN_NAME) :: SMI,SMA,SINC
      character(2*LEN_ID+1) :: SID
      logical :: dbg = .false.
      !  dbg=(trim(SRC%ID).eq.'SGN')
      NARG=0
      SARG=' '
      if (associated(SRC%DEF)) then
        LL=LEN(SARG)
        if (IDX.le.0) then
          imin=1
          imax=SRC%LENGTH
        else if (IDX.gt.SRC%LENGTH) then
          return
        else
          imin=IDX
          imax=IDX
        endif

        select case (SRC%DEF%TID)
        case(ID_FIELD_FLOAT)
          if (allocated(SRC%RFIELD)) then
            do i=imin,imax
              call FLOAT2STR(SRC%RFIELD(i),CNUM)
              if (len_trim(SARG)+len_trim(CNUM).lt.LL) then
                SARG=trim(SARG)//' '//trim(CNUM)
                NARG=NARG+1
              else
                exit
              endif
            enddo
          endif
        case(ID_FIELD_INT)
          if (allocated(SRC%IFIELD)) then
            do i=imin,imax
              call INT2STR(SRC%IFIELD(i),CNUM)
              if (len_trim(SARG)+len_trim(CNUM).lt.LL) then
                SARG=trim(SARG)//' '//trim(CNUM)
                NARG=NARG+1
              else
                exit
              endif
            enddo
          endif
        case(ID_FIELD_STR)
          if (allocated(SRC%SFIELD)) then
            do i=imin,imax
              if (len_trim(SARG)+len_trim(SRC%SFIELD(i)).lt.LL) then
                if (imin.eq.imax) then
                ! no | delimiter for single-value string arrays
                  SARG=trim(SRC%SFIELD(i))
                else
                  call APPENDPAR(SARG,'|',trim(SRC%SFIELD(i)))
                endif
                NARG=NARG+1
              else
                exit
              endif
            enddo
          endif
        case(ID_FIELD_ENUM)
          if (dbg) write(*,*) 'FIELD2STR ',trim(SRC%ID)
          call GetEnumString(SRC%DEF%ENUM,SRC%EORD,CNUM)
          if (len_trim(CNUM).lt.LL) then
            SARG=trim(CNUM)
            NARG=1
            if (dbg) write(*,*) '      strval=',trim(CNUM),' ord=',SRC%EORD
          endif
        case(ID_FIELD_SELECT)
        ! selection list: return selected element (note EORD indexing from 0)
          if (allocated(SRC%SFIELD).and.(SRC%EORD.le.SRC%LENGTH)) then
            SARG=trim(SRC%SFIELD(SRC%EORD+1))
            NARG=1
          else
            SARG='undefined'
          endif
        case(ID_FIELD_RANGE)
        ! only for array elements
          if (IDX.gt.0) then
            SID=SRC%VFIELD(idx)%ID
            call FLOAT2STR(1.D0*SRC%VFIELD(idx)%VMI,SMI)
            call FLOAT2STR(1.D0*SRC%VFIELD(idx)%VMA,SMA)
            call FLOAT2STR(1.D0*SRC%VFIELD(idx)%VINC,SINC)
            CNUM=trim(SID)//' '//trim(SMI)//' '//trim(SMA)//' '//trim(SINC)
            if (len_trim(SARG)+len_trim(CNUM).lt.LL) then
              SARG=trim(CNUM)
              NARG=1
            endif
          endif
        end select
      else
        write(*,*) trim(SRC%ID),' FIELD2STR: FIELD definition not associated'
      endif
      end SUBROUTINE FIELD2STR



C---------------------------------------------------------
      SUBROUTINE ROW2STR(SRC,IROW,IDX,SARG,NARG)
C---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IROW,IDX
      character(*) :: SARG
      integer,intent(out) :: NARG
      character(256) :: CNUM
      integer :: imin,imax,LL,i,ncol
      !character(LEN_NAME) :: SID,SMI,SMA,SINC
      NARG=0
      SARG=' '
      !write(*,*) 'ROW2STR ',SRC%NP,irow,IDX,SRC%MAXROW
      if (associated(SRC%DEF)) then
        LL=LEN(SARG)
        ncol=SRC%DEF%NCOL
        !write(*,*) 'ROW2STR ncol=',LL,ncol
        if (IDX>ncol) return
        if (irow>SRC%MAXROW) return
        if (IDX<=0) then
          imin=(irow-1)*ncol+1
          imax=irow*ncol
        else
          imin=(irow-1)*ncol+IDX
          imax=imin
        endif
        !write(*,*) 'ROW2STR ',irow,imin,imax
        select case (SRC%DEF%TID)
        case(ID_FIELD_TABLE)
          if (allocated(SRC%RFIELD)) then
            do i=imin,imax
              call FLOAT2STR(SRC%RFIELD(i),CNUM)
              if (len_trim(SARG)+len_trim(CNUM).lt.LL) then
                SARG=trim(SARG)//' '//trim(CNUM)
                NARG=NARG+1
              else
                exit
              endif
            enddo
          endif
        end select
      else
        write(*,*) trim(SRC%ID),' FIELD2STR: FIELD definition not associated'
      endif
      end SUBROUTINE ROW2STR


!---------------------------------------------------------
      SUBROUTINE STR2FIELD(SRC,IDX,SARG,NARG)
! store value passed as IDX and SARG in SRC (TFIELD)
! IDX ... array element index.
!         IDX<1 ... read the whole array
!         IDX>0 ... read corresponding array element
!---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IDX
      character(*) :: SARG
      integer,intent(out) :: NARG
      integer :: imin,imax,LL,i,na,IS,IL,ierr,IL1
      REAL(KIND(1.D0)) :: ARG(FIELD_BUFFER_FLOAT)
      character(FIELD_BUFFER_STR) :: SBUF
      type(TENUMDEF) :: EDEF
      character(LEN_LINE) :: LINE
      logical :: IsInteger,isEqual
      logical :: dbg = .false.
      !  dbg=(trim(SRC%ID).eq.'SGN')

      NARG=0
      if (associated(SRC%DEF)) then
        LL=LEN(SARG)
        if (IDX.le.0) then
    ! whole array
          imin=1
          imax=SRC%LENGTH
        else if (IDX.gt.SRC%LENGTH) then
    ! out of range => do nothing
          write(*,*) 'STR2FIELD IDX out of range',IDX,SRC%LENGTH
          return
        else
    ! single element
          imin=IDX
          imax=IDX
        endif
        select case (SRC%DEF%TID)
        case(ID_FIELD_FLOAT)
          if (allocated(SRC%RFIELD)) then
            call GETLINPARG(trim(SARG),ARG,FIELD_BUFFER_FLOAT,NA)
            if (NA.ne.(imax-imin+1)) then
              WRITE(LINE,"(a,' expected: ',I5,' actual: ',I5)") trim(SRC%DEF%ID),imax-imin+1,NA
              call MSG_ERROR('STR2FIELD.FLOAT','Wrong number of input parameters, '//trim(LINE),0,1)
              NARG=0
              return
            endif
            do i=1,NA
              SRC%RFIELD(imin+i-1)=ARG(i)
              NARG=NARG+1
            enddo
          endif
        case(ID_FIELD_INT)
          if (allocated(SRC%IFIELD)) then
            call GETLINPARG(trim(SARG),ARG,FIELD_BUFFER_FLOAT,NA)
            if (NA.ne.(imax-imin+1)) then
              WRITE(LINE,"(a,' expected: ',I5,' actual: ',I5)") trim(SRC%DEF%ID),imax-imin+1,NA
              call MSG_ERROR('STR2FIELD.INT','Wrong number of input parameters, '//trim(LINE),0,1)
              NARG=0
              return
            endif
            do i=1,NA
              SRC%IFIELD(imin+i-1)=NINT(ARG(i))
              NARG=NARG+1
            enddo
          endif
        case(ID_FIELD_STR)
          if (allocated(SRC%SFIELD)) then
          ! allow empty string on input, but convert to 'none'
            if (len_trim(SARG)==0) then
              !write(*,*) 'STR2FIELD ',imin,SRC%LENGTH
              NA=1
              if (SRC%LENGTH.ne.NA) then
                if (.not.ALLOC_FIELD(SRC,NA)) return
              endif
              SRC%SFIELD(imin)='none'
              NARG=NARG+1
              !write(*,*) 'STR2FIELD NARG=',NARG
            else
          ! reallocate field to exactly match the string aray size
              call COUNTPAR(trim(SARG),'|',NA)
              if (NA.ne.SRC%LENGTH) then
                if (.not.ALLOC_FIELD(SRC,NA)) return
              endif
            !if (NA.ne.(imax-imin+1)) then
            !  WRITE(LINE,"(a,' expected: ',I5,' actual: ',I5)") trim(SRC%DEF%ID),imax-imin+1,NA
            !  call MSG_ERROR('STR2FIELD.STR','Wrong number of input parameters, '//trim(LINE),0,1)
            !  NARG=0
            !  return
            !endif
          ! handle strings splitted in multiple lines
              do i=1,NA
                call FINDSTRPAR(SARG,'|',i,IS,IL)
                IL=min(IL,FIELD_BUFFER_STR)
                IL1=len_trim(SRC%SFIELD(imin+i-1))
                if (IL.gt.0) then
                  if (IL1>0) then
                    SRC%SFIELD(imin+i-1)=trim(SRC%SFIELD(imin+i-1))//SARG(IS:IS+IL-1)
                  else
                    SRC%SFIELD(imin+i-1)=SARG(IS:IS+IL-1)
                  endif
                else
                  if (IL1==0) SRC%SFIELD(imin+i-1)=' '
                endif
                NARG=NARG+1
              enddo
            endif
          endif
        case(ID_FIELD_ENUM)
    ! ENUM type: set EORD (selection index) according to the string value
    ! set EORD=-1 if value not defined
    ! Accept also ORD number as string
          call GetEnumDef(SRC%DEF%ENUM,EDEF)
          SBUF=trim(SARG)
        ! exceptions:
        ! accept other answers to NOYES type
          if (trim(EDEF%ID).eq.'NOYES') then
            select case(trim(SBUF))
              case('N','F','FALSE'); SBUF='NO'
              case('Y','T','TRUE'); SBUF='YES'
            end select
          endif
        ! defined string for SIGN is '+1', not '1', but '1' should be legal
          if ((trim(EDEF%ID).eq.'SIGN').and.(trim(SBUF).eq.'1')) SBUF='+1'
          call GetEnumOrd(SRC%DEF%ENUM,trim(SBUF),i)
        ! accept ORD number
          if ((i.lt.0).and.(IsInteger(trim(SBUF),i))) then
            if ((i.lt.0).or.(i.ge.EDEF%DIM)) i=-1
          endif
          SRC%EORD=i
          NARG=1
        case(ID_FIELD_SELECT)
    ! Selection list:
          call COUNTPAR(trim(SARG),'|',NA)
          if (NA.gt.1) then
      ! assumed input format is "%d|%s|%s|..." meaning EORD + array of list values
            call FINDSTRPAR(SARG,'|',1,IS,IL)
            SBUF=SARG(IS:IS+IL-1)
            if (IsInteger(trim(SBUF),i)) then
            ! reallocate field to exactly match the list size
              if ((NA-1).ne.SRC%LENGTH) then
                if (.not.ALLOC_FIELD(SRC,NA-1)) return
              endif
              if ((i.lt.0).or.(i.ge.SRC%LENGTH)) i=-1
              SRC%EORD=i
              do i=2,NA
                call FINDSTRPAR(SARG,'|',i,IS,IL)
                SBUF=SARG(IS:IS+IL-1)
              enddo
              SRC%NP=NA-1
            else
              call MSG_ERROR('STR2FIELD.SELECT','wrong LIST format for '//trim(SRC%ID)//' ('//trim(SARG)//')',0,1)
              return
            endif
      ! assumed input format is "%s" or "%d" ... string value or ORD nuber to be selected
          else if (NA.eq.1) then
            SRC%EORD=i-1
            isEqual=.false.
            i=0
            do while ((i.lt.SRC%LENGTH).and.(.not.isEqual))
              i=i+1
              isEqual=(trim(SARG).eq.trim(SRC%SFIELD(i)))
            !  write(*,*) 'STR2FIELD ',i,' [',trim(SARG),']',' [',trim(SRC%SFIELD(i)),']',isEqual
            enddo
            if (isEqual) then
              SRC%EORD=i-1
            else
          ! try to interpret SARG as ord number
              if (IsInteger(trim(SARG),i)) then
                if ((i.lt.0).or.(i.ge.SRC%LENGTH)) i=-1
                SRC%EORD=i
              endif
            endif
          !  write(*,*) 'STR2FIELD SELECT ',SRC%EORD,' [',trim(SARG),']',SRC%LENGTH,isEqual
          endif
          NARG=1
        case(ID_FIELD_RANGE)
        ! only for array elements
          if (IDX.gt.0) then
            call FINDSTRPAR(SARG,' ',1,IS,IL)
            if (IL.gt.2) then
              read(SARG(IS+IL:),*,IOSTAT=ierr) ARG(1),ARG(2),ARG(3)
              if (ierr.eq.0) then
                SRC%VFIELD(idx)%ID=SARG(IS:IS+IL-1)
                SRC%VFIELD(idx)%VMI=ARG(1)
                SRC%VFIELD(idx)%VMA=ARG(2)
                SRC%VFIELD(idx)%VINC=ARG(3)
                NARG=1
              endif
            endif
         ! otherwise read only one integer = number of range vector elements that should  follow
          else if (len_trim(SARG).gt.0) then
            read(SARG,*,IOSTAT=ierr) i
            if (ierr.eq.0) then
              SRC%NP=min(i,SRC%LENGTH)
              NARG=1
            else
              SRC%NP=0
            endif
          endif
          ! don't modify NP => return now
          return
        end select
        SRC%NP=max(SRC%NP,NARG)
      else
        write(*,*) trim(SRC%ID),' STR2FIELD: FIELD definition not associated'
      endif
      end SUBROUTINE STR2FIELD


!---------------------------------------------------------
      SUBROUTINE STR2ROW(SRC,IROW,IDX,SARG,NARG)
! store value passed as IDX and SARG in SRC (TFIELD)
! version for a table row
! IDX ... array element index.
!         IDX<1 ... read the whole array
!         IDX>0 ... read corresponding array element
!---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IROW,IDX
      character(*) :: SARG
      integer,intent(out) :: NARG
      integer :: imin,imax,LL,i,na,ncol,L
      REAL(KIND(1.D0)) :: ARG(FIELD_BUFFER_FLOAT)
      character(LEN_LINE) :: LINE
      logical :: dbg = .false.
      !  dbg=(trim(SRC%ID).eq.'SGN')

      NARG=0
      if (associated(SRC%DEF)) then
        LL=LEN(SARG)
        ncol=SRC%DEF%NCOL
        if (IDX>ncol) return
        if (irow>SRC%MAXROW) return
        if (IDX<=0) then
          imin=(irow-1)*ncol+1
          imax=irow*ncol
        else
          imin=(irow-1)*ncol+IDX
          imax=imin
        endif
        L=imax-imin+1
        !write(*,*) 'STR2ROW ',irow,imin,imax
        select case (SRC%DEF%TID)
        case(ID_FIELD_TABLE)
          if (allocated(SRC%RFIELD)) then
            call GETLINPARG(trim(SARG),ARG,FIELD_BUFFER_FLOAT,NA)
            if (NA.ne.L) then
              WRITE(LINE,"(a,' expected: ',I5,' actual: ',I5)") trim(SRC%DEF%ID),imax-imin+1,NA
              call MSG_ERROR('STR2ROW.FLOAT','Wrong number of input parameters, '//trim(LINE),0,1)
              NARG=0
              return
            endif
            !write(*,*) 'STR2ROW ',imin,ARG(1)
            do i=1,NA
              SRC%RFIELD(imin+i-1)=ARG(i)
              NARG=NARG+1
            enddo
          endif
        end select
        !SRC%NP=max(SRC%NP,NARG)
      else
        write(*,*) trim(SRC%ID),' STR2ROW: FIELD definition not associated'
      endif
      end SUBROUTINE STR2ROW


!---------------------------------------------------------
      SUBROUTINE ITEM2FIELD(SRC,IDX,SARG)
! store single item to RANGE and SELECT fields
!---------------------------------------------------------
      TYPE(TFIELD) :: SRC
      integer,intent(in) :: IDX
      character(*) :: SARG
      REAL(KIND(1.D0)) :: ARG(3)
      integer :: IS,IL,ierr
      if (associated(SRC%DEF).and.(IDX.gt.0).and.(IDX.le.SRC%LENGTH)) then
        select case (SRC%DEF%TID)
        case(ID_FIELD_SELECT)
          if (allocated(SRC%SFIELD)) then
            SRC%SFIELD(IDX)=trim(adjustl(SARG))
          !  write(*,*) 'ITEM2FIELD ',IDX,trim(SARG)
          else
            write(*,*) 'ITEM2FIELD not allocated SRC%LENGTH=',SRC%LENGTH
          endif
        case(ID_FIELD_RANGE)
          if (allocated(SRC%VFIELD)) then
            call FINDSTRPAR(SARG,' ',1,IS,IL)
            if (IL.gt.2) then
              read(SARG(IS+IL:),*,IOSTAT=ierr) ARG(1),ARG(2),ARG(3)
              if (ierr.eq.0) then
                SRC%VFIELD(idx)%ID=SARG(IS:IS+IL-1)
                SRC%VFIELD(idx)%VMI=ARG(1)
                SRC%VFIELD(idx)%VMA=ARG(2)
                SRC%VFIELD(idx)%VINC=ARG(3)
              endif
            endif
          endif
        end select
      else
        write(*,*) trim(SRC%ID),' ITEM2FIELD: FIELD definition not associated'
      endif
      end SUBROUTINE ITEM2FIELD


C---------------------------------------------------------
      SUBROUTINE FIELD_COPY(SRC,TAR,error)
C---------------------------------------------------------
      TYPE(TFIELD),intent(inout) :: SRC
      TYPE(TFIELD),intent(inout) :: TAR
      logical,intent(out) :: error
      integer :: L
      error=.true.
      if (associated(SRC%DEF).and.associated(TAR%DEF)) then
    !    write(*,*) 'FIELD_COPY owners=',SRC%DEF%OWNER,TAR%DEF%OWNER,' IDX:',SRC%DEF%IDX,TAR%DEF%IDX
        if ((SRC%DEF%OWNER.eq.TAR%DEF%OWNER).and.(SRC%DEF%IDX.eq.TAR%DEF%IDX)) then
            if (SRC%LENGTH.ne.TAR%LENGTH) then
              if (.not.ALLOC_FIELD(TAR,SRC%LENGTH)) then
                write(*,*) trim(SRC%ID),' FIELD_COPY: cannot allocate target '
                return
              endif
            endif
            L=SRC%LENGTH
            TAR%ID=SRC%ID
            TAR%NP=SRC%NP
            TAR%MAXROW=SRC%MAXROW
            select case (SRC%DEF%TID)
            case(ID_FIELD_FLOAT,ID_FIELD_TABLE)
              if ((allocated(SRC%RFIELD)).and.(allocated(TAR%RFIELD))) then
                TAR%RFIELD(1:L)=SRC%RFIELD(1:L)
                error=.false.
              endif
            case(ID_FIELD_INT)
              if ((allocated(SRC%IFIELD)).and.(allocated(TAR%IFIELD))) then
                TAR%IFIELD(1:L)=SRC%IFIELD(1:L)
                error=.false.
              endif
            case(ID_FIELD_STR)
              if ((allocated(SRC%SFIELD)).and.(allocated(TAR%SFIELD))) then
                TAR%SFIELD(1:L)=SRC%SFIELD(1:L)
                error=.false.
              endif
            case(ID_FIELD_RANGE)
              if ((allocated(SRC%VFIELD)).and.(allocated(TAR%VFIELD))) then
                TAR%VFIELD(1:L)=SRC%VFIELD(1:L)
                error=.false.
              endif
            case(ID_FIELD_ENUM)
              TAR%EORD=SRC%EORD
              error=.false.
            case(ID_FIELD_SELECT)
              TAR%EORD=SRC%EORD
              if ((allocated(SRC%SFIELD)).and.(allocated(TAR%SFIELD))) then
                TAR%SFIELD(1:L)=SRC%SFIELD(1:L)
                error=.false.
              endif
              error=.false.
            end select
        else
          write(*,*) trim(SRC%ID),' FIELD_COPY unequal types: [',trim(SRC%DEF%TID),',',trim(TAR%DEF%TID),']'
        endif
      else
        write(*,*) trim(SRC%ID),' FIELD_COPY: field definition not associated'
      endif
      end SUBROUTINE FIELD_COPY


!---------------------------------------------------------
      logical function FIELD_IS_TYPE(FIELD,TID)
! return true if the FIELD has definition associated and is of required type
!---------------------------------------------------------
      TYPE(TFIELD) :: FIELD
      character(*) :: TID
      logical :: LOG1
      LOG1=associated(FIELD%DEF)
      LOG1=(LOG1.and.(FIELD%DEF%TID.eq.trim(TID)))
      FIELD_IS_TYPE=LOG1
      end function FIELD_IS_TYPE


C---------------------------------------------------------
      logical function ALLOC_FIELD(FIELD,isize)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TFIELD) :: FIELD
      integer,intent(in) :: isize
      integer :: ierr,isz
      integer,parameter :: MAX_FIELD_DIM=256
      ierr=0
      isz=min(isize,MAX_FIELD_DIM)
      call DISPOSE_FIELD(FIELD)
      if (associated(FIELD%DEF)) then
        if (FIELD%DEF%TID.eq.ID_FIELD_TABLE) then
          ALLOC_FIELD=ALLOC_TABLE(FIELD,FIELD%DEF%NROW)
          return
        endif
        if (isz.gt.0) then
          select case (FIELD%DEF%TID)
            case(ID_FIELD_FLOAT)
              allocate(FIELD%RFIELD(1:isz),STAT=ierr)
              if (ierr.eq.0) FIELD%RFIELD(1:isz)=0.D0
            case(ID_FIELD_INT)
              allocate(FIELD%IFIELD(1:isz),STAT=ierr)
              if (ierr.eq.0) FIELD%IFIELD(1:isz)=0
            case(ID_FIELD_STR)
              allocate(FIELD%SFIELD(1:isz),STAT=ierr)
              if (ierr.eq.0) FIELD%SFIELD(1:isz)=' '
            case(ID_FIELD_RANGE)
              allocate(FIELD%VFIELD(1:isz),STAT=ierr)
            case(ID_FIELD_ENUM)
              isz=1
              FIELD%EORD=-1
            case(ID_FIELD_SELECT)
              FIELD%EORD=-1
              allocate(FIELD%SFIELD(1:isz),STAT=ierr)
              if (ierr.eq.0) FIELD%SFIELD(1:isz)=' '
            case default
              ierr=-2
          end select
          if (ierr.ne.0) then
            call MSG_ERROR('ALLOCATE_FIELD','cannot allocate for '//trim(FIELD%ID),0,1)
            FIELD%LENGTH=0
          ELSE
            FIELD%LENGTH=isz
            FIELD%ID=FIELD%DEF%ID
            FIELD%NP=FIELD%LENGTH
            if (trim(FIELD%DEF%TID).eq.ID_FIELD_RANGE) FIELD%NP=0
            if (trim(FIELD%DEF%TID).eq.ID_FIELD_SELECT) FIELD%NP=0
          endif
        endif
      else
        ierr=1
      endIF
      ALLOC_FIELD=(ierr.eq.0)
      end function ALLOC_FIELD



!---------------------------------------------------------
      logical function ALLOC_TABLE(FIELD,rows)
! allocate space for table
!---------------------------------------------------------
      TYPE(TFIELD) :: FIELD
      integer,intent(in) :: rows
      integer :: ierr,isz,nrow,ncol
      ierr=0
      ncol=FIELD%DEF%NCOL
      nrow=min(rows,FIELD%DEF%NROW)
      isz=ncol*nrow
      call DISPOSE_FIELD(FIELD)
      allocate(FIELD%RFIELD(1:isz),STAT=ierr)
      if (ierr.eq.0) FIELD%RFIELD(1:isz)=0.D0
      if (ierr.ne.0) then
        call MSG_ERROR('ALLOCATE_FIELD','cannot allocate for '//trim(FIELD%ID),0,1)
        FIELD%LENGTH=0
      ELSE
        FIELD%LENGTH=isz
        FIELD%ID=FIELD%DEF%ID
        FIELD%NP=nrow
        FIELD%MAXROW=nrow
      endIF
      ALLOC_TABLE=(ierr.eq.0)
      end function ALLOC_TABLE



C---------------------------------------------------------
      subroutine DISPOSE_FIELD(FIELD)
C clear buffers of the object
C---------------------------------------------------------
      TYPE(TFIELD) :: FIELD
      integer i
      integer :: ierr
      ierr=0
      FIELD%LENGTH=0
      FIELD%NP=0
      FIELD%MAXROW=0
      FIELD%EORD=-1
      if (.not.associated(FIELD%DEF)) return
      select case(FIELD%DEF%TID)
      case(ID_FIELD_FLOAT,ID_FIELD_TABLE)
        if (allocated(FIELD%RFIELD)) deallocate(FIELD%RFIELD,STAT=ierr)
      case(ID_FIELD_INT)
        if (allocated(FIELD%IFIELD)) deallocate(FIELD%IFIELD,STAT=ierr)
      case(ID_FIELD_STR)
        if (allocated(FIELD%SFIELD)) deallocate(FIELD%SFIELD,STAT=ierr)
      case(ID_FIELD_RANGE)
        if (allocated(FIELD%VFIELD)) deallocate(FIELD%VFIELD,STAT=ierr)
      case(ID_FIELD_SELECT)
        if (allocated(FIELD%SFIELD)) deallocate(FIELD%SFIELD,STAT=ierr)
      end select
      ! DO NOT remove association to DEF. It is used in subsequent FIELD allocation
      ! FIELD%DEF=>null()
      if (ierr.ne.0) then
        write(*,*) 'DISPOSE_FIELD err=',ierr,i,trim(FIELD%ID)
      endif
      end subroutine DISPOSE_FIELD


C---------------------------------------------------------
      subroutine CLEAR_FIELD(FIELD,IDX)
C clear field data without disposing memory
C---------------------------------------------------------
      TYPE(TFIELD) :: FIELD
      integer,intent(in) :: IDX
      integer :: imin,imax
      if (.not.associated(FIELD%DEF)) return
      if (IDX.le.0) then
    ! whole array
        imin=1
        imax=FIELD%LENGTH
      else if (IDX.gt.FIELD%LENGTH) then
    ! out of range => do nothing
        return
      else
    ! single element
        imin=IDX
        imax=IDX
      endif

      select case(FIELD%DEF%TID)
      case(ID_FIELD_FLOAT)
        if (allocated(FIELD%RFIELD)) FIELD%RFIELD(imin:imax)=0.D0
      case(ID_FIELD_INT)
        if (allocated(FIELD%IFIELD)) FIELD%IFIELD(imin:imax)=0
      case(ID_FIELD_STR)
        if (allocated(FIELD%SFIELD)) FIELD%SFIELD(imin:imax)=''
      end select
      end subroutine CLEAR_FIELD

      end module FIELDDATA
