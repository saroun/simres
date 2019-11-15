!//////////////////////////////////////////////////////////////////////
!////  $Id: dialogs.f,v 1.38 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.38 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Basic interactive dialogs for SIMRES console
!////
!//////////////////////////////////////////////////////////////////////
      MODULE DIALOGS
      use CLASSDEF
      USE FIELDDATA
      use VCTABLE
      use COMPONENTS_IO
      use XMLPARSE
      use XMLWRITER
      use IO
      implicit none

      private
      integer,parameter :: MARG=64
      real(KIND(1.D0)) :: DLG_ARG(MARG)
      integer :: DLG_NARG
      character(256) :: DLG_SARG

      logical :: QUIETMODE

      public DLG_STRING
      public DLG_INPUT
      public DLG_PARAM
      public DLG_INTEGER
  !    public QUIETMODE

      contains


C--------------------------------------------------------------
      SUBROUTINE DLG_INPUT(labels,values,idef)
C Dialog for numerical input.
C INPUT:
C   labels  ... a string with value names, items are delimited by :
C   values  ... if idef>0, should contain default values
C   idef    ... if >0, sprompt inlcudes default values accepted by <enter>
C RETURN:
C   values  ... real*8 array with return values
C NOTE! No check is made on values array dimension
C--------------------------------------------------------------
      CHARACTER*(*) labels
      REAL*8 values(*)
      INTEGER*4 idef
      CHARACTER*64 S,S1
      CHARACTER*128 sprompt
      INTEGER*4 IS,IL,IS1,IL1,ILP,IP,itry,ios
      character(LEN_LINE) LINE


1     format(a,' : ',$)
3     format(a)
4     format('invalid number format, ',$)
6     format('try again')
8     format('no input')

20    format(a)
21    format(a,' [',a,']')
22    format(G12.4)

      IL=1
      IP=0
      DO WHILE (IL.GE.0)
c// get next value name
        IP=IP+1
        CALL FINDSTRPAR(labels,':',IP,IS,IL)
c// format sprompt
        IF (IL.GT.0) THEN
          S=labels(IS:IS+IL-1)
          ILP=IL
        ELSE
          S='input number'
          ILP=12
        ENDIF
        IF (idef.gt.0) THEN
          write(S1,22) values(IP)
          CALL BOUNDS(S,IS1,IL1)
          write(sprompt,21) S(1:ILP),S1(IS1:IS1+IL1-1)
        ELSE
          write(sprompt,20) S(1:ILP)
        ENDIF
        ILP=LEN_TRIM(sprompt)
        IF (ILP.GT.128) ILP=128
        IF (IL.GE.0) THEN
c// read data and check validity
           itry=0
10         itry=itry+1
           WRITE(SOUT,1) sprompt(1:ILP)
           call RSXREAD(LINE)
           READ(LINE,3) S
           IL1=LEN_TRIM(S)
           IF (IL1.GT.0) READ(S,*,iostat=ios,err=11) values(IP)
c// validate input
11         IF (ios.ne.0) THEN ! format error
             write(sout,4)
           ENDIF
           IF (itry.lt.5.AND.ios.ne.0) then ! 5 attempts to enter a valid number
             write(sout,6)
             goto 10
           ELSE IF (ios.ne.0) then
             write(sout,8)
           ENDIF
        ENDIF
      ENDDO
      END SUBROUTINE DLG_INPUT


C--------------------------------------------------------------
      SUBROUTINE DLG_INTEGER(label,ivalue,idef,imin,imax)
C Single integer number input, with range checking.
C INPUT:
C   label   ... a string with value name
C   ivalue  ... if idef>0, should contain default value
C   idef    ... if >0, sprompt inlcudes the default value accepted by <enter>
C   imin,imax  ... limits (inclusive)
C RETURN:
C   ivalue  ... integer*4 return value
C--------------------------------------------------------------
      CHARACTER*(*) label
      INTEGER*4 ivalue,idef,imin,imax
      CHARACTER*64 S,S1,S2,S3
      character(LEN_LINE) :: LINE
      CHARACTER*128 sprompt
      INTEGER*4 IS1,IS2,IS3,IL,IL1,IL2,IL3,ILP,itry,ios


1     format(a,' : ',$)
3     format(a)
4     format('invalid number format, ',$)
5     format('value outside limits, ',$)
6     format('try again')
8     format('no input')

20    format(a,' (',a,' .. ',a,') ')
21    format(a,' (',a,' .. ',a,') [',a,']')
22    format(I8)

      sprompt=' '
c// format sprompt
      WRITE(S1,22) imin
      WRITE(S2,22) imax
      CALL BOUNDS(S1,IS1,IL1)
      CALL BOUNDS(S2,IS2,IL2)
      IL=LEN_TRIM(label)
      IF (IL.GT.0) THEN
        S=label(1:IL)
      ELSE
        S='input number'
        IL=12
      ENDIF
      IF (idef.gt.0) THEN
        write(S3,22) ivalue
        CALL BOUNDS(S3,IS3,IL3)
        write(sprompt,21) S(1:IL),S1(IS1:IS1+IL1-1),
     &      S2(IS2:IS2+IL2-1),S3(IS3:IS3+IL3-1)
      ELSE
        write(sprompt,20) S(1:IL),S1(IS1:IS1+IL1-1),
     &      S2(IS2:IS2+IL2-1)
      ENDIF
      ILP=LEN_TRIM(sprompt)
      IF (ILP.GT.128) ILP=128
c// read data and check validity
      itry=0
10    itry=itry+1
      WRITE(SOUT,1) sprompt(1:ILP)
      call RSXREAD(LINE)
      READ(LINE,3) S
      IL1=LEN_TRIM(S)
      ios=0
      IF (IL1.GT.0) READ(S,*,iostat=ios,err=11) ivalue
c// validate input
11    IF (ios.ne.0) THEN ! format error
          write(sout,4)
      ELSE IF (ivalue.LT.imin.OR.ivalue.GT.imax) THEN   ! range error
          ios=1
          write(sout,5)
      ENDIF
      IF (itry.lt.5.AND.ios.ne.0) then ! 5 attempts to enter a valid number
          write(sout,6)
          goto 10
      ELSE IF (ios.ne.0) then
          write(sout,8)
      ENDIF
      END SUBROUTINE DLG_INTEGER

C--------------------------------------------------------------
      SUBROUTINE DLG_DOUBLE(label,vnum,idef,dmin,dmax)
C Single real*8 number input, with range checking.
C INPUT:
C   label   ... a string with vnum name
C   vnum   ... if idef>0, should contain default vnum
C   idef    ... if >0, sprompt includes the default vnum accepted by <enter>
C   dmin,dmax  ... limits (inclusive)
C RETURN:
C   vnum  ... real*8 return vnum
C--------------------------------------------------------------
      CHARACTER*(*) label
      INTEGER*4 idef
      REAL*8 vnum,dmin,dmax
      CHARACTER*64 S1,S2,S3,SRANGE,SDEF
      CHARACTER(LEN_LINE) sprompt,LINE
      INTEGER*4 IL1,itry,ios


1     format(a,': ',$)
3     format(a)
4     format('invalid number format, ',$)
5     format('value outside limits, ',a,', ',$)
6     format('try again')
8     format('invalid input')

22    format('(',a,' .. ',a,')')
23    format('[',a,']')
24    format(a,' ',a,' ',a)

      sprompt=' '
      S1=' '
      S2=' '
      S3=' '
      SRANGE=' '
      SDEF=' '
c// format sprompt
      IF (idef.gt.0) THEN
        CALL FLOAT2STR(vnum,S3)
        write(SDEF,23) trim(S3)
      endif
      if (dmin.gt.-1.D10) then
        CALL FLOAT2STR(dmin,S1)
      endif
      if (dmax.lt.1.D10) then
        CALL FLOAT2STR(dmax,S2)
      endif
      if (len_trim(S1//S2).gt.0) then
        write(SRANGE,22) trim(S1),trim(S2)
      endif
      S1='input number'
      IF (label.ne.' ') S1=trim(label)
      write(sprompt,24) trim(S1),trim(SRANGE),trim(SDEF)
c// read data and check validity
      itry=0
10    itry=itry+1
      WRITE(SOUT,1) trim(sprompt)
      call RSXREAD(LINE)
      READ(LINE,3) S1
      IL1=LEN_TRIM(S1)
      IOS=0
      if (SRANGE.EQ.' ') then
        SRANGE='(-1E10 .. +1E10)'
      endif
      IF (IL1.GT.0) READ(S1,*,iostat=ios,err=11) vnum
c// validate input
11    IF (ios.ne.0) THEN ! format error
          write(sout,4)
      ELSE IF (vnum.LT.dmin.OR.vnum.GT.dmax) THEN   ! range error
          ios=1
          write(sout,5) trim(SRANGE)
      ENDIF
      IF ((itry.lt.5).AND.(ios.ne.0)) then ! 5 attempts to enter a valid number
          write(sout,6)
          goto 10
      ELSE IF (ios.ne.0) then
          write(sout,8)
      ENDIF
      END SUBROUTINE DLG_DOUBLE

C--------------------------------------------------------------
      SUBROUTINE DLG_ARRAY(label,ARR,NARR,dmin,dmax)
C Array real*8 number input, with range checking.
C INPUT:
C   NARR    ... array dimension
C   label   ... a string with value name
C   value   ... if idef>0, should contain default value
C   dmin,dmax  ... limits (inclusive)
C RETURN:
C   value  ... real*8 return value
C--------------------------------------------------------------
      CHARACTER*(*) label
      INTEGER*4 NARR
      REAL*8 ARR(NARR),dmin,dmax
      CHARACTER(LEN_LINE) LINE
      CHARACTER*64 S,S1,S2
      CHARACTER*128 sprompt
      INTEGER*4 I,IS1,IS2,IL,IL1,IL2,ILP,itry,ios


1     format(a,' : ',$)
3     format(a)
4     format('invalid format, ',$)
5     format('value outside limits, ',$)
6     format('try again')
8     format('no input')

20    format(a,' (',a,' .. ',a,')')
21    format(a,' (',a,' .. ',a,') [',a,']')
c22    format(F10.4)

      sprompt=' '
c// format sprompt
      CALL FLOAT2STR(dmin,S1)
      CALL FLOAT2STR(dmax,S2)
      CALL BOUNDS(S1,IS1,IL1)
      CALL BOUNDS(S2,IS2,IL2)
      IL=LEN_TRIM(label)
      IF (IL.GT.0) THEN
        S=label(1:IL)
      ELSE
        S='input array of numbers'
        IL=LEN_TRIM(S)
      ENDIF
      write(sprompt,20) S(1:IL),S1(IS1:IS1+IL1-1),S2(IS2:IS2+IL2-1)

      ILP=LEN_TRIM(sprompt)
      IF (ILP.GT.128) ILP=128
c// read data and check validity
      itry=0
10    itry=itry+1
      WRITE(SOUT,1) sprompt(1:ILP)
      call RSXREAD(LINE)
      IL1=LEN_TRIM(LINE)
      IOS=0
      IF (IL1.GT.0) READ(LINE,*,iostat=ios,err=11) ARR
c// validate input
11    IF (ios.ne.0) THEN ! format error
          write(sout,4)
      ELSE
          ios=0
          do i=1,narr
            if ((ARR(i).lt.dmin).or.(arr(i).gt.dmax)) ios=1
          enddo
          IF (ios.ne.0) write(sout,5)
      ENDIF
      IF ((itry.lt.5).AND.(ios.ne.0)) then ! 5 attempts to enter a valid number
          write(sout,6)
          goto 10
      ELSE IF (ios.ne.0) then
          write(sout,8)
      ENDIF
      END SUBROUTINE DLG_ARRAY

C--------------------------------------------------------------
      SUBROUTINE DLG_STRING(label,value,idef)
C Single string input
C INPUT:
C   label   ... a string with value name
C   value   ... if idef>0, should contain default value
C   idef    ... if >0, sprompt inlcudes the default value accepted by <enter>
C RETURN:
C   value  ...  return string value
C--------------------------------------------------------------
      CHARACTER*(*) label,value
      INTEGER*4 idef
      CHARACTER*128 S,sprompt
      CHARACTER(LEN_LINE) LINE
      INTEGER*4 IS1,IL,IL1,ILP,LL

1     format(a,' : ',$)
2     format(a,' [',a,']')
3     format(a)

      if (QUIETMODE) then
        call RSXREAD(LINE)
        LL=len_trim(LINE)
        if (LL.gt.0) then
          call BOUNDS(LINE,IS1,IL1)
          IL1=min(IL1,LEN(VALUE))
          VALUE=LINE(IS1:IS1+IL1-1)
        else
          VALUE=' '
        endif
        return
      endif
      sprompt=' '
c// format sprompt
      IL=LEN_TRIM(label)
      IF (IL.GT.0) THEN
        S=label(1:IL)
      ELSE
        S='input string'
        IL=12
      ENDIF
      if (idef.gt.0) then
        CALL BOUNDS(value,IS1,IL1)
        write(sprompt,2) S(1:IL),value(IS1:IS1+IL1-1)
      else
        sprompt=S(1:IL)
      endif
      ILP=LEN_TRIM(sprompt)
      IF (ILP.GT.128) ILP=128
c// read value
      WRITE(SOUT,1) sprompt(1:ILP)
      call RSXREAD(LINE)
      IL1=LEN_TRIM(LINE)
      IF (IL1.GT.LEN(value)) IL1=LEN(value) ! check for value size
      IF (IL1.GT.0) value=LINE(1:IL1)
      END SUBROUTINE DLG_STRING


!--------------------------------------------------------------
      subroutine PARSE_CMD(cmdline,CNAME,PNAME,ND,IX,IY,PDATA)
! parse command string and return
! CNAME   .. class ID or ?
! PNAME   .. parameter ID or ?
! ND .. dimension of the array
! IX,IY  ... indices of the element
! PDATA   .. the rest of line is interpreted as numerical or string data
C--------------------------------------------------------------
      CHARACTER(*),intent(in) :: cmdline
      CHARACTER(*),intent(out) :: CNAME,PNAME,PDATA
      integer,intent(out) :: ND,IX,IY
      integer :: IS,IL,LL
      character(LEN_LINE) :: LINE
      character(32) :: ID
      LINE=adjustl(cmdline)
      CNAME='?'
      PNAME='?'
      PDATA=' '
      IX=0
      IY=0
      ND=0
      LL=len_trim(LINE)
  ! read class name
      IS=1
      CALL FINDPAR(LINE,1,IS,IL)
      if (IL.le.0) return
      CNAME=LINE(IS:IS+IL-1)
      call MKUPCASE(CNAME)
  ! read parameter name
      IS=1
      CALL FINDPAR(LINE,2,IS,IL)
      if (IL.le.0) return
    ! try to intepret PNAME(i) or PNAME(i,j) as an array element
      call PARSE_ID(LINE(IS:IS+IL-1),ID,ND,IX,IY)
      call MKUPCASE(ID)
      PNAME=ID
      if (LL.gt.IS+IL) PDATA=trim(adjustl(LINE(IS+IL:)))
      END subroutine PARSE_CMD



!-----------------------------------------------------------------
      SUBROUTINE PARSE_ID(PNAME,ID,ND,IX,IY)
! Separate ID and indices from array ID of the form ID(IX) or ID(IX,IY)
! OUTPUT:
! ID   ... output string
! ND   ... array dimension (1 or 2)
! IX,IY ... dimensions
!-----------------------------------------------------------------
      CHARACTER(*) :: PNAME
      CHARACTER(*) :: ID
      integer,intent(out) :: ND,IX,IY
      character(LEN_LINE) :: SARG
      integer :: IA,IB,LP,LB,LC,LS,ie
      LP=len_trim(PNAME)
      IX=0
      IY=0
      ND=0
      SARG=''
      ie=0
      ID=trim(PNAME)
      !write(*,*) 'PARSE_ID (',trim(ID),')'
      if (LP.gt.3) then
        if (PNAME(LP:LP).eq.')') then
          LB=index(PNAME,'(')
          !write(*,*) 'PARSE_ID LB=',LB
          if (LB.gt.1) then
            if (LP.gt.LB+1) then
              SARG=PNAME(LB+1:LP-1)
              LS=len_trim(SARG)
              LC=index(SARG,',')
              !write(*,*) 'PARSE_ID (',trim(SARG),')',LC
              if (LC>1) then

                READ(SARG(1:LC-1),*,iostat=ie) IA
                if (ie==0) READ(SARG(LC+1:LS),*,iostat=ie) IB
                if (ie==0) then
                  IX=IA
                  IY=IB
                  ND=2
                endif
              else
                READ(SARG,*,iostat=ie) IA
                !write(*,*) 'PARSE_ID IA=',IA,ie
                if (ie==0) then
                  IX=IA
                  ND=1
                endif
              endif
            endif
            ID=PNAME(1:LB-1)
          endif
        endif
        !write(*,*) 'PARSE_ID (',trim(ID),') end: ',IX,IY
      endif

      end SUBROUTINE PARSE_ID

!--------------------------------------------------------------
      logical function DLG_GET_CLASS(isCMD,CNAME,PNAME,PDATA,IX,IY,OBJ)
! Get object reference from CNAME
! Interactive input is called if CNAME does not point to an existing object
! PNAME,PDATA,IA then return newly entered parameter name, data and array index
! if isCMD, look only for commands
! RETURN true if object is found
C--------------------------------------------------------------
      logical,intent(in) :: isCMD
      CHARACTER(*),intent(inout) :: CNAME,PNAME,PDATA
      integer,intent(inout) :: IX,IY
      type(TCLASS),intent(out) :: OBJ
      CHARACTER(LEN_LINE) :: LINE
      character(LEN_NAME) :: ITEMSTR
      character(64) :: PSTR
      INTEGER :: NLOOP,ND
      integer,parameter :: MAXLOOP=10
!      write(*,*) 'DLG_GET_CLASS [',trim(cmdline),']'

1     format('Give <',a,'> [parameter] [value], or "?", or "Q"')
4     FORMAT(a,' ',a,' does not exist')

! format prompt string
      if (isCMD) then
        ITEMSTR='command'
      else
        ITEMSTR='component'
      endif
      write(PSTR,1) trim(ITEMSTR)
! clear OBJ
      call PCLASS_CLEAR(OBJ)
      if (len_trim(CNAME).le.0) CNAME='?'
! try to identify OBJ
      NLOOP=0
      do while ((OBJ%ICLS.le.0).and.(NLOOP.lt.MAXLOOP))
        NLOOP=NLOOP+1
    ! if CNAME=?, get user input
        do while ((trim(CNAME).eq.'?').and.(.not.QUIETMODE))
          if (isCMD) then
            if (SINP.eq.5) call COMMANDS_PRN(SOUT)
          else
            if (SINP.eq.5) call BEAM_PRN(SOUT)
          endif
          CALL DLG_STRING(trim(PSTR),LINE,0)
          call PARSE_CMD(LINE,CNAME,PNAME,ND,IX,IY,PDATA)
          !call PARSE_CMD(LINE,CNAME,PNAME,IA,PDATA)
        enddo
  !      if (isCMD) write(*,*) 'DLG_GET_CLASS CNAME=[',trim(CNAME),']'
    ! evaluate input
        if (trim(CNAME).eq.'Q') exit
        if (isCMD) then
          OBJ=PCLASS(getCmdObjByID(CNAME))
  !        write(*,*) 'DLG_GET_CLASS OBJ=[',trim(OBJ%ID),OBJ%TCLS,OBJ%ICLS,']'
        else
          OBJ=getClassByID(CNAME)
        !  write(*,*) 'DLG_GET_CLASS OBJ=[',trim(OBJ%ID),OBJ%TCLS,OBJ%ICLS,']'
        endif
    ! OBJ not identified => try again interactively
        if (OBJ%ICLS.le.0) then
          write(sout,4) trim(ITEMSTR),trim(CNAME)
          CNAME='?'
        else
          CNAME=OBJ%ID
        endif
        if (QUIETMODE) NLOOP=MAXLOOP
      enddo
      DLG_GET_CLASS=(OBJ%ICLS.gt.0)
!      write(*,*) 'DLG_GET_CLASS [',trim(OBJ%ID),'] [',trim(inpline),'] ',DLG_GET_CLASS
      END function DLG_GET_CLASS

!--------------------------------------------------------------
      integer function DLG_GET_PNAME(OBJ,PNAME,PDATA,IX,IY)
! Get valid parameter ID for given object
! Interactive input is called if PNAME is not a valid parameter ID
! PNAME,PDATA,IA then return newly entered parameter name, data and array index
! RETURN parameter index in classes definition
C--------------------------------------------------------------
      type(TCLASS),intent(in) :: OBJ
      CHARACTER(*),intent(inout) :: PNAME,PDATA
      integer,intent(inout) :: IX,IY
      CHARACTER(LEN_LINE) :: LINE
      character(LEN_NAME) :: CNAME
      INTEGER :: NLOOP,IP,ICHILD,ND
      integer,parameter :: MAXLOOP=5
      logical :: dbg = .false.
5     FORMAT('parameter ',a,' is not defined for ',a,' ICLS=',I5)
      if (dbg) write(*,*) 'DLG_GET_PNAME [',trim(pname),']'
      IP=0
      DLG_GET_PNAME=0
      if (len_trim(PNAME).le.0) PNAME='?'
! try to identify parameter
      NLOOP=0
      do while ((IP.le.0).and.(NLOOP.lt.MAXLOOP))
        NLOOP=NLOOP+1
    ! if PNAME=?, get user input
        do while ((trim(PNAME).eq.'?').and.(.not.QUIETMODE))
          if (SINP.eq.5) CALL CLASS_PRN(SOUT,OBJ)
          LINE=' '
          CALL DLG_STRING('Give <parameter> [value], or "?", or "Q"',LINE,0)
         ! to separate index (i) from parameter name:
          call PARSE_CMD(trim(OBJ%ID)//' '//trim(LINE),CNAME,PNAME,ND,IX,IY,PDATA)
          !call PARSE_CMD(trim(OBJ%ID)//' '//trim(LINE),CNAME,PNAME,IA,PDATA)
        enddo
    ! evaluate input
        select case (trim(PNAME))
        case('.')
            NLOOP=MAXLOOP
        case('Q','q')
            NLOOP=MAXLOOP
            PNAME='Q'
        case default
          if (len_trim(PNAME).gt.0) then
            IP=GetFieldIndex(OBJ%ICLS,PNAME,ICHILD)
            if (IP.le.0) then
              write(sout,5) trim(PNAME),trim(OBJ%ID),OBJ%ICLS
              PNAME='?'
            endif
          endif
        end select
        if (QUIETMODE) NLOOP=MAXLOOP
        if (dbg)   READ(*,*)
      enddo
      DLG_GET_PNAME=IP
      if (dbg) write(*,*) 'DLG_GET_PNAME [',trim(PNAME),'] [',trim(PDATA),'] IX,IY=',IX,IY
      END function DLG_GET_PNAME



!--------------------------------------------------------
      subroutine list_enum(FIELD,IU)
! report class definition
!--------------------------------------------------------
      type(TFIELD) :: FIELD
      integer,intent(in) :: IU
      type(TENUMDEF) :: EDEF
      integer :: NP,j
      if (associated(FIELD%DEF).and.(FIELD%DEF%ENUM.gt.0)) then
        call GetEnumDef(FIELD%DEF%ENUM,EDEF)
        NP=EDEF%DIM
        do j=1,NP
          write(IU,"(2x,I4,5x,a)") j-1,trim(EDEF%ITEMS(j))
        enddo
      endif
      end subroutine list_enum

!--------------------------------------------------------------
      integer function READ_FIELD(FIELD,IX,PDATA)
!--------------------------------------------------------------
      type(TFIELD) :: FIELD
      integer,intent(in) :: IX
      CHARACTER(*),intent(inout) :: PDATA
      CHARACTER(LEN_LINE) :: LINE
      character(LEN_NAME) :: CNUM
      integer,parameter :: MAXLOOP=5
      integer :: NA,NEX,IRES,NLOOP,i
      NA=0
  ! get expected number of elements
      select case(trim(FIELD%DEF%TID))
        case(ID_FIELD_RANGE,ID_FIELD_ENUM,ID_FIELD_SELECT)
          NEX=1
        case default
          NEX=FIELD%DEF%DIM
      end select
      if (IX.gt.0) NEX=1
      IRES=0
      NLOOP=0
      !write(*,*) 'READ_FIELD ',trim(FIELD%ID),' ',IX
      do while((IRES.le.0).and.(NLOOP.lt.MAXLOOP))
        NLOOP=NLOOP+1
  ! try to read input from PDATA
        ! if (len_trim(PDATA).gt.0) then
          call CLEAR_FIELD(FIELD,IX)
          call STR2FIELD(FIELD,IX,PDATA,NA)
        ! endif
  ! if not successfult, get input interactively
        if (NA.lt.NEX) then
        if (XMLOUT==0) then
      ! read value from console
          select case(trim(FIELD%DEF%TID))
          case(ID_FIELD_RANGE)
            CALL DLG_STRING('Give number of range items, or "Q"',LINE,0)
          case(ID_FIELD_ENUM)
            call list_enum(FIELD,6)
            CALL DLG_STRING('Give ord. number, value, or "Q"',LINE,0)
          case default
            CALL DLG_STRING('Give <value>, or "Q"',LINE,0)
          end select
      ! interpret input
          select case (trim(LINE))
          case('.')
            IRES=-1
            NLOOP=MAXLOOP
          case('Q','q')
            IRES=-2
            NLOOP=MAXLOOP
          case default
            PDATA=trim(LINE)
          end select
        else
          ires=-1
          NLOOP=MAXLOOP
          write(*,*) 'READ_FIELD: wrong field data format, NA=',NA,' NEX=',NEX
        endif
        else
        ! read items of range list
          select case(trim(FIELD%DEF%TID))
          case(ID_FIELD_RANGE)
            NEX=FIELD%NP
            NA=0
            call INT2STR(FIELD%NP,CNUM)
            if (.not.QUIETMODE) write(*,*) 'give '//trim(CNUM)//' items in format: OBJECT.PARAM MIN MAX INCREMENT'
            do i=1,FIELD%NP
              call INT2STR(i,CNUM)
              CALL DLG_STRING('item '//trim(CNUM),LINE,0)
              call MKUPCASE(LINE)
              if (trim(LINE).eq.'Q') exit
              call STR2FIELD(FIELD,i,trim(LINE),ires)
              if (ires.le.0) then
                if (.not.QUIETMODE) call MSG_WARN('wrong input format for range item',1)
                exit
              endif
              NA=NA+1
            enddo
            ires=NA
            if (NA.ne.NEX) then
              NA=0
              NEX=1
              PDATA=' '
              ires=0
            endif
          case default
            ires=NA
          end select
        endif
        if (QUIETMODE) NLOOP=MAXLOOP
      enddo
      if (IRES.ne.NEX) IRES=0
      READ_FIELD=IRES
      end function READ_FIELD

!--------------------------------------------------------------
      integer function DLG_SET_VALUE(OBJ,PNAME,IX,IY,PDATA)
! Get parameter value(s) for OBJ%PNAME(IA) from PDATA string.
! The data from PDATA are stored in DLG_ARG and DLG_SARG fields of this module.
! Get values interactively if PDATA does not contain valid data.
! Sets new values to OBJ using PCLASS_INP
! RETURN
!   -2  ... QUIT ('Q' answer)
!   -1  ... step back to parameter selection dialog is required ('.' answer)
!    0  ... no valid data entered
!   >0  ... data where successfuly set to OBJ
C--------------------------------------------------------------
      type(TCLASS),intent(in) :: OBJ
      CHARACTER(*),intent(in) :: PNAME
      integer,intent(in) :: IX,IY
      CHARACTER(*),intent(inout) :: PDATA
      integer :: NA,IRES,IEXP
      TYPE(TFIELD) :: FIELD
      TYPE(TFIELDDEF),TARGET :: PDEF
      logical :: ValidValues
    !  character(FIELD_BUFFER_STR) :: SARG
    !  write(*,*) 'DLG_SET_VALUE ',trim(PNAME)
      IRES=0
      call GetFieldDef(OBJ%ICLS,PNAME,PDEF)
      if (PDEF%IDX.gt.0) then
        FIELD%DEF=>PDEF
      !  write(*,*) 'DLG_SET_VALUE PDEF=',PDEF%ID
        if (PDEF%TID.eq.ID_FIELD_TABLE) then
          if (IX>0) then
            if (ALLOC_TABLE(FIELD,PDEF%NROW)) then
              call PCLASS_OUT(OBJ,PNAME,FIELD,NA)
              call STR2ROW(FIELD,IX,IY,PDATA,NA)
              IEXP=1
              if (IY==0) IEXP=FIELD%DEF%NCOL
              ValidValues=(NA==IEXP)
              if (ValidValues) then
                call PCLASS_INP(OBJ,PNAME,FIELD,IRES)
                call BEAM_OUT ! update BPAR_ARRAY
              endif
              call DISPOSE_FIELD(FIELD)
            endif
          endif
        else if (ALLOC_FIELD(FIELD,PDEF%DIM)) then
          call PCLASS_OUT(OBJ,PNAME,FIELD,NA)
          IRES=READ_FIELD(FIELD,IX,PDATA)
          ValidValues=(IRES.gt.0)
          if (ValidValues) then
            call PCLASS_INP(OBJ,PNAME,FIELD,IRES)
            call BEAM_OUT ! update BPAR_ARRAY
          endif
          call DISPOSE_FIELD(FIELD)
        endif
      endif
      DLG_SET_VALUE=IRES
      END function DLG_SET_VALUE

!--------------------------------------------------------------
      integer function DLG_PARAM(cmdline,setpar)
! Get absolute index of an object parameter specified on a command line
! if setpar=true, try to read values from cmdline and set them to object
C--------------------------------------------------------------
      logical,intent(in) :: setpar
      CHARACTER(*),intent(in) :: cmdline
      type(TCLASS) :: OBJ
      CHARACTER(LEN_LINE) :: LINE,PDATA
      CHARACTER(LEN_NAME) :: CNAME,PNAME
      INTEGER :: NLOOP,IP,IX,IY,ND,IRES,LEVEL,IS,IL
      integer,parameter :: MAXLOOP=5
      logical :: QUIT,isCMD
      logical :: dbg=.false.

      if(dbg) write(*,*) 'DLG_PARAM [',trim(cmdline),']'
! default result
      DLG_PARAM=0
      PDATA=' '
      PNAME=' '
      CNAME=' '
      LINE=trim(adjustl(cmdline))
      isCMD=.false.

! search for 'CMD' string indicating command input
      IS=1
      CALL FINDPAR(cmdline,1,IS,IL)
      if (IL.eq.3) then
        CNAME=cmdline(IS:IS+IL-1)
        call MKUPCASE(CNAME)
        if (CNAME.eq.'CMD') then
          isCMD=.true.
          LINE=trim(adjustl(cmdline(IS+IL:)))
        endif
      endif

      call PARSE_CMD(LINE,CNAME,PNAME,ND,IX,IY,PDATA)
      !call PARSE_CMD(LINE,CNAME,PNAME,IA,PDATA)
    ! no interactive mode if cmdline contains full entry
      QUIETMODE=((setpar).and.(len_trim(PDATA).gt.0))
    ! no interactive mode if XMLOUT>0 (GUI mode)
      if (XMLOUT>0) QUIETMODE=.true.
      QUIT=.false.
      LEVEL=0
      NLOOP=0
      if(dbg) write(*,*) 'DLG_PARAM ',level,' [',trim(CNAME),'] [',trim(PNAME),'] [',trim(PDATA),']',QUIETMODE
! level=0 loop
      do while((LEVEL.eq.0).and.(NLOOP.lt.MAXLOOP).and.(.not.QUIT))
! get valid OBJ or die
        if(dbg) write(*,*) 'DLG_PARAM CNAME=','[',trim(CNAME),']'
        if (.not.DLG_GET_CLASS(isCMD,CNAME,PNAME,PDATA,IX,IY,OBJ)) return
        if(dbg) write(*,*) 'DLG_PARAM A',level,' [',trim(CNAME),'] [',trim(PNAME),'] [',trim(PDATA),']'
! level=1 loop: we have valid OBJ
        LEVEL=1
        do while ((LEVEL.eq.1).and.(NLOOP.lt.MAXLOOP).and.(.not.QUIT))
          NLOOP=NLOOP+1
      ! get valid PNAME
          IP=DLG_GET_PNAME(OBJ,PNAME,PDATA,IX,IY)
          if(dbg) write(*,*) 'DLG_PARAM B',level,' [',trim(PNAME),'] [',trim(PDATA),']',IX,IY
          if (IP.le.0) QUIT=.false.
      ! evaluate PNAME
          select case(trim(PNAME))
          case('Q') ; QUIT=.true.  ! quit
          case('.')           ! return one level back
            CNAME='?'
            LEVEL=0
          end select
! we have valid OBJ and PNAME => try to get values
          if (IP.gt.0) then
            if (.not.isCmd) DLG_PARAM=1
            if (setpar) then
              IRES=DLG_SET_VALUE(OBJ,PNAME,IX,IY,PDATA)
              if (dbg) write(*,*) 'DLG_PARAM C',level,' [',trim(PNAME),'] [',trim(PDATA),']',IX,IY
              select case(IRES)
              case(-2)
                QUIT=.true.   ! quit
              case(-1)
                CNAME='?' ! return one level back and get interactive input
                LEVEL=0
              case default    ! stay in this level and get interactive input
                PNAME='?'
                if (IRES.gt.0) NLOOP=0
              end select
            endif
          endif
          if (QUIETMODE) QUIT=.true.
        enddo
      enddo
      END function DLG_PARAM

      end MODULE DIALOGS