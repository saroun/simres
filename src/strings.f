C///////////////////////////////////////////////////////////////////////////
C////  $Id: strings.f,v 1.50 2019/08/16 17:16:25 saroun Exp $
C////
C////  strings.f  (c) J.Saroun, 1999-2011
C////
C////  String handling subroutines for RESTRAX
C////
C///////////////////////////////////////////////////////////////////////////

!-----------------------------------------------------------------
      SUBROUTINE DELNUL(LINE)
! replace all NULL characters with spaces
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) :: LINE
      integer :: i, L
      L = LEN(LINE)
      if (L>0) then
        do i=1,L
          if (ichar(LINE(i:i))==0) then
            LINE(i:i) = ' '
          endif
        enddo
      endif
      end SUBROUTINE DELNUL

!-----------------------------------------------------------------
      SUBROUTINE READ_ID_R8(LINE, ID, PARAM, MA, NA,IERR)
! Read REAL*8 number and its identifier from the LINE
! suppose format "NAME=number " or "NAME number "
! IERR=-1 ... wrong syntax
! IERR=-2 ... wrong number format
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: LINE
      CHARACTER*(*),intent(out) :: ID
      integer, intent(in) :: MA
      real(kind(1.D0)),intent(out) :: PARAM(MA)
      integer, intent(out) :: NA,IERR
      INTEGER :: LL,IEQ,INAME,LNAME,IVAL,LVAL,NP,ITAB
      real(kind(1.D0)) :: Z(64)
      CHARACTER*(16) :: STRID
      CHARACTER*(32) :: STRVAL
      NA=0
      NP=SIZE(PARAM)
      IERR=-1
      if (NP>64) then
        write(*,*) 'ERROR in READ_ID_R8: Can''t read more than 64 elements per line.'
        goto 98
      endif
      LL=len_trim(LINE)
      IEQ=index(LINE,"=")
      if ((IEQ<2).or.(LL<4)) goto 99
      ! no tabs after =
      ITAB=index(LINE,char(9))
      if (ITAB>IEQ) LL = ITAB-1
      CALL BOUNDS(LINE(1:IEQ-1),INAME,LNAME)
      STRID=LINE(INAME:INAME+LNAME-1)
      CALL BOUNDS(LINE(IEQ+1:LL),IVAL,LVAL)
      STRVAL=LINE(IVAL+IEQ:IVAL+IEQ+LVAL-1)

      call READ_ARRAY8(trim(STRVAL),64,Z,NA)
      NA=MIN(NA,NP)
      if (NA<=0) goto 98

      ! READ(trim(STRVAL),*,ERR=98) Z
      IERR=0
      !PARAM=Z
      PARAM(1:NA)=Z(1:NA)
      ID=trim(STRID)
      return
98    IERR=-2   ! cannot read value
99    RETURN
      END SUBROUTINE READ_ID_R8



C     -----------------------------------------------------------------
      SUBROUTINE READ_R8(NAME, LINE, RESULT, IERR)
C     Read REAL*8 number from the LINE identified by NAME
C     suppose format "NAME=number " or "NAME number "
C     -----------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) LINE,NAME
      INTEGER*4 L,ISTART,ILEN,I,INAME,LNAME,IERR,IS
      REAL*8 Z,RESULT
      CHARACTER*1 CH

      IERR=-1
      CALL BOUNDS(NAME,INAME,LNAME)
      IF (LNAME.LE.0) GOTO 99
      L=LEN_TRIM(LINE)
      IS=1
10    i=INDEX(LINE(IS:L),NAME(INAME:INAME+LNAME-1))
      if ((i.gt.1).AND.(LINE(is+i-2:is+i-2).NE.' ')) THEN  ! space delimiter must precede the name
        IS=IS+i-1+LNAME ! try other occurences after this one
        GOTO 10
      ENDIF
C// identifier NAME found
      IF(i.GT.0) THEN
        i=i+IS-1+LNAME  ! i=1st character after NAME
C// name must be followed by = or space
        IF (LINE(i:i).NE.'='.AND.LINE(i:i).NE.' ') RETURN
        i=i+1
        ISTART=1
        CALL FINDPAR(LINE(i:L),1,ISTART,ILEN) ! find next substring
        ISTART=ISTART+i-1
        CH=LINE(ISTART:ISTART)
C// skip the '=' character
        IF ((ILEN.GE.1).AND.(CH.EQ.'=')) ISTART=ISTART+1
C// exclude T and F, which might be interpreted as a valid number on some systems
        IF (CH.EQ.'T'.OR.CH.EQ.'F') RETURN
C// try to read number from the rest of the line
        READ(LINE(ISTART:L),*,ERR=99) Z
        IERR=0
        RESULT=Z
c        write(*,*) LINE(ISTART:L)//'=',Z
c        pause
      ENDIF

      RETURN
99    IERR=-2   ! cannot read value
      END SUBROUTINE READ_R8

C     -----------------------------------------------------------------
      SUBROUTINE READ_I4(NAME,LINE,RESULT,IERR)
C     Read INTEGER*4 number from the LINE identified by NAME
C     suppose format "NAME=number " or "NAME number "
C     -----------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER*(*) LINE,NAME
      INTEGER*4 L,ISTART,ILEN,I,INAME,LNAME,IERR,IS
      INTEGER*4 Z,RESULT
      CHARACTER*1 CH

      IERR=-1
      CALL BOUNDS(NAME,INAME,LNAME)
      IF (LNAME.LE.0) GOTO 99
      L=LEN_TRIM(LINE)
      IS=1
10    i=INDEX(LINE,NAME(INAME:INAME+LNAME-1))
      if (i.gt.1) THEN  ! space delimiter must precede the name
            if (LINE(is+i-2:is+i-2).NE.' ') then
                  IS=IS+i-1+LNAME ! try other occurences after this one
                  GOTO 10
            endif
      ENDIF
C// identifier NAME found
      IF(i.GT.0) THEN
            i=i+IS-1+LNAME  ! i=first character after NAME
C// name must be followed by = or space
            IF (LINE(i:i).NE.'='.AND.LINE(i:i).NE.' ') RETURN
            i=i+1
            ISTART=1
            CALL FINDPAR(LINE(i:L),1,ISTART,ILEN) ! find next substring
            ISTART=ISTART+i-1
            CH=LINE(ISTART:ISTART)
C// skip the '=' character
            IF ((ILEN.GE.1).AND.(CH.EQ.'=')) ISTART=ISTART+1
C// try to read a number from the rest of the line
            READ(LINE(ISTART:L),*,ERR=99) Z
            IERR=0
            RESULT=Z
      ENDIF

      RETURN
99    IERR=-2   ! cannot read value
      END SUBROUTINE READ_I4

C     -----------------------------------------------------------------
      SUBROUTINE READ_STR(NAME,LINE,RSLT,IERR)
C     Read STRING from the LINE identified by NAME
C     suppose format "NAME=string"
C     -----------------------------------------------------------------

      IMPLICIT NONE
      CHARACTER(LEN=*) LINE,NAME
      INTEGER :: L,ISTART,ILEN,I,INAME,LNAME,IERR,LR,IS
      CHARACTER(LEN=*) RSLT
      CHARACTER(LEN=1) :: CH
      IERR=-1
      CALL BOUNDS(NAME,INAME,LNAME)
      L=LEN_TRIM(LINE)
      IF (LNAME.LE.0.OR.L.LE.0) GOTO 99
      LR=LEN(RSLT)
      IS=1
10    i=INDEX(LINE(IS:L),NAME(INAME:INAME+LNAME-1)//'=')
 ! space delimiter must precede the name
      if (i.gt.1) THEN 
            if (LINE(is+i-2:is+i-2).NE.' ') then
                  IS=IS+i-1+LNAME ! try other occurences after this one
                  GOTO 10
            endif
      ENDIF
! identifier NAME found
      IF(i.GT.0) THEN
        i=i+IS-1+LNAME  ! i=first character after NAME
! name must be followed by = or space
        IF (LINE(i:i).NE.'='.AND.LINE(i:i).NE.' ') goto 98
        i=i+1
        ISTART=1
        CALL FINDPAR(LINE(i:L),1,ISTART,ILEN) ! find next substring
        ISTART=ISTART+i-1
        CH=LINE(ISTART:ISTART)
C// skip the '=' character
        IF((ILEN.GT.0).AND.(CH.EQ.'=')) ISTART=ISTART+1
        IF(ILEN.GT.LR) ILEN=LR
c        write(*,*) LINE(1:i-1)//'='//LINE(ISTART:ISTART+ILEN-1)
        IF(ILEN.GT.0) THEN
C// RSLT is all after the delimiter (=)
          RSLT=LINE(ISTART:ISTART+ILEN-1)
        ELSE
          RSLT=' '
        ENDIF
        IERR=0
      ELSE
      ENDIF
      RETURN

98    IERR=-2   ! cannot read value
      return

99    IERR=-3   ! cannot read value
      END



C     ---------------------------------------------------
      CHARACTER*(*) FUNCTION STRIP(line)
C     removes initial and terminal spaces
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) line
      INTEGER*4 L,I
      L=LEN(line)
      IF (L.EQ.0) THEN
         STRIP=' '
         RETURN
      ENDIF
      DO WHILE (line(L:L).EQ.' '.AND.L.GT.0)
         L=L-1
      ENDDO
      I=1
      DO WHILE (line(I:I).EQ.' '.AND.I.LT.L)
         I=I+1
      ENDDO
      IF (I.GT.L) THEN
         STRIP=' '
      ELSE
         STRIP=line(I:L)
      ENDIF
      END

C     ---------------------------------------------------
      CHARACTER*(*) FUNCTION CONCAT(str1,str2)
C     connects 2 strings without spaces in between
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) str1,str2
      INTEGER*4 I1,I2,J1,J2
      CALL BOUNDS(str1,I1,J1)
      CALL BOUNDS(str2,I2,J2)
      CONCAT=str1(I1:I1+J1-1)//str2(I2:I2+J2-1)
      END

C     ---------------------------------------------------
      SUBROUTINE SCPY(str1,is,il,str2)
C     return str1(is:is+il-1) in str2, with range checking
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) str1,str2
      INTEGER*4 is,il,L
      IF (IL.GT.0) THEN
        L=MIN(il,LEN(str2))
        str2=str1(is:is+L-1)
      ELSE
        str2=' '
      ENDIF
      END


C     ---------------------------------------------------
      SUBROUTINE BOUNDS(line,ISTART,ILEN)
C     get position of string by stripping off the surrounding spaces
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) line
      INTEGER*4 L,I,ISTART,ILEN
      L=LEN_TRIM(line)
      ISTART=1
      IF (L.EQ.0) THEN
         ILEN=0
      else
        DO WHILE ((L.GT.0).AND.(line(L:L).EQ.CHAR(9)))
          L=L-1
        ENDDO
      endif
      IF(L.gt.0) then
        I=1
        DO WHILE (I.LE.L.AND.(line(I:I).EQ.' '.OR.line(I:I).EQ.CHAR(0).OR.line(I:I).EQ.CHAR(9)))
           I=I+1
        ENDDO
        IF (I.GT.L) THEN
           ILEN=0
        ELSE
           ISTART=I
           ILEN=L-I+1
        ENDIF
      else
        ILEN=0
      ENDIF
      END

c     ---------------------------------------------------
      SUBROUTINE STRCOMPACT(line,ISTART,ILEN)
C     remove multiplied spaces, return start (ISTART) and length (ILEN) of the resulting string
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) line
      INTEGER*4 L,I,J,ISTART,ILEN
      L=LEN_TRIM(line)
      I=1
      DO WHILE (I.LT.L)
        IF (LINE(I:I+1).EQ.'  ') THEN
           DO J=I+1,L-1
             LINE(J:J)=LINE(J+1:J+1)
           ENDDO
           LINE(L:L)=' '
           L=L-1
        ELSE
           I=I+1
        ENDIF
      ENDDO
      CALL BOUNDS(line,ISTART,ILEN)

      END

C     -------------------------------------------------------------------
      INTEGER*4 FUNCTION TRUELEN(line)
C     return true length of the string
C     (without trailing spaces or NULL characters)
C     -------------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) line
      INTEGER*4 L
      TRUELEN=0
      L=LEN(line)
      IF (L.EQ.0) THEN
         TRUELEN=0
         RETURN
      ENDIF
      DO WHILE (L.GT.0. AND.(line(L:L).EQ.' '.OR.line(L:L).EQ.char(0)))
         L=L-1
      ENDDO
      TRUELEN=L
      END

C     ---------------------------------------------------
      SUBROUTINE WRITELINE(LINE,IU)
C     writes string without surrounding spaces
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER*4 L,IS,IU
      CALL BOUNDS(LINE,IS,L)
1     FORMAT(a)
      WRITE(IU,1) LINE(IS:IS+L-1)
      END


C     ---------------------------------------------------
      CHARACTER*(*) FUNCTION UPCASE(LINE)
C     converts LINE to uppercase
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER*4 L,I
      L=LEN_TRIM(LINE)
      DO i=1,L
         IF((LINE(i:i).GE.'a').AND.(LINE(i:i).LE.'z')) THEN
            UPCASE(i:i)=CHAR(ICHAR(LINE(i:i))-32)
         ELSE
            UPCASE(i:i)=LINE(i:i)
         ENDIF
      ENDDO
      END

C     ---------------------------------------------------
      SUBROUTINE MKUPCASE(LINE)
C     converts LINE to uppercase
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      INTEGER*4 L,I
      L=LEN_TRIM(LINE)
      DO i=1,L
         IF((LINE(i:i).GE.'a').AND.(LINE(i:i).LE.'z')) THEN
            LINE(i:i)=CHAR(ICHAR(LINE(i:i))-32)
         ENDIF
      ENDDO
      END

!-----------------------------------------------------
      integer function IDCOMPARE(STR1,STR2)
! compares two ID strings, case non sensitive
! compares max. 16 characters
! RETURN
!    -1 .. not equal or substring
!     0 .. exact match
!     1 .. STR1 matches the initial part of STR2
!-----------------------------------------------------
      character(*),intent(in) :: STR1,STR2
      character(32) :: CSTR1,CSTR2
  !    write(*,*) 'IDCOMPARE: ',trim(STR1),' ',trim(STR2)
      CSTR1=adjustl(STR1)
      CSTR2=adjustl(STR2)
  !    write(*,*) '           ',trim(STR1),' ',trim(STR2)
      call MKUPCASE(CSTR1)
      call MKUPCASE(CSTR2)
      if (trim(CSTR1).eq.trim(CSTR2)) then
        IDCOMPARE=0
      else if (INDEX(CSTR2,trim(CSTR1)).eq.1) then
        IDCOMPARE=1
      else
        IDCOMPARE=-1
      endif
      end function IDCOMPARE

!     ---------------------------------------------------
      LOGICAL*4 FUNCTION IsFloat(LINE,RNUM)
!     true if LINE is a number
!     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      REAL(KIND(1.D0)),intent(out) :: RNUM
      INTEGER :: L,IERR
      L=LEN_TRIM(LINE)
      IERR=0
      IF ((L.LE.0).OR.(INDEX('1234567890.+-',LINE(1:1)).le.0)) IERR=1
      if (IERR.eq.0) READ(LINE(1:L),*,IOSTAT=IERR)  RNUM
      IsFloat=(IERR.EQ.0)
      END FUNCTION IsFloat

!     ---------------------------------------------------
      LOGICAL*4 FUNCTION IsInteger(LINE,INUM)
!     true if LINE is an integer
!     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      integer,intent(out) :: INUM
      INTEGER :: L,IERR
      L=LEN_TRIM(LINE)
      IERR=0
      IF ((L.LE.0).OR.(INDEX('1234567890.+-',LINE(1:1)).le.0)) IERR=1
      if (IERR.eq.0) READ(LINE(1:L),*,IOSTAT=IERR)  INUM
      IsInteger=(IERR.EQ.0)
    !  write(*,*) 'IsInteger ',LINE(1:L),' ',L,INUM,IsInteger
      END FUNCTION IsInteger

!-----------------------------------------------------------------
      SUBROUTINE READ_ARRAY4(LINE,A,MA,NA)
! Read REAL*4 array A(MA) from the LINE
! NA is actually read number of data found on line
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: LINE
      INTEGER,intent(in) :: MA
      INTEGER,intent(out) :: NA
      REAL,intent(out) :: A(MA)
      REAL(KIND(1.D0)) :: Z
      INTEGER*4 IS,IL
      logical :: IsFloat
      IS=1
      NA=0
      CALL FINDPAR(LINE,1,IS,IL)
      DO WHILE ((IL.GT.0).and.(NA.lt.MA))
        if (IsFloat(LINE(IS:IS+IL-1),Z)) THEN
          NA=NA+1
          A(NA)=Z
          IS=IS+IL
          call FINDPAR(LINE,1,IS,IL)
        else
          exit
        endif
      ENDDO
      end SUBROUTINE READ_ARRAY4


!-----------------------------------------------------------------
      SUBROUTINE STR2ARRAY4(LINE,DEL,A,MA,NA)
! Read REAL*4 array A(MA) from the LINE
! NA is actually read number of data found on line
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: LINE
      character(1),intent(in) :: DEL
      INTEGER,intent(in) :: MA
      INTEGER,intent(out) :: NA
      REAL,intent(out) :: A(MA)
      REAL(KIND(1.D0)) :: Z
      INTEGER*4 IS,IL,IP,IS1,IL1,IS0
      logical :: IsFloat
      IS=1
      NA=0
      IP=1
      call BOUNDS(LINE,IS0,IL)
      CALL FINDSTRPAR(LINE(IS0:),DEL,IP,IS,IL)
      IS=IS+IS0-1
      DO WHILE ((IL.GT.0).and.(NA.lt.MA))
        call BOUNDS(LINE(IS:IS+IL-1),IS1,IL1)
        IS1=IS+IS1-1
        if (IsFloat(LINE(IS1:IS1+IL1-1),Z)) THEN
          NA=NA+1
          A(NA)=Z
          IP=IP+1
          call FINDSTRPAR(LINE(IS0:),DEL,IP,IS,IL)
          IS=IS+IS0-1
        else
          exit
        endif
      ENDDO
      end SUBROUTINE STR2ARRAY4

!-----------------------------------------------------------------
      SUBROUTINE READ_ARRAY8(LINE,MA,A,NA)
! Read REAL*8 array A(MA) from the LINE
! NA is actually read number of data found on line
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: LINE
      integer, intent(in) :: MA
      REAL(KIND(1.D0)),intent(inout) :: A(MA)
      INTEGER,intent(out) :: NA
      REAL(KIND(1.D0)) :: Z
      INTEGER :: IS,IL
      logical :: IsFloat
      IS=1
      NA=0
      CALL FINDPAR(LINE,1,IS,IL)
      DO WHILE ((IL.GT.0).and.(NA.lt.MA))
        if (IsFloat(LINE(IS:IS+IL-1),Z)) THEN
          NA=NA+1
          A(NA)=Z
          IS=IS+IL
          call FINDPAR(LINE,1,IS,IL)
        else
          exit
        endif
      ENDDO
      end SUBROUTINE READ_ARRAY8


!-----------------------------------------------------------------
      SUBROUTINE READ_STRARRAY(LINE,DLM,STRARRAY,NA)
! Read space-delimited strings from the LINE
! and convert it into a string delimited by DLM
! NA is the number of substrings (array elements) found in line
! e.g. LINE="a dfg vbn" and DLM="|" gives STRARRAY="a|dfg|vbn" and NA=3
!-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*),intent(in) :: LINE
      character*1,intent(in) :: DLM
      CHARACTER*(*),intent(out) :: STRARRAY
      INTEGER,intent(out) :: NA
      INTEGER*4 IS,IL,LL,LS
      IS=1
      NA=0
      CALL FINDPAR(LINE,1,IS,IL)
      LL=LEN(STRARRAY)
      LS=0
      STRARRAY=''
      DO WHILE ((IL.GT.0).and.(LS.lt.LL))
        LS=len_trim(STRARRAY)
        if ((LS+IL+1).gt.LL) then
          EXIT
        else if (LS.eq.0) then
          STRARRAY=LINE(IS:IS+IL-1)
        else
          STRARRAY=trim(STRARRAY)//DLM//LINE(IS:IS+IL-1)
        endif
        NA=NA+1
        IS=IS+IL
        call FINDPAR(LINE,1,IS,IL)
      ENDDO
      end SUBROUTINE READ_STRARRAY


!---------------------------------------------------
      SUBROUTINE GETLINPARG(LINE,ARG,MA,NA)
! Get numerical arguments from the line and store them in ARG
! Takes only consecutive series of space-delimited numbers,
! stops enumerating at a non-number parameter or EOL
!---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      integer,intent(in) :: MA
      INTEGER,intent(out) :: NA
      REAL(KIND(1.D0)),intent(out) :: ARG(MA)
      ARG(1:MA)=0.D0
      !if (MA>1) write(*,*) 'GETLINPARG ',trim(LINE),' MA=', MA
      call READ_ARRAY8(LINE,MA,ARG,NA)
      END SUBROUTINE GETLINPARG

!     -----------------------------------------------------------------
      logical function isIndent(CH)
!     -----------------------------------------------------------------
       character*1 CH
        isIndent=((CH.EQ.' ').or.(CH.EQ.CHAR(9)))
      end function isIndent

!     -----------------------------------------------------------------
      SUBROUTINE FINDPAR(LINE,IPAR,ISTART,ILEN)
!     Finds IPAR-th parameter found on LINE, starting from
!     ISTART-th character
!     returns starting position in ISTART and length in ILEN
!     -----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      integer,intent(in) :: IPAR
      integer, intent(out) :: ILEN
      integer, intent(INout) :: ISTART
      INTEGER:: L,i1,i2,K,INQ
      logical :: isIndent
      logical :: notIndent
  !    write(*,*) 'FINDPAR [',trim(LINE),'] ipar=',IPAR
      ILEN=0
      i1=MAX(1,ISTART)
      L=LEN_TRIM(LINE)
      if (L.LT.i1) then ! no parameters after ISTART
        ISTART=0
        return
      endif
      i2=i1
      DO K=1,IPAR
! skip initial spaces and tabs
        DO WHILE ((i1.le.L).AND.isIndent(LINE(i1:i1)))
          i1=i1+1
        END DO
        i2=i1
        if (i1.gt.L) exit
        INQ=0
      !  write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
! find i2 = (last non-space character after i1) + 1
        DO WHILE (i2.le.L)
            !if(i2==29) then
            !      WRITE(*,*) 'FINDPAR ['//TRIM(line)//'] ',i1, i2, L
            !endif
            notIndent = (.NOT.isIndent(LINE(i2:i2)))
            if (notIndent.or.INQ.EQ.1) then
                  If (LINE(i2:i2).EQ.'"'.AND.i1.eq.i2) INQ=1 ! first quotes
                  if (LINE(i2:i2).EQ.'"'.AND.i1.LT.i2) INQ=0 ! second quotes
                  i2=i2+1
            ELSE
                  EXIT
            endif
        END DO
      !  write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
        IF (LINE(i1:i1).EQ.'"') i1=i1+1
        IF ((i2.gt.1).and.(LINE(i2-1:i2-1).EQ.'"')) i2=i2-1
      !  write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
      !  write(*,*) k,' FINDPAR ',LINE(i2:i2)
        IF (K.LT.IPAR) then
          if (i2.gt.L) then
            i2=i1
            exit
          else if (LINE(i2:i2).EQ.'"') then
            I1=I2+1 ! skip quotes
          else
            I1=I2 ! try next parameter
          endif
          if (I1.gt.L) then
            i2=i1
            exit
          endif
        endif
      ENDDO
      ILEN=I2-I1
      ISTART=I1
  !    write(*,*) k,' FINDPAR result=',I1,I2,' [',LINE(ISTART:ISTART+ILEN-1),']'
      END SUBROUTINE FINDPAR

!     -----------------------------------------------------------------
      SUBROUTINE FINDPAR_NEW(LINE,IPAR,ISTART,ILEN)
!     Finds IPAR-th parameter found on LINE, starting from
!     ISTART-th character
!     returns starting position in ISTART and length in ILEN
!     -----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      integer,intent(in) :: IPAR
      integer, intent(out) :: ISTART,ILEN
      INTEGER*4 L,i1,i2,K,INQ
1     format(I3,$)
      write(*,*) 'FINDPAR [',trim(LINE),'] ipar=',IPAR
      ILEN=0
      i1=MAX(1,ISTART)
      L=LEN_TRIM(LINE)
      if (L.LT.i1) then ! no parameters after ISTART
        ISTART=0
        return
      endif
      i2=i1
      DO K=1,IPAR
! skip initial spaces and tabs
        DO WHILE ((i1.le.L).AND.(LINE(i1:i1).EQ.' '))
        ! DO WHILE ((i1.le.L).AND.(LINE(i1:i1).EQ.' '.or.LINE(i1:i1).EQ.achar(9)))
          i1=i1+1
        END DO
        i2=i1
        if (i1.gt.L) exit
        INQ=0
        write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
! find i2 = (last non-space character after i1) + 1
        DO WHILE ((i2.le.L).AND.(LINE(i2:i2).NE.' '.OR.INQ.EQ.1))
        ! DO WHILE ((i2.le.L).AND.(LINE(i2:i2).NE.' '.or.LINE(i2:i2).NE.achar(9).OR.INQ.EQ.1))
          if (LINE(i2:i2).EQ.'"'.AND.i1.eq.i2) INQ=1 ! first quotes
          if (LINE(i2:i2).EQ.'"'.AND.i1.LT.i2) INQ=0 ! second quotes
          write(*,1) IACHAR(LINE(i2:i2))
          i2=i2+1
        END DO
        write(*,*)
        write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
        IF (LINE(i1:i1).EQ.'"') i1=i1+1
        IF ((i2.gt.1).and.(LINE(i2-1:i2-1).EQ.'"')) i2=i2-1
        write(*,*) k,' FINDPAR i1=',i1,' i2=',i2
        write(*,*) k,' FINDPAR ',LINE(i2:i2)
        IF (K.LT.IPAR) then
          if (i2.gt.L) then
            i2=i1
            exit
          else if (LINE(i2:i2).EQ.'"') then
            I1=I2+1 ! skip quotes
          else
            I1=I2 ! try next parameter
          endif
          if (I1.gt.L) then
            i2=i1
            exit
          endif
        endif
      ENDDO
      ILEN=I2-I1
      ISTART=I1
      write(*,*) k,' FINDPAR result=',I1,I2,' [',LINE(ISTART:ISTART+ILEN-1),']'
      END SUBROUTINE FINDPAR_NEW


!------------------------ procedures handling STRING BUFFERS (strings delimited by a character, like 'asdf:ghjk:' ------


C-----------------------------------------------------------------
      SUBROUTINE ALLOCSTRPAR(LINE,DLM,NARG)
C Defines buffer with NARG items, delimitted by DLM
C filled by '0'
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*),intent(out) :: LINE
      CHARACTER*1,intent(in) :: DLM
      integer,intent(in) :: NARG
      integer :: i,LL
        LL=LEN(LINE)
        LINE=' '
        if (LL.le.NARG*2) then
          do i=1,NARG
            call APPENDPAR(LINE,DLM,'0')
          enddo
        else
          write(*,*) 'Error on ALLOCSTRPAR, buffer is too small'
        endif
      end SUBROUTINE ALLOCSTRPAR

C-----------------------------------------------------------------
      SUBROUTINE FINDSTRPAR(LINE,DLM,IPAR,ISTART,ILEN)
C Finds IPAR-th string parameter found on LINE
C DLM ... delimiter character
C returns starting position in ISTART and length in ILEN
C ILEN<0 ... not found
C ILEN=0 ... empty string (e.g. the 2nd item in  "abc::def:qwerty:")
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      CHARACTER*1 DLM
      INTEGER*4 L,i1,i2,IPAR,ISTART,ILEN,K
      
      call BOUNDS(LINE,i1,L)

      ! L=LEN_TRIM(LINE)
      IF (L.LT.1) GOTO 99 ! no characters
      K=1
      !i1=1
      i2=INDEX(LINE(i1:L),DLM)
  !    if (DLM.eq.' ') write(*,*) 'FINDSTRPAR [',LINE(i1:L),']',i1,i2,L,IPAR
      IF (i2.le.0) THEN ! no delimiter
        IF (IPAR.GT.1) GOTO 99
        ISTART=i1
        ILEN=L-i1+1
        RETURN
      ENDIF
      DO WHILE ((i2.GT.0).AND.(K.LT.IPAR).AND.(i1+i2.LE.L))
        i1=i1+i2
        i2=INDEX(LINE(i1:L),DLM)
        K=K+1
      ENDDO
      IF (K.lt.IPAR) GOTO 99 ! no IPAR-th item on the list
      IF (i2.LE.0) THEN
         i2=L+1-i1 ! last string does not end with DLM
      ELSE
         i2=i2-1   ! item ends before next DLM
      ENDIF
      ISTART=i1
      ILEN=i2
      RETURN

99    ILEN=-1
      END SUBROUTINE FINDSTRPAR

C-----------------------------------------------------------------
      SUBROUTINE COUNTPAR(LINE,DLM,NP)
C return number of parameters in LINE delimited by DLM
C a trailing delimiter is not counted
C INPUT:
C    LINE ... list of parameter names delimited by :
C    DLM  ... delimiter character
C OUTPUT:
C   NP    ... number of parameters
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) LINE
      CHARACTER*1 DLM
      INTEGER*4 NP,L,i1,i2,K

      NP=0
      L=LEN_TRIM(LINE)
      IF (L.LE.0) RETURN
      i1=1
      i2=1
      K=1
C how many delimiters
      i2=INDEX(LINE(i1:L),DLM)
      DO WHILE ((i2.gt.0).AND.(i1+i2.le.L))
        K=K+1
        i1=i1+i2
        i2=INDEX(LINE(i1:L),DLM)
      ENDDO
      NP=K
      END SUBROUTINE COUNTPAR

C-----------------------------------------------------------------
      SUBROUTINE GETPARINDX(LINE,DLM,SPAR,ID)
C Finds index of the string SPAR as found on the list in LINE (inverse to FINDSTRPAR)
C is case sensitive!
C NOTE: SPAR=' ' results in ID=0
C INPUT:
C    LINE ... list of parameter names delimited by :
C    DLM  ... delimiter character
C    SPAR ... searched string
C OUTPUT:
C   ID    ... string position on the list (=0 if not found)
C-----------------------------------------------------------------
      IMPLICIT NONE
      integer, intent(out) :: ID
      CHARACTER(1),intent(in) :: DLM
      CHARACTER*(*),intent(in) :: LINE,SPAR
      CHARACTER(1) CH
      INTEGER :: LS,LP,LL
      ID=0
      LS=LEN_TRIM(SPAR)
      if (LS.le.0) return
      LL=LEN_TRIM(LINE)
      CH=' '
      if (LL.le.0) return
      if (LINE(LL:LL).ne.DLM) CH=DLM
      LP=INDEX(DLM//trim(LINE)//CH,DLM//TRIM(SPAR)//DLM)
      if (LP.LE.0) return
      CALL COUNTPAR(LINE(1:LP-1+LS),DLM,ID)

      END SUBROUTINE GETPARINDX

C-----------------------------------------------------------------
      SUBROUTINE APPENDPAR(LINE,DLM,SPAR)
C Appends SPAR to the list of parameters LINE delimited by DLM
C INPUT:
C    LINE ... list of parameter names delimited by :
C    DLM  ... delimiter character
C    SPAR ... searched string
C OUTPUT:
C    LINE  ... input with appended SPAR, delimited by DLM
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) LINE,SPAR
      CHARACTER*1 DLM
      character(1024) AUX
      INTEGER*4 I,J,LS
      LS=LEN(LINE)
      I=LEN_TRIM(LINE)
      J=LEN_TRIM(SPAR)
      IF (LS.GT.I+J)THEN
        IF (I.LE.0) THEN
          LINE=trim(trim(adjustl(SPAR))//DLM)
        ELSE
          AUX=trim(LINE)
          if (DLM.ne.' ') then
            LINE=trim(AUX)//trim(adjustl(SPAR))//DLM
          else
            LINE=trim(AUX)//trim(adjustl(SPAR))
          endif
        ENDIF
      else
        write(*,*) 'APPENDPAR: full buffer [',trim(LINE),'] ',trim(SPAR),LS
      ENDIF
      END SUBROUTINE APPENDPAR

C-----------------------------------------------------------------
      SUBROUTINE PREPENDPAR(LINE,DLM,SPAR)
C Pre-pends SPAR to the list of parameters LINE delimited by DLM
C INPUT:
C    LINE ... list of parameter names delimited by :
C    DLM  ... delimiter character
C    SPAR ... searched string
C OUTPUT:
C    LINE  ... input with appended SPAR, delimited by DLM
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*) LINE,SPAR
      character(1024) AUX
      CHARACTER*1 DLM
      INTEGER*4 I,J,LS
      LS=LEN(LINE)
      I=LEN_TRIM(LINE)
      J=LEN_TRIM(SPAR)
!      write(*,*) 'PREPENDPAR: ',LS,I,J
      IF (LS.GT.I+J)THEN
        IF (I.LE.0) THEN
          LINE=trim(trim(adjustl(SPAR))//DLM)
        ELSE
          if (DLM.ne.' ') then
            AUX=DLM//trim(LINE)
          else
            AUX=trim(LINE)
          endif
          LINE=trim(adjustl(SPAR))//trim(AUX)
        ENDIF
      else
        write(*,*) 'PREPENDPAR: full buffer [',trim(LINE),'] ',LS
      ENDIF
      END SUBROUTINE PREPENDPAR

C-----------------------------------------------------------------
      SUBROUTINE SETSTRPAR(LINE,DLM,INDX,SPAR)
C Replace INDXth item of the LINE buffer with SPAR
C INPUT:
C    LINE ... old buffer (DLM-delimited strings)
C    DLM  ... delimiter character
C    INDX ... item index (1-based)
C    SPAR ... searched string
C OUTPUT:
C    LINE  ... new buffer
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*),intent(inout) :: LINE
      CHARACTER(*),intent(in) :: SPAR
      CHARACTER*1,intent(in) :: DLM
      integer,intent(in) :: INDX
      character(1024) AUX
      INTEGER :: IS,IL,IS1,LL,LS,NP
      if (INDX.le.0) return
      LL=LEN(LINE)
      LS=len_trim(LINE)
      call COUNTPAR(LINE,DLM,NP)
      if (NP.ge.INDX) then
        call FINDSTRPAR(LINE,DLM,INDX,IS,IL)
        if ((LS+len_trim(SPAR)-IL).le.LL) then
          IS1=MAX(IS,1)
          if (IS1.gt.1) then
            AUX=LINE(1:IS1-1)//trim(SPAR)
          else
            AUX=trim(SPAR)
          endif
          LINE=trim(AUX)//trim(LINE(IS1+IL:))
        else
          write(*,*) 'Error on SETSTRPAR, buffer is too small ['//trim(LINE)//'] '//trim(SPAR)
        endif
      else
        write(*,*) 'SETSTRPAR: item ',INDX,' is not allocated in buffer ['//trim(LINE)//']'
      endif
      END SUBROUTINE SETSTRPAR

C-----------------------------------------------------------------
      SUBROUTINE GETSTRPAR(LINE,DLM,INDX,SARG)
C Return IPAR-th string parameter found on LINE
C-----------------------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*),intent(inout) :: LINE
      CHARACTER*1,intent(in) :: DLM
      integer,intent(in) :: INDX
      CHARACTER(*),intent(out) :: SARG
      INTEGER :: IS,IL,LL
      LL = LEN(SARG)
      call FINDSTRPAR(LINE,DLM,INDX,IS,IL)
      if (IL>LL) IL=LL
      if ((IS.ge.1).and.(IL.GE.1)) then
        SARG=LINE(IS:IS+IL-1)
      else
        SARG=' '
      endif
      end SUBROUTINE GETSTRPAR


C-----------------------------------------------------------------
      SUBROUTINE ARRAY2STR(ARG,NARG,STR)
C Converts real number into a string with as short format as possible
C INPUT:
C ARG     ... real*8 array
C NARG    ... number of elements in ARG
C OUTPUT:
C STR   ... output string
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL(KIND(1.D0)) :: ARG(*)
      integer :: NARG
      CHARACTER*(*) STR
      character(32) :: CNUM
      INTEGER*4 LS,LA,LL,i
      LL=LEN(STR)
      STR=' '
      do i=1,NARG
        LS=LEN_TRIM(STR)
        call FLOAT2STR(ARG(i),CNUM)
        LA=LEN_TRIM(CNUM)
        if ((LS+LA+1).le.LL) then
          STR=trim(STR)//' '//trim(CNUM)
        endif
      enddo
      end SUBROUTINE ARRAY2STR


C-----------------------------------------------------------------
      SUBROUTINE IARRAY2STR(ARG,NARG,STR)
C Converts real number into a string with as short format as possible
C INPUT:
C ARG     ... integer array
C NARG    ... number of elements in ARG
C OUTPUT:
C STR   ... output string
C-----------------------------------------------------------------
      IMPLICIT NONE
      integer :: ARG(*)
      integer :: NARG
      CHARACTER*(*) STR
      character(32) :: CNUM
      INTEGER*4 LS,LA,LL,i
      LL=LEN(STR)
      STR=' '
      do i=1,NARG
        LS=LEN_TRIM(STR)
        call INT2STR(ARG(i),CNUM)
        LA=LEN_TRIM(CNUM)
        if ((LS+LA+1).le.LL) then
          STR=trim(STR)//' '//trim(CNUM)
        endif
      enddo
      end SUBROUTINE IARRAY2STR

C-----------------------------------------------------------------
      SUBROUTINE FLOAT2STR(Z,STR)
C Converts real number into a string with as short format as possible
C INPUT:
C Z     ... real*8 number
C OUTPUT:
C STR   ... output string
C-----------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 Z,AZ
      CHARACTER*(*) STR
      CHARACTER*64 S
      INTEGER*4 IS,IL, LSTR

1     format(G12.6)
2     format(I9)
3     format(F10.5)
4     format(F10.7)

      AZ=ABS(Z)
      LSTR = LEN(STR)
    ! integer
      IF ((1.D0*NINT(AZ).EQ.AZ).AND.(AZ.LT.1E9)) THEN
        WRITE(S,2) NINT(Z)
    ! number outside fixed formad
      ELSE IF ((AZ.LT.0.001).OR.(AZ.GE.10000)) THEN
        WRITE(S,1) Z
    ! fixed format
      ELSE IF (AZ.LT.1.0) THEN
        WRITE(S,4) Z
      ELSE
        WRITE(S,3) Z
      ENDIF
      CALL BOUNDS(S,IS,IL)

C// strip unsignificant zeros off
      IF ((INDEX(S,'.').GT.0).AND.(INDEX(S,'E').LE.0)) THEN
        DO WHILE (S(IS+IL-1:IS+IL-1).EQ.'0')
          IL=IL-1
        ENDDO
        IF (S(IS+IL-1:IS+IL-1).EQ.'.') IL=IL-1
      ENDIF
      if (IL>LSTR) IL=LSTR
      STR=S(IS:IS+IL-1)
      END SUBROUTINE FLOAT2STR

!-----------------------------------------------------------------
      SUBROUTINE INT2STR(INUM,STR)
! Converts integer number into a string without trailing spaces etc.
! INPUT:
! INUM    ... integer number
! OUTPUT:
! STR   ... output string
!-----------------------------------------------------------------
      IMPLICIT NONE
      integer,intent(in) :: INUM
      CHARACTER*(*) STR
      CHARACTER*64 S
      INTEGER*4 IS,IL,LSTR

      WRITE(S,*) INUM
      LSTR = LEN(STR)
      CALL BOUNDS(S,IS,IL)
      IL=MIN(LEN(STR),IL)
      if (IL>LSTR) IL=LSTR
      STR=S(IS:IS+IL-1)
      END SUBROUTINE INT2STR

C-----------------------------------------------------------------
      SUBROUTINE LASTSUBSTR(STR,SUBSTR,IPOS)
C Find last occurence of the SUBSTR in STR
C Return position in IPOS (=0 if not found)
C-----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 IPOS
      CHARACTER*(*) STR,SUBSTR
      INTEGER*4 IS,LS,ISS,LSS,LL,I

      CALL BOUNDS(STR,IS,LS)
      LL=IS+LS-1
      CALL BOUNDS(SUBSTR,ISS,LSS)
      I=1
      IPOS=0
      IF (LS.GT.0.AND.LSS.GT.0.AND.I.GT.0.AND.IS.LE.LL) THEN
        I=INDEX(STR(IS:LL),SUBSTR(ISS:ISS+LSS-1))
        IF (I.GT.0) THEN
          IPOS=IS+I-1
          IS=IPOS+LSS
        ENDIF
      ENDIF
      END


!-----------------------------------------------------------------
      SUBROUTINE PARSE_ENTITY(LINE)
! parses some basic special characters to HTML entities
!-----------------------------------------------------------------
      use IO
      IMPLICIT NONE
      integer*4 DIM
      parameter(DIM=3) ! nuber of entities
      character*(*) LINE
      character*256 S
      character*8 HTML(DIM)
      character*1 ASCI(DIM)
      integer*4 LL,LS,LH,I,J,K,IRES
C HTML entities and corresponding ASCII characters to be parsed
C Must start with & !!
      DATA HTML/'&amp;','&lt;','&gt;'/
      DATA ASCI/'&'    ,'<'   ,'>'   /
      LL=LEN(LINE)
      S=LINE
      do J=1,DIM
        IRES=1
        K=1
        do while (IRES.GT.0)
          IRES=0
          LS=LEN_TRIM(S)
          I=INDEX(S(K:LS),ASCI(J))
          LH=LEN_TRIM(HTML(J))
          if (I.GT.0.AND.LS.LE.LL-LH) then
            if (K+I-2.LE.0) then
              LINE=trim(HTML(J))//S(K+I:LS)
            else
              LINE=S(1:K+I-2)//trim(HTML(J))//S(K+I:LS)
            endif
            K=K+I+LH-1
            IRES=IRES+I
            S=LINE
          endif
        enddo
      enddo
      END SUBROUTINE PARSE_ENTITY

!-----------------------------------------------------------------
      SUBROUTINE PARSE_INPUT(LINE)
! replace tabs and null characters with spaces
!-----------------------------------------------------------------
      IMPLICIT NONE
      integer, parameter :: ME=2 ! nuber of entities
      character*(*) LINE
      integer,parameter :: ICH(ME)=(/9,0/)
      character*1 CH(ME)
      integer*4 LL,J,K
    !  DATA ICH/9,0/
      LL=len_trim(LINE)
      do K=1,ME
        CH(K)=CHAR(ICH(K))
      enddo
      do J=1,LL
        DO K=1,ME
          if (LINE(J:J).eq.CH(k)) LINE(J:J)=' '
        enddo
      enddo
      END SUBROUTINE PARSE_INPUT
