!//////////////////////////////////////////////////////////////////////
!////  $Id: xmlinfo.f,v 1.13 2019/08/15 17:24:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.13 $
!////     $Date: 2019/08/15 17:24:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Various XML output, e.g. for communication with GUI
!////  Some output procedures handle also non-XML (console) output (XMLOUT=0)
!////
!//////////////////////////////////////////////////////////////////////
      MODULE XMLINFO
      USE IO
      implicit NONE

      integer :: XMLOUT=0

      contains

!--------------------------------------------------
      SUBROUTINE XML_OPENWSTART(FNAME)
! Send start XML tags to STDOUT before passing data to GUI
! for savig in a file
!--------------------------------------------------
      CHARACTER(*),intent(in) :: FNAME
1     FORMAT('<OPENW file="',a,'">',$)
      call XML_RSXDUMP(smes,' ',1)
      write(smes,1) trim(fname)
      write(smes,*) ' <![CDATA['
      end SUBROUTINE XML_OPENWSTART


!--------------------------------------------------
      SUBROUTINE XML_OPENWEND
! close file for write and return unit number
! Only send closing XML tags to STDOUT. File save is handled by GUI !!
!--------------------------------------------------
1     FORMAT('</OPENW>')
      write(smes,*) ']]>'
      write(smes,1)
      call XML_RSXDUMP(smes,' ',0)
      end SUBROUTINE XML_OPENWEND

C--------------------------------------------------------------------
      SUBROUTINE MSG_WARN(MSG,IRSX)
C Print error message to SMES
C INPUT:
C   MSG         ... message text
C   IRSX        ... if > 0, print also the RSXDUMP tags
C--------------------------------------------------------------
      integer :: IRSX
      CHARACTER(LEN=*) :: MSG
      CHARACTER(LEN=256) :: S

10    FORMAT('<WARNING priority="',a,'">',$)
11    FORMAT(a,$)
14    FORMAT('</WARNING>')
15    FORMAT('Warning: ',a)
      IF (XMLOUT.GT.0) THEN
        S=MSG
        call PARSE_ENTITY(S)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',1)
        WRITE(SMES,10) 'low'
        WRITE(SMES,11) TRIM(S)
        WRITE(SMES,14)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',0)
      ELSE
        WRITE(SMES,15) TRIM(MSG)
      ENDIF
      END SUBROUTINE MSG_WARN


C--------------------------------------------------------------------
      SUBROUTINE MSG_ERROR(SRC,MSG,IPR,IRSX)
C Print error message
C INPUT:
C   MSG         ... message text
C   IRSX        ... if > 0, print also the RSXDUMP tags
C   IPR         ... priority
C--------------------------------------------------------------
      integer :: IRSX,IPR
      CHARACTER(LEN=*) :: SRC,MSG
      CHARACTER(LEN=256) :: S

10    FORMAT('<ERROR src="',a,'" priority="',a,'">',$)
11    FORMAT(a,$)
14    FORMAT('</ERROR>')
15    FORMAT('Error in ',a,': ',a)

      IF (XMLOUT.GT.0) THEN
        S=MSG
        call PARSE_ENTITY(S)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',1)
        if (IPR.EQ.0) then
          WRITE(SMES,10) TRIM(SRC),'low'
        else
          WRITE(SMES,10) TRIM(SRC),'high'
        endif
        WRITE(SMES,11) TRIM(S)
        WRITE(SMES,14)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',0)
      ELSE
        WRITE(SMES,15) TRIM(SRC),TRIM(MSG)
      ENDIF
      END SUBROUTINE MSG_ERROR

C--------------------------------------------------------------------
      SUBROUTINE MSG_INFO(MSG,IRSX)
C Print info message
C INPUT:
C   MSG         ... message text
C   IRSX        ... if > 0, print also the RSXDUMP tags
C--------------------------------------------------------------
      integer :: IRSX
      CHARACTER(LEN=*) :: MSG
      CHARACTER(LEN=256) :: S

10    FORMAT('<INFO priority="',a,'">',$)
11    FORMAT(a,$)
14    FORMAT('</INFO>')
15    FORMAT('Info: ',a)

      IF (XMLOUT.GT.0) THEN
        S=MSG
        call PARSE_ENTITY(S)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',1)
        WRITE(SMES,10) 'low'
        WRITE(SMES,11) TRIM(S)
        WRITE(SMES,14)
        if (IRSX.GT.0) call XML_RSXDUMP(SMES,' ',0)
      ELSE
        WRITE(SMES,15) TRIM(MSG)
      ENDIF
      END SUBROUTINE MSG_INFO

C--------------------------------------------------------------------
      SUBROUTINE XML_RSXDUMP(IU,VALNAME,START)
C Main document tag for RESTRAX output
C format is:
C <RSXDUMP>...</RSXDUMP> for VALNAME=' '
C <RSXDUMP><VALNAME>...</VALNAME></RSXDUMP> for VALNAME=<>' '
C <RSXDUMP><VALNAME raytracing="false|true">...</VALNAME></RSXDUMP>
C     for VALNAME in 'BRAG|FWHM|RES|FIT|PHON|'
C INPUT:
C   IU          ... output unit
C   VALNAME     ... command name
C   START       ... START(1) or END(0) tag
C--------------------------------------------------------------
      integer :: IU,START
      CHARACTER(LEN=*) :: VALNAME
      CHARACTER(LEN=5) :: BOOL(0:1)
      integer :: LL
c commands with "raytracing" attribute
      DATA BOOL/'false','true'/

7     FORMAT('<RSXDUMP>')
8     FORMAT('</RSXDUMP>')
c 9     FORMAT('<RSXDUMP>',/,'<',a,' raytracing="',a,'">')
10    FORMAT('<RSXDUMP>',/,'<',a,'>')
11    FORMAT('</',a,'>',/,'</RSXDUMP>')
19    FORMAT('M',a,', Monte Carlo (NESS):')
c20    FORMAT(a,', Analytical (TRAX):')
c21    FORMAT(72('-'))
21    format(a)

      LL=LEN_TRIM(VALNAME)
c only RSXDUMP tag
      if (LL.le.0) then
        IF (XMLOUT.GT.0) THEN
          if (start.gt.0) then
            WRITE(IU,7)
          else
            WRITE(IU,8)
          endif
        else if (start.le.0) then
          WRITE(IU,21) ' '
        endif
      else
c XML output
        IF (XMLOUT.GT.0) THEN
          if (start.gt.0) then
            WRITE(IU,10) TRIM(VALNAME)
          else
            WRITE(IU,11) TRIM(VALNAME)
          endif
        else
c normal output
          if (start.gt.0) then
            WRITE(IU,21) TRIM(VALNAME)
          else
            WRITE(IU,21) ' '
          endif
        endif
      endif
      END SUBROUTINE XML_RSXDUMP

C--------------------------------------------------------------------
      SUBROUTINE XML_ARRAY(IU,NCOL,CAPTION,HCOL,AMAT)
C Print XML output tag with a matrix of float values
C as XML_MATRIX, but takes only 1-dim array
C INPUT:
C   IU          ... output unit
C   NCOL        ... number of columns
C   CAPTION     ... table caption
C   HCOL        ... column headers (delimited by |)
C   AMAT        ... array with matrix values
C--------------------------------------------------------------
      INTEGER, INTENT(in) :: IU,NCOL
      CHARACTER(LEN=*) , intent(in) :: CAPTION,HCOL
      REAL(KIND(1.D0)), INTENT(in) :: AMAT(NCOL)
      CHARACTER, PARAMETER :: TB=CHAR(9)
      integer i,j,is,ilen
      CHARACTER(LEN=16) :: CNUM
      CHARACTER(LEN=256) :: S,SC,SR

1     FORMAT('<MATRIX ncol="',I3,'" nrow="1">')
2     FORMAT('<CAPTION>',a,'</CAPTION>')
3     FORMAT('<TH>',a,'</TH>')
4     FORMAT('<TD>',G12.4,'</TD>')
10    FORMAT('</MATRIX>')
11    format(a)
12    format(G12.4)

c  XML OUTPUT
      IF (XMLOUT.GT.0) THEN
        S=CAPTION
        SC=HCOL
        SR=' '
        call PARSE_ENTITY(S)
        call PARSE_ENTITY(SC)
        call PARSE_ENTITY(SR)
        write(IU,1) NCOL
c caption
        write(IU,2) TRIM(S)
c column headers
        if (HCOL.NE.' ') THEN
          DO I=1,NCOL
            call FINDSTRPAR(SC,'|',I,IS,ILEN)
            IF (ILEN.GT.0) THEN
              write(IU,3) SC(IS:IS+ILEN-1)
            ELSE
              write(IU,3) ' '
            ENDIF
          ENDDO
        endif
        write(IU,11) '<TR>'
c row header
        call FINDSTRPAR(SR,'|',I,IS,ILEN)
        IF (ILEN.GT.0) THEN
          write(IU,3) SR(IS:IS+ILEN-1)
        ELSE
          write(IU,3) ' '
        ENDIF
c row values
        DO J=1,NCOL
          write(IU,4) AMAT(J)
        ENDDO
        write(IU,11) '</TR>'
        write(IU,10)

c  NORMAL OUTPUT
      ELSE

101     format(a14,$)
102     format(a)
c caption
        write(IU,11) TRIM(caption)//':'
c column headers
        DO I=1,NCOL
          CNUM=' '
          call FINDSTRPAR(HCOL,'|',I,IS,ILEN)
          IF (ILEN.GT.0) CNUM=HCOL(IS:IS+ILEN-1)
          write(IU,101) CNUM
        ENDDO
        write(IU,*)
c row values
        DO J=1,NCOL
          write(CNUM,12) AMAT(j)
          write(IU,101) CNUM
        ENDDO
        write(IU,*)
      ENDIF
      END SUBROUTINE XML_ARRAY

C--------------------------------------------------------------------
      SUBROUTINE XML_IARRAY(IU,NCOL,CAPTION,HCOL,AMAT)
C Print XML output tag with a matrix of int values
C as XML_MATRIX, but takes only 1-dim array
C INPUT:
C   IU          ... output unit
C   NCOL        ... number of columns
C   CAPTION     ... table caption
C   HCOL        ... column headers (delimited by |)
C   AMAT        ... array with matrix values
C--------------------------------------------------------------
      INTEGER, INTENT(in) :: IU,NCOL
      CHARACTER(LEN=*) , intent(in) :: CAPTION,HCOL
      integer, INTENT(in) :: AMAT(NCOL)
      CHARACTER, PARAMETER :: TB=CHAR(9)
      integer i,j,is,ilen
      CHARACTER(LEN=16) :: CNUM
      CHARACTER(LEN=256) :: S,SC,SR

1     FORMAT('<MATRIX ncol="',a,'" nrow="1">')
2     FORMAT('<CAPTION>',a,'</CAPTION>')
3     FORMAT('<TH>',a,'</TH>')
4     FORMAT('<TD>',a,'</TD>')
10    FORMAT('</MATRIX>')
11    format(a)
12    format(I8)

c  XML OUTPUT
      IF (XMLOUT.GT.0) THEN
        S=CAPTION
        SC=HCOL
        SR=' '
        call PARSE_ENTITY(S)
        call PARSE_ENTITY(SC)
        call PARSE_ENTITY(SR)
        call INT2STR(NCOL,CNUM)
        write(IU,1) trim(CNUM)
c caption
        write(IU,2) TRIM(S)
c column headers
        if (HCOL.NE.' ') THEN
          DO I=1,NCOL
            call FINDSTRPAR(SC,'|',I,IS,ILEN)
            IF (ILEN.GT.0) THEN
              write(IU,3) SC(IS:IS+ILEN-1)
            ELSE
              write(IU,3) ' '
            ENDIF
          ENDDO
        endif
        write(IU,11) '<TR>'
c row header
        call FINDSTRPAR(SR,'|',I,IS,ILEN)
        IF (ILEN.GT.0) THEN
          write(IU,3) SR(IS:IS+ILEN-1)
        ELSE
          write(IU,3) ' '
        ENDIF
c row values
        DO J=1,NCOL
          call INT2STR(AMAT(J),CNUM)
        !  write(CNUM,*) AMAT(J)
          write(IU,4) trim(CNUM)
        ENDDO
        write(IU,11) '</TR>'
        write(IU,10)

c  NORMAL OUTPUT
      ELSE

101     format(a9,$)
102     format(a)
c caption
        write(IU,11) TRIM(caption)//':'
c column headers
        DO I=1,NCOL
          CNUM=' '
          call FINDSTRPAR(HCOL,'|',I,IS,ILEN)
          IF (ILEN.GT.0) CNUM=HCOL(IS:IS+ILEN-1)
          write(IU,101) trim(CNUM)
        ENDDO
        write(IU,*)
c row values
        DO J=1,NCOL
          write(CNUM,*) AMAT(j)
          write(IU,101) trim(CNUM)
        ENDDO
        write(IU,*)
      ENDIF

      END SUBROUTINE XML_IARRAY

C--------------------------------------------------------------------
      SUBROUTINE XML_MATRIX(IU,NCOL,NROW,CAPTION,HCOL,HROW,AMAT,NMAT)
C Print XML output tag with a matrix of float values
C INPUT:
C   IU          ... output unit
C   NCOL        ... number of columns
C   NROW        ... number of rows
C   CAPTION     ... table caption
C   HCOL        ... column headers (delimited by |)
C   HROW        ... row headers (delimited by |)
C   AMAT        ... array with matrix values
C   NMAT        ... leading (column) dimension of AMAT
C--------------------------------------------------------------
      INTEGER, INTENT(in) :: IU,NCOL,NROW,NMAT
      CHARACTER(LEN=*) , intent(in) :: CAPTION,HCOL,HROW
      REAL(KIND(1.D0)), INTENT(in) :: AMAT(NMAT,NCOL)
      CHARACTER, PARAMETER :: TB=CHAR(9)
      integer i,j,is,ilen
      CHARACTER(LEN=16) :: CNUM
      CHARACTER(LEN=256) :: S,SC
      CHARACTER(LEN=2048) :: SR

1     FORMAT('<MATRIX ncol="',I3,'" nrow="',I3,'">')
2     FORMAT('<CAPTION>',a,'</CAPTION>')
3     FORMAT('<TH>',a,'</TH>')
4     FORMAT('<TD>',G12.4,'</TD>')
10    FORMAT('</MATRIX>')
11    format(a)
12    format(G12.4)

c      write(*,*) 'XML_MATRIX [',trim(CAPTION),']'
c  XML OUTPUT
      IF (XMLOUT.GT.0) THEN
        S=CAPTION
        SC=HCOL
        SR=HROW
        call PARSE_ENTITY(S)
        call PARSE_ENTITY(SC)
        call PARSE_ENTITY(SR)
        write(IU,1) NCOL,NROW
c caption
        write(IU,2) TRIM(S)
c column headers
        if (HCOL.NE.' ') THEN
          DO I=1,NCOL
            call FINDSTRPAR(SC,'|',I,IS,ILEN)
            IF (ILEN.GT.0) THEN
              write(IU,3) SC(IS:IS+ILEN-1)
            ELSE
              write(IU,3) ' '
            ENDIF
          ENDDO
        endif
c rows
        DO I=1,NROW
          write(IU,11) '<TR>'
c row header
          call FINDSTRPAR(SR,'|',I,IS,ILEN)
          IF (ILEN.GT.0) THEN
            write(IU,3) SR(IS:IS+ILEN-1)
          ELSE
            write(IU,3) ' '
          ENDIF
c row values
          DO J=1,NCOL
            write(IU,4) AMAT(I,J)
          ENDDO
          write(IU,11) '</TR>'
        ENDDO
        write(IU,10)

c  NORMAL OUTPUT
      ELSE

101     format(a14,$)
102     format(a)
c caption
        write(IU,11) TRIM(caption)//':'
c column headers
        DO I=1,NCOL
          CNUM=' '
          call FINDSTRPAR(HCOL,'|',I,IS,ILEN)
          IF (ILEN.GT.0) CNUM=HCOL(IS:IS+ILEN-1)
          write(IU,101) CNUM
        ENDDO
        write(IU,*)
c rows
        DO I=1,NROW
c row values
          DO J=1,NCOL
            write(CNUM,12) AMAT(I,J)
            write(IU,101) CNUM
          ENDDO
c row header
          CNUM=' '
          call FINDSTRPAR(HROW,'|',I,IS,ILEN)
          IF (ILEN.GT.0) write(IU,102) HROW(IS:IS+ILEN-1)
        ENDDO
        write(IU,*)
      ENDIF

      END SUBROUTINE XML_MATRIX

C----------------------------------------------------------------------
      SUBROUTINE XML_VALUE(IU,VALNAME,VALTYP,VALUNI,INTVAL,FLOATVAL,STRVAL)
C Print XML output tag with signgle value data on one row
C format <VALUE name="" units="" >value</VALUE>
C INPUT:
C   IU          ... output unit
C   VALNAME     ... data name
C   VALTYP      ... data type (float,int,string)
C   VALUNI      ... data units (any string, e.g. meV, deg, ...)
C   INTVAL      ... integer value
C   FLOTVAL     ... float value
C   STRVAL      ... string value
C----------------------------------------------------------------------
      integer :: IU
      CHARACTER(LEN=*) :: VALNAME,VALTYP,VALUNI,STRVAL
      real(KIND(1.D0)) :: FLOATVAL
      integer :: INTVAL
      CHARACTER(LEN=6) :: VT
      CHARACTER(LEN=256) :: CHVAL
      integer :: LL
      CHARACTER(LEN=6) :: UPCASE

1     FORMAT(a1,$)
8     FORMAT(a,'=',a,' [',a,']')
9     FORMAT(a,'=',a)
c10    FORMAT('<VALUE name="',a,'" type="',a,'" units="',a'">',$)
c no type attribute
10    FORMAT('<VALUE name="',a,'" units="',a,'">',$)
11    FORMAT(G12.5)
12    FORMAT(G12.5)
13    FORMAT(a)
14    FORMAT(a,$)
15    FORMAT('</VALUE>')

C format value string
      LL=MIN(5,LEN(VALTYP))
      VT=UPCASE(VALTYP(1:LL))
      CHVAL=' '
      IF (VT(1:3).EQ.'INT') THEN
        write(CHVAL,11) INTVAL
      ELSE IF (VT(1:5).EQ.'FLOAT') THEN
        write(CHVAL,12) FLOATVAL
      ELSE IF (VT(1:6).EQ.'STRING') THEN
        LL=MIN(256,LEN(STRVAL))
        write(CHVAL,13) TRIM(STRVAL(1:LL))
      ENDIF

      IF (XMLOUT.GT.0) THEN
        call PARSE_ENTITY(CHVAL)
        WRITE(IU,10) TRIM(VALNAME),TRIM(VALUNI)
        WRITE(IU,14) TRIM(CHVAL)
        write(IU,15)
      ELSE
        IF (VALUNI.EQ.' ') THEN
          write(IU,9) TRIM(VALNAME),TRIM(CHVAL)
        ELSE
          write(IU,8) TRIM(VALNAME),TRIM(CHVAL),TRIM(VALUNI)
        ENDIF
      ENDIF
      END SUBROUTINE XML_VALUE

C----------------------------------------------------------------------
      SUBROUTINE XML_FVALUE(IU,VALNAME,VALUNI,FLOATVAL,ERRVAL)
C Print XML output tag: floating value with error
C <FVALUE name="" units="" >value</FVALUE>
C or if ERRVAL<>0
C <FVALUE name="" units="" >
C    <ERROR>error</ERROR>
C    value
C </FVALUE>
C INPUT:
C   IU          ... output unit
C   VALNAME     ... data name
C   VALUNI      ... data units (any string, e.g. meV, deg, ...)
C   FLOTVAL     ... float value
C   ERRVAL      ... error value
C----------------------------------------------------------------------
      INTEGER :: IU
      CHARACTER(*) :: VALNAME,VALUNI
      REAL(KIND(1.D0)) :: FLOATVAL,ERRVAL
      CHARACTER(256) :: CHVAL,CHERR

1     FORMAT(G12.5)

5     FORMAT(a,'=',a,$)
6     FORMAT(' +- ',a,$)
7     FORMAT(' [',a,']')

10    FORMAT('<FVALUE name="',a,'" units="',a,'">',$)
11    FORMAT(a,$)
12    FORMAT('<ERROR>',a,'</ERROR>',$)
13    FORMAT('</FVALUE>')

C format value string
      CHVAL=' '
      CHERR=' '
      write(CHVAL,1) FLOATVAL
      IF (ERRVAL.GT.0.D0) write(CHERR,1) ERRVAL

      IF (XMLOUT.GT.0) THEN
        call PARSE_ENTITY(CHVAL)
        call PARSE_ENTITY(CHERR)
        WRITE(IU,10) TRIM(VALNAME),TRIM(VALUNI)
        WRITE(IU,11) TRIM(CHVAL)
        IF (ERRVAL.GT.0.D0) WRITE(IU,12) TRIM(CHERR)
        write(IU,13)
      ELSE
        write(IU,5) TRIM(VALNAME),TRIM(CHVAL)
        if (ERRVAL.GT.0.D0) write(IU,6) TRIM(CHERR)
        IF (VALUNI.NE.' ') then
          write(IU,7) TRIM(VALUNI)
        else
          write(IU,*)
        endif
      ENDIF
      END SUBROUTINE XML_FVALUE

C----------------------------------------------------------------------
      SUBROUTINE XML_TITLE(IU,XMLDATA)
C Print XML output tag: header line
C <TITLE>xmldata</TITLE>
C----------------------------------------------------------------------
      INTEGER :: IU
      CHARACTER(*) :: XMLDATA
11    FORMAT(a)
12    FORMAT('<TITLE>',/,a,/,'</TITLE>')

      IF (XMLOUT.GT.0) THEN
        call PARSE_ENTITY(XMLDATA)
        WRITE(IU,12) TRIM(XMLDATA)
      ELSE
        WRITE(IU,11) TRIM(XMLDATA)
      ENDIF
      END SUBROUTINE XML_TITLE

C--------------------------------------------------------------------
      SUBROUTINE XML_TAG(IU,VALNAME,VALSTR,START)
C Print start or end tag
c format for start=2:
c  <VALNAME>valstr</VALNAME> for start=2
c format for start=1:
c   <VALNAME>
c   valstr
c format for start=0:
c   </VALNAME>
C INPUT:
C   IU          ... output unit
C   VALNAME     ... command name
C   VALSTR      ... contents
C   START       ... START(1) or END(0) tag  or both (2)
C--------------------------------------------------------------
      integer :: IU,START,LL
      CHARACTER(LEN=*) :: VALNAME,VALSTR
      CHARACTER(LEN=256) :: S

1     FORMAT(a)
8     format('<',a,' />')
9     format('<',a,'>',a,'</',a,'>')
10    FORMAT('<',a,'>')
11    FORMAT('</',a,'>')
19    FORMAT(a)
20    FORMAT(a,'=',a)

      LL=LEN_TRIM(VALSTR)
      IF (XMLOUT.GT.0) THEN
        S=VALSTR
        call PARSE_ENTITY(S)
        select case(start)
        case(2)
          if (VALSTR.EQ.' ') then
            WRITE(IU,8) TRIM(VALNAME)
          else
            WRITE(IU,9) TRIM(VALNAME),TRIM(S),TRIM(VALNAME)
          endif
        case(1)
          WRITE(IU,10) TRIM(VALNAME)
          if (LL.GT.0) WRITE(IU,1) TRIM(S)
        case(0)
          WRITE(IU,11) TRIM(VALNAME)
        end select
      else if ((start.gt.0).AND.(LL.GT.0)) then
          WRITE(IU,20) TRIM(VALNAME),TRIM(VALSTR)
      endif
      END SUBROUTINE XML_TAG

C--------------------------------------------------------------------
      SUBROUTINE XML_SWITCH(IU,VALNAME,PAR)
C Print a switch tag
c format <VALNAME state="on|off"/>
C INPUT:
C   IU          ... output unit
C   VALNAME     ... switch name
C   PAR         ... (0|1)
C--------------------------------------------------------------
      integer :: IU,PAR
      CHARACTER(LEN=*) :: VALNAME

10    FORMAT('<',a,' state="on" />')
11    FORMAT('<',a,' state="off" />')
20    FORMAT(a,' switched on')
21    FORMAT(a,' switched on')

      IF (XMLOUT.GT.0) THEN
        call XML_RSXDUMP(IU,' ',1)
        if (par.gt.0) then
          WRITE(IU,10) TRIM(VALNAME)
        else
          WRITE(IU,11) TRIM(VALNAME)
        endif
        call XML_RSXDUMP(IU,' ',0)
      else
        if (par.gt.0) then
          WRITE(IU,20) TRIM(VALNAME)
        else
          WRITE(IU,21) TRIM(VALNAME)
        endif
      endif
      END SUBROUTINE XML_SWITCH

C--------------------------------------------------------------------
      SUBROUTINE XML_TASPAR(IU,VALNAME,VALTYP,VALUNI,INTVAL,FLOATVAL,STRVAL)
C Print XML output a TAS parameter
C format <VALNAME units="" >value</VALNAME>
C INPUT:
C   IU          ... output unit
C   VALNAME     ... data name
C   VALTYP      ... data type (float,int,string)
C   VALUNI      ... data units (any string, e.g. meV, deg, ...)
C   INTVAL      ... integer value
C   FLOTVAL     ... float value
C   STRVAL      ... string value
C--------------------------------------------------------------
      integer :: IU
      CHARACTER(LEN=*) :: VALNAME,VALTYP,VALUNI,STRVAL
      real(KIND(1.D0)) :: FLOATVAL
      integer :: INTVAL
      CHARACTER(LEN=6) :: VT
      CHARACTER(LEN=256) :: CHVAL
      integer :: LL
      CHARACTER(LEN=6) :: UPCASE

1     FORMAT(a1,$)
8     FORMAT(a,'=',a,' [',a,']')
9     FORMAT(a,'=',a)
10    FORMAT('<',a,' units="',a,'">',$)
11    FORMAT(G12.5)
12    FORMAT(G12.5)
13    FORMAT(a)
14    FORMAT(a,$)
15    FORMAT('</',a,'>')

C format value string
      LL=MIN(5,LEN(VALTYP))
      VT=UPCASE(VALTYP(1:LL))
      CHVAL=' '
      IF (VT(1:3).EQ.'INT') THEN
        write(CHVAL,11) INTVAL
      ELSE IF (VT(1:5).EQ.'FLOAT') THEN
        write(CHVAL,12) FLOATVAL
      ELSE IF (VT(1:6).EQ.'STRING') THEN
        LL=MIN(256,LEN(STRVAL))
        write(CHVAL,13) TRIM(STRVAL(1:LL))
      ENDIF

      IF (XMLOUT.GT.0) THEN
        WRITE(IU,10) TRIM(VALNAME),TRIM(VALUNI)
        WRITE(IU,14) TRIM(CHVAL)
        write(IU,15) TRIM(VALNAME)
      ELSE
c units are shown only in the XML output
        write(IU,9) TRIM(VALNAME),TRIM(CHVAL)
      ENDIF
      END SUBROUTINE XML_TASPAR

C--------------------------------------------------------------------


C--------------------------------------------------------------------
      SUBROUTINE XML_PROGRESS(IU,STATE,NMAX,NPOS,CAPTION)
C Print XML tag: information about MC simuleation progress
C <PROGRESS><START nmax="100" />Caption</PROGRESS> .. okno se otevre modalne
C <PROGRESS><STEP  n="5" /></PROGRESS> .. nastavi se pozice
C <PROGRESS><STOP /></PROGRESS> .. okno se zavre

C INPUT:
C   IU          ... output unit
C   STATE   ... (0) start (1) step (2) end
C   TESTIM  ... estimeted time
C   TPASSED ... passed time
C   NCMAX   ... requested events
C   NC      ... passed events
C   NEVENT  ... total attempts
C------------------------------------------------------------------
      integer :: IU,STATE,NMAX,NPOS
      CHARACTER(LEN=*) :: CAPTION
      CHARACTER(LEN=10) :: CNC,CNCMAX
10    FORMAT('<START nmax="',a,'">',a,'</START>')
11    FORMAT('<STEP n="',a,'" />')
12    format('<STOP />')
1     FORMAT('.',$)
3     FORMAT(' Wait please, ',a,' in progress',$)
4     FORMAT(' finished')

c XML output
      if (XMLOUT.GT.0) then
c format numbers
        call FLOAT2STR(1.D0*NPOS,CNC)
        call FLOAT2STR(1.D0*NMAX,CNCMAX)
c initial tag
        call XML_RSXDUMP(IU,'PROGRESS',1)
        select case (STATE)
        case(0)
          write(IU,10) trim(CNCMAX),trim(CAPTION)
        case(1)
          write(IU,11) trim(CNC)
        case(2)
          write(IU,12)
        end select
        call XML_RSXDUMP(IU,'PROGRESS',0)
      else
        select case (STATE)
        case(0)
          write(IU,3) trim(CAPTION)
        case(1)
          write(IU,1)
        case(2)
          write(IU,4)
        end select
      endif
      END SUBROUTINE XML_PROGRESS

!--------------------------------------------------------------------
      SUBROUTINE XML_STATUS(IU, STATE, KEY, MSG)
! Print XML tag: information about programe state
! <STATUS value=STATE>message</STATUS> .. okno se otevre modalne
C------------------------------------------------------------------
      INTEGER :: IU
      CHARACTER(LEN=*) :: STATE, KEY, MSG
10    FORMAT('<STATUS value="',a,'" key="',a,'">',a,'</STATUS>')
11    FORMAT('<STATUS value="',a,'">',a,'</STATUS>')
      if (XMLOUT.GT.0) then
        call XML_RSXDUMP(IU,' ',1)
        if (LEN_TRIM(KEY).eq.0) then
          write(IU,11) trim(STATE), trim(MSG)
        else
          write(IU,10) trim(STATE), trim(KEY), trim(MSG)
        endif

        call XML_RSXDUMP(IU,' ',0)
      endif
      END SUBROUTINE XML_STATUS



      end MODULE XMLINFO
