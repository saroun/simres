C///////////////////////////////////////////////////////////////////////////
C////                                                                   ////
C////  LINP - v.1.2, (c) J.Saroun, 1999-2001                            ////
C////                                                                   ////
C////  Universal command line interpreter                               ////
C////                                                                   ////
C////                                                                   ////
C///////////////////////////////////////////////////////////////////////////
C////
C////
C////  Usage:
C////
C////  CALL LINPSET(NLINES,' PROMPT',COMMANDS,HINTS)
C////   ..... to set the prompt, commands and hints
C////
C////  LINE=LINPEXEC(ICOM,NPAR)
C////   ..... reads the input and returns command number(ICOM),
C/////        number of following parameters (NPAR)
C////         and the rest of the input string after the command (LINE)
C////
C////  LINE=FUNCTION LINPEXECSTR(SCOMM,ICOM,NPAR)
C////  .....  as LINPEXEC, but treats the string SCOMM instead of std. input
C////
C////  CALL LINPGETIO(IN,OUT,ERR), LINPSETIO(IN,OUT,ERR)
C////   .... get or set the inout, output and error unit numbers
C////
C////
C///////////////////////////////////////////////////////////////////////////

C     ---------------------------------------------------
      INTEGER*4 FUNCTION ORDCOM(what,commands,ncmd)
C     returns ordinal number of command
C     copy of GETICOM from "linp", but with var 'commands' argument
C     ---------------------------------------------------
      IMPLICIT NONE
      CHARACTER*(*) what
      INTEGER*4 i,L,ipos,nc,ncmd
      CHARACTER*(*) commands(ncmd)

      nc=0
      ORDCOM=0
      ipos=1
      CALL FINDPAR(what,1,ipos,L)
      if (L.LE.0) RETURN
      CALL MKUPCASE(what(ipos:ipos+L-1))
      DO i=1,ncmd
        IF (INDEX(commands(i),what(ipos:ipos+L-1)).eq.1) THEN
          ORDCOM=I
          nc=nc+1
        ENDIF
        IF (nc.gt.1) THEN
           ORDCOM=-3    ! ambiguous command
           RETURN
        ENDIF
      ENDDO

      END


C---------------------------------------------------
      INTEGER*4 FUNCTION GETICOM(what)
C returns ordinal number of command
C command not found  ... return 0
C ambiguous command  ... return -3
C---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      CHARACTER*(*) what
      INTEGER*4 i,L,ipos,nc,i1,i2
      CHARACTER*5 cm

      nc=0
      GETICOM=0
      ipos=1
      CALL FINDPAR(what,1,ipos,L)
      IF (L.GT.5) L=5

      if (L.LE.0) RETURN
      CALL MKUPCASE(what(ipos:ipos+L-1))
      DO i=1,linp_nc
        cm=linp_c(i)
        CALL MKUPCASE(cm)
        IF (INDEX(cm,what(ipos:ipos+L-1)).eq.1) THEN
          GETICOM=I
          nc=nc+1
          i1=1
          CALL FINDPAR(cm,1,i1,i2)
          IF (cm(i1:i1+i2-1).EQ.what(ipos:ipos+L-1)) RETURN
        ENDIF
      ENDDO
      IF (nc.gt.1) THEN
           GETICOM=-3    ! ambiguous command
           RETURN
      ENDIF

      END

C     ---------------------------------------------------
      SUBROUTINE LINPGETIO(IN,OUT,ERR)
C     ---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 IN,OUT,ERR
        IN=linp_in
        OUT=linp_out
        ERR=linp_err
      END


C     ---------------------------------------------------
      SUBROUTINE LINPSETIO(IN,OUT,ERR)
C     ---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 IN,OUT,ERR
        linp_in=IN
        linp_out=OUT
        linp_err=ERR
        linp_eof=0
      END

C     ---------------------------------------------------
      SUBROUTINE LINPSET(nc,prompt,commands,hints)
C     ---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 nc,I,dimc,dimh
      CHARACTER*(*) commands(nc)
      CHARACTER*(*) hints(nc),prompt

      CALL BOUNDS(prompt,I,linp_np)
      IF (I+linp_np-1.GT.20) linp_np=21-I
      linp_p=prompt(1:I+linp_np-1)
C if nc=0, only prompt can be set up
      if (nc.le.0) RETURN
      dimc=LEN(commands(1))
      dimh=LEN(hints(1))
      linp_nc=nc
      IF (dimc.GT.5) dimc=5
      IF (dimh.GT.60) dimh=60
      IF (linp_nc.gt.linp_dim) linp_nc=linp_dim
      DO I=1,linp_nc
        linp_c(I)=commands(I)(1:dimc)
        linp_h(I)=hints(I)(1:dimh)
c        write(*,*) i,' ',linp_c(i),'  ',linp_h(I)
      ENDDO

      END

C---------------------------------------------------
      CHARACTER*(*) FUNCTION LINPEXEC(ICOM,NPAR)
c read a string from input and treat through LINPEXECSTR
c output:
c    ICOM: command ID
c    NPAR: number of command arguments
c return: command with arguments
c NOTE:
c if end of input file => ICOM=-4
C---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 ICOM,NPAR,LL
      CHARACTER*1024 LINE,LINPEXECSTR
      CHARACTER*1 CH

1     FORMAT(a)
2     FORMAT(a,$)
10    IF (linp_in.EQ.5) WRITE(linp_out,2) linp_p(1:linp_np)//'> '

      IF (linp_eof.GT.0) GOTO 20
      READ(linp_in,1,END=20) LINE
      CH=LINE(1:1)
      IF (CH.EQ.'#'.OR.CH.EQ.' '.OR.CH.EQ.char(0)) GOTO 10
      IF (LINE(1:1).EQ.'#') GOTO 10
      LL=LEN(LINE)
c      IF (linp_in.NE.5) WRITE(linp_out,1) linp_p(1:linp_np)//
c     &   '> '//LINE(1:LL)
      LINPEXEC=LINPEXECSTR(LINE(1:LL),ICOM,NPAR)
      IF (ICOM.EQ.-3) GOTO 10
      RETURN

20    LINPEXEC='EOF'
c      write(*,*) 'linpexec: EOF'
      NPAR=0
      ICOM=-4
      END

c---------------------------------------------------
      CHARACTER*(*) FUNCTION LINPEXECSTR(SCOMM,ICOM,NPAR)
c  Treat command string
c  input:
c    SCOMM ... command string
c  output:
c    ICOM: command ID
c    NPAR: number of command arguments
c  return: command with arguments
c  NOTE:
c  - ICOM=-5 ... command is an integer => return this integer as NPAR
c  - ICOM=-4 ... end of input file
c  - ICOM=-3 ... ambiguous command
c  - ? gives a list of commands with hints
c  -
c---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 ICOM,I,J,K,L,IP,LL,GETICOM,NPAR
      CHARACTER*1024 LINE
      CHARACTER*(*) SCOMM
      LOGICAL*4 IsInteger

      DATA linp_in,linp_out,linp_err /5,6,7/
      DATA linp_np,linp_p /4,'LINP'/

1     FORMAT(a)
c3     FORMAT(I)
      L=LEN_TRIM(SCOMM)
      LINE=SCOMM(1:L)

C      READ(linp_in,1) LINE

c// try first, whether the input is integer
C// if yes, return it as string, and in NPAR as value
      IF (IsInteger(LINE(1:L),I)) THEN
        LINPEXECSTR=SCOMM(1:L)
        ICOM=-5
        NPAR=I
        RETURN
      ENDIF

20    ICOM=GETICOM(line)     ! find command number
      LL=LEN_TRIM(line)
      IF (ICOM.EQ.-3) WRITE(linp_err,1) 'Ambiguous command !'
      IF (ICOM.LT.0) RETURN

      K=1
  !    write(*,*) 'linp 1 [',trim(line),']'
      CALL FINDPAR(line,1,K,L) ! command starts at K
  !    write(*,*) 'linp 2 ',ICOM
      IF (L.EQ.0) RETURN  ! line was empty
      IF (ICOM.GT.0) THEN
c        CALL BOUNDS(line(K+L:LL),IP,I)
        IP=1
        CALL FINDPAR(line,2,IP,L) ! arguments, if any, start at IP
  !      write(*,*) 'linp 3'
c        IP=IP+K+L-1
        LINPEXECSTR=line(IP:LL)   ! return string with arguments

        NPAR=0
        DO WHILE (L.GT.0) ! get nummber of arguments
          K=1
  !        write(*,*) 'linp 4 IP=',IP,' LL=',LL,' NPAR=',NPAR,' K=',K
          CALL FINDPAR(line(IP:LL),NPAR+1,K,L)
  !        write(*,*) 'linp 4 '
          IF (L.GT.0) NPAR=NPAR+1
        ENDDO
      ELSE IF (ICOM.EQ.0) THEN
         IF (line(K:K+L-1).eq.'?') THEN
           I=0
           K=1
           CALL FINDPAR(line,2,K,L) ! is there a command name as an argument?
  !         write(*,*) 'linp 5'
           IF (L.GT.0) I=GETICOM(line(K:K+L-1))
           IF (I.GT.0) THEN ! if yes, print the hint for it
              IF (linp_h(I) .ne.' ') then
               WRITE(linp_out,*) linp_c(I)//'  '//linp_h(I)
              ENDIF
           ELSE ! otherwise print hints for all
               DO J=1,linp_nc
                 IF (linp_h(J).ne.' ') then
                   WRITE(linp_out,*) linp_c(J)//'  '//linp_h(J)
                 ENDIF
               ENDDO
               WRITE(linp_out,*) 'LIST   list values'
               WRITE(linp_out,*) 'QUIT   quit this menu'
           ENDIF
         ELSE
           NPAR=0
           CALL MKUPCASE(line(K:K+L-1))
           IF (INDEX('QUIT',line(K:K+L-1)).GT.0) THEN
             LINPEXECSTR='QUIT'
           ELSE IF (INDEX('LIST',line(K:K+L-1)).GT.0) THEN
             LINPEXECSTR='LIST'
           ELSE IF (INDEX('PAUSE',line(K:K+L-1)).GT.0) THEN
             LINPEXECSTR='PAUSE'
             CALL DOPAUSE
           ELSE IF (INDEX('ECHO',line(K:K+L-1)).GT.0) THEN
             LINPEXECSTR='ECHO'
             CALL DOECHO
           ELSE
              LINPEXECSTR=' '
              if (len_trim(SCOMM)>0) then
                WRITE(linp_out,*) 'Unknown command: '//trim(SCOMM)
              endif
              WRITE(linp_out,*) 'type ? for help'
           ENDIF
         ENDIF
       ENDIF
       END

C--------------------------------
      SUBROUTINE DOPAUSE
c stop & wait for enter
c--------------------------------
      implicit none
      INCLUDE 'linp.inc'
      CHARACTER*1 CH

1     FORMAT(' press ENTER ...',$)
2     format(A)
      write(linp_out,1)
      read(*,2) CH
      write(linp_out,*)

      if (CH.eq.'q'.OR.CH.EQ.'Q') linp_eof=1
      end

C--------------------------------
      SUBROUTINE DOECHO
c copy input to output until END - useful in batch files
c--------------------------------
      implicit none
      INCLUDE 'linp.inc'
      CHARACTER*1024 COMM,S
      INTEGER*4 IS,IL
      LOGICAL*4 BK

2     format(A)

      BK=.FALSE.
      write(linp_out,*)
      DO WHILE(.NOT.BK)
        read(linp_in,2,END=20) COMM
        CALL BOUNDS(COMM,IS,IL)
        S=COMM(IS:IS+IL-1)
        CALL MKUPCASE(S)
        BK=(IL.EQ.3.AND.S(1:IL).EQ.'END')
        IF (.NOT.BK) write(linp_out,*) COMM(IS:IS+IL-1)
      ENDDO
      RETURN
20    write(linp_out,*) 'ECHO command: END statement is missing'
      CALL DOPAUSE
      END

C     ---------------------------------------------------
      CHARACTER*(*) FUNCTION GETCOM(ICOM)
c     returns ICOM-th command name
C     ---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 ICOM
      IF (ICOM.GT.0.AND.ICOM.LE.linp_nc) THEN
         GETCOM=linp_c(ICOM)
      ELSE IF (ICOM.EQ.-1) THEN
         GETCOM='QUIT'
      ELSE IF (ICOM.EQ.-2) THEN
         GETCOM='?'
      ELSE
         GETCOM=' '
      ENDIF
      END

C     ---------------------------------------------------
      CHARACTER*(*) FUNCTION GETHINT(ICOM)
C     ---------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      INTEGER*4 ICOM
      IF (ICOM.GT.0.AND.ICOM.LE.linp_nc) THEN
         GETHINT=linp_h(ICOM)
      ELSE IF (ICOM.EQ.-1) THEN
         GETHINT='Quit interpreter'
      ELSE IF (ICOM.EQ.-2) THEN
         GETHINT='Show hints'
      ELSE IF (ICOM.EQ.-3) THEN
         GETHINT='Ambiguous command'
      ELSE
         GETHINT=' '
      ENDIF
      END





