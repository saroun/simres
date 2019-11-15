         PROGRAM RESTRAX

C////////////////////////////////////////////////////////////////////////
C
C                 ************************************ 
C                 ***                              *** 
C                 ***        S I M R E S           *** 
C                 ***            (PWD)             ***
C                 ***   (C) J.Saroun & J.Kulda     *** 
C                 ***        ILL Grenoble          ***
C                 ***      evaluation version      ***
C                 ************************************ 
C
C Program for:
C
C * Calculation of three-axis spectrometers (TAS) resolution matrices,
C * Monte Carlo simulation of the TAS resolution functions,
C * Calculation and graphical representation of projections the 
C   resolution functions and ellipsoids and widths of vanadium, 
C   Bragg and phonon scans.  
C 
C This program was written by J.Saroun and J.Kulda (ILL Grenoble
C      and NPI Rez near Prague) combining two older programs:  
C 
C 1) RESCAL - provides all except of the graphics subroutines and 
C    the subroutine AFILL, which calculates the resolution matrix.  
C 
C    - written by M.  Hergreave & P.  Hullah, P.C.L.  London, July 1979 
C    - updated by P.  Frings, ILL Grenoble, October 1986 
C 
C 2) TRAX - provides a kernel replacing the subroutine AFILL and taking 
C    into account real dimensions of the spectrometer components and 
C    focusing by curved perfect and mosaic crystals.  
C 
C    - written by (*)M.Popovici, A.D.Stoica and I.Ionita, Institute for 
C    Nuclear Power Reactors, Pitesti, 1984-1986 
C    (*) now University of Missouri
C    
C    - described in:
C
C    [1] M.Popovici, A.D.Stoica and I.Ionita, J.Appl.Cryst. 20 (1987),90.
C    [2] M.Popovici et al., Nucl.Instrum. Methods A338 (1994), 99.
C 
C and adding new code
C 
C 3) M.C. integration routine, interface routines and PGPLOT graphics 
C 
C 4) M.C. simulation of the TAS resolution functions
C 
C    - written by J.Saroun and J.Kulda, Institute Laue-Langevin, Grenoble, 
C      and Nuclear Physics Institute, Rez near Prague, November 1995 
C      
C****************************************************************************
C *** For all additional information contact the authors:                 ***
C ***                                                                     ***
C ***                 kulda@ill.fr        saroun@ujf.cas.cz                           ***
C ***                                                                     ***
C****************************************************************************

C***********************************************************
C
C ONLY M.C. SIMULATION OF NEUTRON FLUX IN THIS VERSION !!!
C
C***********************************************************


C//////////////////////////////////////////////////////////////////////
C////
C////  This file contains:
C////  
C////  * MAIN block
C////  * SUBROUTINE LOGO
C////  * SUBROUTINE INTRAC(NVARS,NCMDS,ICMD)
C////  * SUBROUTINE FNDVAR(LINE,NUM,NPARS,IVNO)
C////  * SUBROUTINE GETNAM(LINE,NUM,IVAR)
C////  * SUBROUTINE CHRNUM(LINE,NUM)
C////  * INTEGER FUNCTION CHRNO(CHAR)
C////  * SUBROUTINE UNITS(F,CUNIT,IER)
C////  * SUBROUTINE BRAG(IOCONS)
C////  * SUBROUTINE PHON(IOCONS)
C////  * SUBROUTINE RESOL(IOCONS)
C////  * SUBROUTINE RESCON(IOCONS)
C////  * SUBROUTINE GETRO(ICONS)
C////  * BLOCK DATA
C//// 
C//////////////////////////////////////////////////////////////////////
                         

        PArameter(iocons=5,NNVARS=49)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
      
      INCLUDE 'ness_common.inc'
        
      REAL*8 KFIX
      DIMENSION BBB(100),A(4,4),AAA(NNVARS)
C
      CHARACTER CUNIT*5
      CHARACTER*72 ShellArg
C        CHARACTER*20 DEVSTR

C
        COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1,ALF2,
     1  ALF3,ALF4,BET1,BET2,BET3,BET4,AS,BS,CS,AA,BB,CC,AX,AY,AZ,BX,
     2  BY,BZ,QH,QK,QL,EN,DH,DK,DL,DE,GH,GK,GL,GMOD,
     3  ROMH,ROMV,ROAH,ROAV,SDI,SHI,DTHETA
        COMMON /MATRIX/A
        COMMON /ERROR/IERR
        COMMON /UNIT/CUNIT
        COMMON /EXTRA/XXC,YYC,F,W,WI,WF
        COMMON /LOOP/ILOOP,JLOOP,IUSE
        COMMON /DATS/BBB
        COMMON /ARG/AGM(40)
        COMMON /BATCH/ IUEXT,IUOUT
        INTEGER*4 ISEED,IRND
        COMMON /RNDGEN/ ISEED,IRND
        character*72 extname,outname,S,S1
      CHARACTER*63 HLP(52),HLPOPT(18)      
      COMMON/HELP/ HLP,HLPOPT
      CHARACTER*30 DEVSTR
      COMMON /GRAPH/ DEVSTR
      LOGICAL*4 CHECKPARAM
      
      CHARACTER*60 CFGNAME
      CHARACTER*128 CFGPATH
      COMMON/CONFNAME/ CFGNAME,CFGPATH
        
      EQUIVALENCE (AAA(1),DM)
      REAL*8 HMOS_DIST !,ERF
      EXTERNAL HMOS_DIST

      CALL ERF_INIT(HMOS_DIST,-6.D+0,6.D+0)
      CALL GETENV('PGPLOT_DEV',DEVSTR)

c        DO WHILE(.TRUE.)
c          READ(*,*) Z
c          Z1=ERF(Z,0)
c          IF(Z1.GT.0.5) THEN
c             Z2=-1.D0*ERF(ERF(-Z,0),1)
c
c          ELSE
c             Z2=ERF(Z1,1)
c          ENDIF
c          WRITE(*,*) Z,Z1,Z2
c        ENDDO
c        STOP



C-----------------------------------------------------------------------
C
        
C Handle command-line options
        idbg=0
        IRND=0
        IOPT=1  ! automatic optimization by default
        NORMMON=0 ! constant monitor efficiency (NOT ~ 1/k)
        MDIST=0
        M=IARGC()
c        I=IGETARG(1,S,30)
        DO I=1,M
           call GETARG(I,S)
c           write(*,*) S
           if (S(1:5).eq.'-dir=') THEN
             CALL BOUNDS(S,IS,IL)
             IS=IS+5
             IL=IL-5
             IF(IL.GT.0) THEN
               CFGPATH=S(IS:IS+IL-1)
               write(*,*) 'dir='//CFGPATH(1:IL)
             ELSE
               CFGPATH=' '
             ENDIF  
           else if (S(1:2).eq.'-d') then
              IDBG=2
              read(S(3:3),*,err=11) J
              if(J.NE.0) IDBG=J
           endif   
11           if (S(1:2).eq.'-s') THEN
              read(S(3:30),*) J
              if(J.NE.0) ISEED=ABS(J) 
           endif
           if (S(1:2).eq.'-t') THEN
              read(S(3:30),*) J
              CALL RAN1SEED(ISEED)
              write(*,*) 'Test of the random number generator:'
              call RAN1TEST(J,1000000*J)
              STOP
           endif
           if (S(1:4).eq.'-flx') THEN
              S1=S(5:30)
              CALL READ_FLUX(S1)
           endif
           if (S(1:4).eq.'-flh') THEN
              read(S(5:30),*) FLXH
              FLXH=FLXH*PI/180
           endif
           if (S(1:4).eq.'-flv') THEN
              read(S(5:30),*) FLXV
              FLXV=FLXV*PI/180
           endif
           if (S(1:3).eq.'-RB') THEN
              read(S(4:30),*) CBAR
              write(*,*) 'Right barrier [mm]: ',CBAR
           endif
           if (S(1:3).eq.'-MX') THEN
              read(S(4:30),*) CMX
              write(*,*) 'crystal shift x [mm]: ',CMX
           endif
           if (S(1:6).eq.'-Voigt') THEN
              MDIST=1
              write(*,*) 'pseudo-Voigt mosaic distribution'
           else if (S(1:7).eq.'-Lorenz') THEN
              MDIST=2
              write(*,*) 'Lorenzian mosaic distribution'
           else if (S(1:4).eq.'-Uni') THEN
              MDIST=3
              write(*,*) 'Uniform mosaic distribution'
           else 
              MDIST=0
              write(*,*) 'Gaussian mosaic distribution'
           endif
           if (S(1:5).eq.'-ran1') THEN
              IRND=1
              write(*,*) 'Numerical Recipes RAN1 generator'
           endif
           if (S(1:5).eq.'-rand') THEN
              IRND=2
              write(*,*) 'System random number generator'
           endif
           if (S(1:6).eq.'-noopt') THEN
              IOPT=0
              write(*,*) 'No automatic sampling optimization'
           endif
           if (S(1:5).eq.'-nmon') THEN
              NORMMON=1
              write(*,*) 'Incident intensities ~ 1/ki'
           endif
           if (S(1:6).eq.'-cross') then
              ISQOM=0
             read(S(7:7),*,err=12) J
              if(J.GT.0.AND.J.LE.3) ISQOM=J
              if (ISQOM.EQ.1) then
                write(*,*) 'SCAN with antif. magnon cross-section'
              else if (ISQOM.EQ.2) then
                write(*,*) 'SCAN with Vanadium sample'
              endif  
           endif   
12           if (S(1:4).eq.'-log') THEN
              LOGFILE=S(5:30)
              write(*,*) 'Events logged in '//LOGFILE
10            format('Log events between [min max]: ',$)
              write(*,10)
              read(*,*) LOGMIN,LOGMAX
           endif
           if (S(1:5).eq.'-help'.or.S(1:1).eq.'?') THEN
              DO J=1,18
                 write(*,*) HLPOPT(J)
              ENDDO
              STOP            
           endif
        ENDDO   
        
C Initialize random number generator        
        CALL RAN1SEED(ISEED)
        
        CALL SETRESPATH(CFGPATH)  ! set default path for configuration
        
        NVARS=NNVARS
        NCMDS=43
        CALL LOGO        
2000    format(a30)
2001    format(' batch file  : ',$)
2004    format(' output file : ',$) 
        write(*,2001)   
        read(*,2000) extname
        if((extname(1:1).ne.' ').and.(extname(1:1).ne.char(0))) THEN
           IUEXT=10 
           
           CALL OPENRESFILE(extname,IUEXT,IERR,1) 
           IF (IERR.NE.0) GOTO 2002
c           OPEN (UNIT=IUEXT,ERR=2002,NAME=extname, READONLY,TYPE='OLD')
           write(*,2004)
           read(*,2000) outname
           if((outname(1:1).ne.' ').and.(outname(1:1).ne.char(0))) THEN             
             IUOUT=11
             OPEN(UNIT=IUOUT,ERR=2005,NAME=outname,STATUS='UNKNOWN')
             write(IUOUT,*) 'RESTRAX - batch job '//extname
c             CLOSE(IUOUT)
           else
              IUOUT=6
           endif   
        else
           IUEXT=5
           IUOUT=6
        endif 
        goto 2003     
2002    write(*,*) 'Cannot open batch file !'
        IUEXT=5
        IUOUT=6
        goto 2003
2005    write(*,*) 'Cannot open output file !'
        IUOUT=6        
2003    continue           
 
        CALL SETCFG(0)
        CALL SET_CRYST('Ge 111  ','Ge 111  ')
        if(IUOUT.NE.6) CLOSE(IUOUT)
        
C-----------------------------------------------------------------------
C   OBTAIN UNITS (MEV OR THZ) TO OBTAIN VALUE OF "F".
C
        CALL UNITS(F,CUNIT,IER)
C *** Set IUSE=1 to avoid asking for a *.RES file
C   IUSE=1
        IUSE=-1
C-----------------------------------------------------------------------
C
C   INTERACTIVE COMMAND AND PARAMETER INTERPRETER.
C
1     CONTINUE
      if(IUOUT.NE.6) then
        OPEN(UNIT=IUOUT,NAME=outname,ACCESS='APPEND',STATUS='UNKNOWN')
      endif


      IERR=0                        
      DO 3 I=1,NVARS
3           BBB(I)=AAA(I)
      CALL INTRAC(NVARS,NCMDS,ICMD)
      DO 2 I=1,NVARS
        AAA(I)=BBB(I)
2     CONTINUE
c      write(*,*) 'ICMD: ',ICMD      
      
      IF (.NOT.CHECKPARAM()) GOTO 111
      
      CALL RECLAT
      IF (IERR.NE.0) THEN
        WRITE (IOCONS+1,901) IERR, 'RECLAT'
        GOTO 111
      ENDIF
c      write(*,*) 'RECLAT O.K.'      
C   RECIPROCAL LATTICE.
        CALL TRNVCT
      IF (IERR.NE.0) THEN
        WRITE (IOCONS+1,901) IERR, 'TRNVCT'
        GOTO 111
      ENDIF
c      write(*,*) 'TRNVCT O.K.'      
     

C///  calls subroutine THRAX to calculate the resolution matrix  A ///  
    
      CALL THRAX(1) 
      IF (IERR.NE.0) THEN
        WRITE (IOCONS+1,901) IERR, 'THRAX'
        IF (IERR.NE.1) GOTO 111
      ENDIF
c      write(*,*) 'THRAX O.K.'      
      
      CALL TRANSMAT            ! calculates transformation marices

c      write(*,*) 'TRANSMAT O.K.'      

        
C   TRANSFORM VECTORS Q, D & G.
C   CALCULATE  RESOLUTION MATRICES WITHOUT OUTPUT ON DISPLAY
        IO = 0
        CALL RESOL(IO)
        
c      write(*,*) 'RESOL O.K.' , ICMD    
        
C
901        FORMAT( ' RESTRAX: ERROR ',I4,' IN ',A)
  
  
c        WRITE(*,*) 'COMMAND: ',ICMD
c        IF (ICMD.EQ.4) CALL NESS_LOOP
        IF (ICMD.EQ.2)  CALL PLOTOUT(5)
        IF (ICMD.EQ.3)  CALL SET_3AX(12)
        IF (ICMD.EQ.4) CALL SET_3AX(13)
        IF (ICMD.EQ.5) CALL PLOTOUT(0)
        IF (ICMD.EQ.6) CALL PRINTOUT
        IF (ICMD.EQ.7) CALL GETRO(IUOUT-1)
        IF (ICMD.EQ.8) CALL SET_3AX(1,IUOUT-1)
        IF (ICMD.EQ.10) CALL SAM_FLUX(5)
        IF (ICMD.EQ.11) CALL SAM_FLUX(1)
        IF (ICMD.EQ.12) CALL SETCFG(1)
        IF (ICMD.EQ.13) CALL TYPECFG
        IF (ICMD.EQ.14) CALL SET_3AX(2)
        IF (ICMD.EQ.15) CALL SET_3AX(3)                 
        IF (ICMD.EQ.16) CALL SET_3AX(4)                 
        IF (ICMD.EQ.17) CALL SET_3AX(5)                 
        IF (ICMD.EQ.18) CALL SET_3AX(6)                 
        IF (ICMD.EQ.19) CALL SET_DEVICE                 
        IF (ICMD.EQ.20) CALL SETVAR(1)                 
        IF (ICMD.EQ.21) CALL SETVAR(2)                
        IF (ICMD.EQ.22) CALL SET_3AX(7)                 
        IF (ICMD.EQ.23) CALL SET_3AX(8)                 
        IF (ICMD.EQ.24) THEN
            SHELLARG=' '
            CALL DOSHELL(SHELLARG)
        ENDIF        
        IF (ICMD.EQ.25) CALL SAM_FLUX(0)
        IF (ICMD.EQ.26) CALL BENCH
        IF (ICMD.EQ.27) CALL SET_3AX(9)
        IF (ICMD.EQ.28) CALL PLOTOUT(1)
        IF (ICMD.EQ.29) CALL PLOTOUT(2)
        IF (ICMD.EQ.30) CALL PLOTOUT(3)
        IF (ICMD.EQ.31) CALL SAM_FLUX(2)
        IF (ICMD.EQ.32) CALL ROCK(1)
        IF (ICMD.EQ.33) CALL GETROOPT(IUOUT-1)
        IF (ICMD.EQ.34) CALL GETROOPTMC(IUOUT-1)
        IF (ICMD.EQ.35) CALL SAM_FLUX(3)
        IF (ICMD.EQ.36) CALL SAM_FLUX(4)
        IF (ICMD.EQ.37) CALL PLOTOUT(4)
        IF (ICMD.EQ.38) CALL SCAN_CHI(0)
        IF (ICMD.EQ.39) CALL SET_3AX(10)
        IF (ICMD.EQ.40) CALL SET_3AX(11)
        IF (ICMD.EQ.41) CALL SCAN_THETA
        IF (ICMD.EQ.42) CALL ROCK(2)
        IF (ICMD.EQ.43) CALL SET_3AX(21)
                    
        IF (ICMD.NE.5.AND.ICMD.NE.6) WRITE(IUOUT,501)
501        FORMAT(1X,70('-'))

111     if(IUOUT.NE.6) then
        CLOSE(IUOUT)

c        pause
      endif
        GOTO 1
        END

C-----------------------------
      SUBROUTINE LOGO
C----------------------------- 
      IMPLICIT NONE
      INCLUDE 'config.inc'
      INTEGER*4 IUEXT,IUOUT
      COMMON /BATCH/ IUEXT,IUOUT

1     FORMAT(2x,'-----------------------------------------------',/,
     &       2x,'S I M R E S  - Monte Carlo ray-tracing ',/,
     &       2x,'Version: ',a40,/,
     &       2x,'Build:   ',a40,/,
     &       2x,'-----------------------------------------------',/,
     &       2x,'(C) J.Saroun & J.Kulda',/,
     &       2x,'ILL Grenoble, NPI Rez near Prague',/,
     &       2x,' ',/,
     &       2x,'Using (in part):',/,
     &       2x,'RESCAL by M.Hargreave, P.Hullah & P. Frings',/,
     &       2x,'TRAX   by M.Popovici, A.D.Stoica & I.Ionita',/,
     &       2x,'-----------------------------------------------',/,
     &       2x,'type ? for command list',/,/)
     

      write(iuout,1) PACKAGE_VERSION,PACKAGE_DATE  
      return
      end

C

C
        SUBROUTINE INTRAC(NVARS,NCMDS,ICMD)
C***************************************************************
C*                                                             *
C*   INTRAC:                                                   *
C*      INTERACTIVE COMMAND INTERPRETER ROUTINE.               *
C*      THIS ROUTINE ACCEPTS CHARACTER DATA TYPED IN AT THE    *
C*      USER'S TERMINAL AND INTERPRETS  THESE CHARACTERS AS    *
C*      ALTERATIONS  TO  PARAMETER VALUES, INBUILT COMMANDS    *
C*      OR  EXTERNALLY  DEFINED  COMMANDS.  IT CHANGES  THE    *
C*      VALUES OF PARAMETERS, PERFORMS INBUILT COMMANDS AND    *
C*      EXITS BACK TO THE CALLING PROGRAMME UPON REACHING A    *
C*      COMMAND WHICH WAS DEFINED IN THE BLOCK DATA ROUTINE    *
C*      IN WHICH ARE ALSO DEFINED THE PARAMETER NAMES.         *
C*                                                             *
C*    WRITTEN BY    PETER H.C. HULLAH,                         *
C*                     DIVN. OF SCIENCE,                       *
C*                       P. C. L.,                             *
C*                         115, NEW CAVENDISH ST.,             *
C*                           LONDON. W1M 8JS.                  *
C*                               (TEL 01-637 7790)             *
C*                                                             *
C*                                                JULY 1978    *
C*                                    UPDATED NOVEMBER 1979    *
C*                                 2ND UPDATE  OCTOBER 1980    *
C***************************************************************
        PArameter(iocons=5)
C
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL*8 NAM(100),IVAR
        INTEGER U
        DIMENSION DAT(100),LINE(120),LLINE(9,120),RET(40)
        CHARACTER*60 FLNAM,DATNAME,DIRNAME
        character*120 CLINE
        CHARACTER*79 STRDUM
        LOGICAL USED(40),NOSAVE,SCAN
        COMMON/CHRT/NOS,RET
        COMMON/LOOP/ILOOP,JLOOP,IUSE
        COMMON/CHNM/NAM
        COMMON/DATS/DAT
        COMMON/ARG/AGM(40)
      CHARACTER*63 HLP(53),HLPOPT(18)      
      COMMON/HELP/ HLP,HLPOPT
        COMMON/USED/ USED
        COMMON/BATCH/ IUEXT,IUOUT        
        DATA SCAN/.FALSE./
        DATA FLNAM,DATNAME,DIRNAME/' ',' ',' '/
        
        
        U=IUEXT
        NPARS=NVARS+NCMDS
              IF(IUSE.EQ.-1)GO TO 40
        IF(IUSE.EQ.0)GO TO 1
        IF(ILOOP.EQ.JLOOP)GO TO 6
        KLOOP=ILOOP+1
        GO TO 11
C---------------------------------------------------------------
C     INITIALISATION.
C       SCAN=.TRUE.

40         WRITE(IOCONS+1,72)
        READ(U,73) DIRNAME
        IF(IUOUT.NE.6) THEN
           WRITE(IUOUT,72)
           WRITE(IUOUT,*) DIRNAME
        ENDIF 
C          U=23
C        JLOOP=4
1         IUSE=1

        WRITE(IOCONS+1,2)
2        FORMAT(' Name of a parameter or data file - ',$)
        READ(U,3)DATNAME        
3        FORMAT(A)
        IF(IUOUT.NE.6) THEN
           WRITE(IUOUT,2)
           WRITE(IUOUT,*) DATNAME
        ENDIF   
101        DO 102 I=1,100
102        DAT(I)=0.0
        NOSAVE=.FALSE.

C////// added to read spectra files:  ///////
        call ReadDatFile(DIRNAME,DATNAME,ie)
            
        if(ie.eq.0) then
           FLNAM=' '
        else if(index(DATNAME,'.scn').ne.0) then
           write(IOCONS+1,*) ' Data file not found'
           goto 1
        else      
        FLNAM=DATNAME
C///////////////////////////////////////////    

        IF(FLNAM(1:1).EQ.' ') GOTO 6
        
        CALL OPENRESFILE(FLNAM,24,IST,1)
        IF (IST.NE.0) GOTO 601
c       OPEN (UNIT=24,ERR=601,IOSTAT=IST,NAME=FLNAM, READONLY,TYPE='OLD')
        GOTO 603
601        IF (IST.EQ.29) THEN
        WRITE(IOCONS+1,501) FLNAM
        WRITE(IOCONS+1,*) ' INTRAC: This File Doesn''t Exist ...' 
     *   //'Try Another '
501        FORMAT(1X,20A)
        ENDIF
        GOTO 1
603        CONTINUE
        DO 4 I=1,NVARS
4        READ(24,*,ERR=66,END=66) DAT(I)
C5        FORMAT(G)
66      CLOSE(UNIT=24)

        endif

C---------------------------------------------------------------
6        KLOOP=1
        CALL SET_CRYST(' ',' ')

C    LOOP IS INITIALLY 1.
        IF(SCAN) GOTO 7
C        U=IOCONS
        JLOOP=1
C---------------------------------------------------------------
C    LOOP FOR GETTING  NEW LINES OF DATA.
7        DO 10 ILOOP=KLOOP,JLOOP
        IF (SCAN) GOTO 41
907        WRITE(IUOUT,8)
8        FORMAT(' ResTrax>',$)
41      READ(U,73) CLINE
C        WRITE(IUOUT,*) CLINE(1:60) 
        IF (IUOUT.NE.6) WRITE(IUOUT,*) CLINE(1:60) 
C          WRITE(IUOUT,*) '********************************************'
C        ENDIF   
        READ(CLINE,9)(LLINE(ILOOP,J),J=1,120)
9        FORMAT(120A1)
        IF(U.NE.IOCONS) GOTO 10
        IF(LLINE(ILOOP,1).NE.' ') GOTO 10
        READ(IOCONS,9)SCR
        GOTO 907
10        CONTINUE
        SCAN=.FALSE.
C===============================================================
C    LOOP FOR PROCESSING DATA.
11        DO 14 ILOOP=KLOOP,JLOOP
        DO 12 J=1,120
        LINE(J) = LLINE(ILOOP,J)
12        CONTINUE
        NUM=1
C     NEW NAME.
13        CALL FNDVAR(LINE,NUM,NPARS,IVNO)
        IF(IVNO.LE.0)GO TO 14
        IF(IVNO.LE.NVARS)GO TO 17
        IF(IVNO.LE.NPARS)GO TO 15
        GO TO (21,25,27,31,32,34,35,71)IVNO-NPARS
14        CONTINUE
        GO TO 6
C===============================================================
C    NAME IS COMMAND.
15        ICMD=IVNO-NVARS
        CALL CHRNUM(LINE,NUM)
C     OBTAIN ARGUMENT VALUES.
        DO 16 I=1,40
        AGM(I)=RET(I)
16        CONTINUE
        RETURN
C---------------------------------------------------------------
C     NAME IS PARAMETER NAME
17        CALL CHRNUM(LINE,NUM)
C     'CHRNUM' GETS THE VALUES OF THE DATA FOLLOWING THE NAME.
C      THESE DATA ARE PLACED IN THE ARRAY RET, THEIR NUMBER
C      BEING PLACED IN NOS.
        NOSMIN=MIN0(NOS,NVARS+1-IVNO)
        IF(NOS.EQ.0)NOSMIN=1
        DO 18 I=1,NOSMIN
        IF(.NOT.USED(I))GO TO 18
        J=IVNO+I-1
        DAT(J)=RET(I)
18        CONTINUE
        WRITE(IUOUT,19)(NAM(IVNO+I-1),DAT(IVNO+I-1),I=1,NOSMIN)
19        FORMAT(3(' ',A5,' = ',G14.7,1X))
        IF(NOSMIN.LT.NOS)WRITE(IOCONS+1,20)
20        FORMAT(' INTRAC: Too Many Data Items. Excess Ignored')
        NOSAVE=.TRUE.
        IF(NUM.GE.121)GO TO 14
        GO TO 13
C---------------------------------------------------------------
C    NAME IS "LOOP" COMMAND.
21        WRITE(IOCONS+1,*) 'Command not available'
        goto 6
c        CALL CHRNUM(LINE,NUM)
c        IF(RET(1).GT.9)GO TO 23
c        IF(RET(1).LT.1)GO TO 6
c        JLOOP=RET(1)
c        WRITE(IOCONS+1,22)JLOOP
c22        FORMAT(' INTRAC: Input ',I1,' Lines')
c        GO TO 7
c23        WRITE(IOCONS+1,24)
c24        FORMAT(' INTRAC: Max Loop Of 9 - Retype Line')
c        GO TO 6
C---------------------------------------------------------------
C     NAME IS "LIST" COMMAND.
25        I1=1
        I2=NVARS
        CALL FNDVAR(LINE,NUM,NVARS,IVNO)
        IF(IVNO.LE.0)GO TO 26
        IF(IVNO.GT.NVARS)GO TO 26
        I1=IVNO
        I2=IVNO
        CALL FNDVAR(LINE,NUM,NVARS,IVNO)
        IF(IVNO.LE.0)GO TO 26
        IF(IVNO.GT.NVARS)GO TO 26
        I1=MIN0(IVNO,I1)
        I2=MAX0(IVNO,I2)
C  **** NEXT LINES SPECIFIC TO USE OF INTRAC WITH RESCAL
26    IF(I1.EQ.1.AND.I2.GE.NVARS) GO TO 50
          WRITE(IUOUT,19)(NAM(I),DAT(I),I=I1,I2)
        GO TO 14

C///  writes also sample position (important for Monte Carlo)   
   50 nos=0
      WRITE(IUOUT,51) (NAM(I),DAT(I),I=1,NVARS)           !,(RET(I),I=1,3)
   51 FORMAT(' ',2(A4,' = ',F10.5,1X)/(' ',3(A4,' = ',F10.2,1X)/),
     1 (' ',3(A4,' = ',F10.0,1X)/),(' ',(A4,' = ',F10.5,1X,
     2 A4,' = ',F10.0,1X)/),2(' ',4(A4,' = ',F10.2,1X)/),
     3 4(' ',3(A4,' = ',F10.4,1X)/),2(' ',4(A4,' = ',F10.4,1X)/),
     4 ' ',4(A4,' = ',F10.4,1X)/(' ',4(A4,' = ',F10.4,1X)/),
     5 ' ',2(A4,' = ',F10.2,1X)/(' ',A4,' = ',F10.2,1X))        
     
      !/' sample position: ',3(2x,F5.1),' mm')
C///  writes also sample position etc. (important for Monte Carlo)   
       nos=0
       call SET_3AX(1)
C       call SET_3AX(2)
C       call SET_3AX(3)
       call SET_3AX(4)
       call SET_3AX(5)
       call SET_3AX(6)
      GO TO 14
C---------------------------------------------------------------
C     NAME IS "SAVE" COMMAND.
27        CALL GETNAM(LINE,NUM,IVAR)
        IF(IVAR.EQ.0.AND.FLNAM(1:1).NE.' ') GO TO 28
        WRITE(IOCONS+1,2002)
2002    FORMAT(' Give Filename To Save Data -',$)
        READ(U,2003)FLNAM
2003    FORMAT(60A)
28        OPEN(UNIT=24,NAME=FLNAM,STATUS='UNKNOWN',ERR=2004)
          DO 29 I=1,NVARS
29        WRITE(24,*) DAT(I)
        CLOSE(UNIT=24)
        WRITE(IOCONS+1,30)FLNAM
30        FORMAT('  INTRAC : Parameters Saved In File ',20A/)
        NOSAVE=.FALSE.
        GO TO 14
2004        WRITE(IOCONS+1,*) 'Cannot save parameters in '//FLNAM
        GO TO 14
C---------------------------------------------------------------
C    NAME IS "FILE" COMMAND.
31        IF (NOSAVE) GOTO 67
        CALL GETNAM(LINE,NUM,IVAR)
        GO TO 1
C---------------------------------------------------------------
C   NAME IS "EXIT" COMMAND.
32        IF(.NOT.NOSAVE) THEN
         IF(IUEXT.NE.5) CLOSE(IUEXT)
         IF(IUOUT.NE.6) CLOSE(IUOUT)
         CALL NESSEND                !  NESSEND must be called to deallocate
         STOP '  -> End of ResTrax'  !  dynamical memory used to store M.C.events
        ENDIF    
67        NOSAVE=.FALSE.
        WRITE(IOCONS+1,33)           ! 
33        FORMAT(1X,'!!  INTRAC: New Parameter Values Not SAVED !!'/
     1  ' Retype Command Or Save')
        GO TO 14
C-----------------------------------------------------------------------------
C    NAME IS "EXFF" COMMAND.
34        CALL NESSEND               !  NESSEND must be called to deallocate
        IF(IUEXT.NE.5) CLOSE(IUEXT)
        IF(IUOUT.NE.6) CLOSE(IUOUT)
        STOP ' End of ResTrax'     !  dynamical memory used to store M.C.events
C------------------------------------------------------------------
C    NAME IS "HELP" COMMAND.
35      CONTINUE
        IHUNIT=22
        IEROP=0
        OPEN(UNIT=IHUNIT,ERR=940,NAME='simres.hlp',
     1  STATUS='OLD',ACTION='READ',IOSTAT=IOSTO)
        GOTO 37
940     IHUNIT=0
c        WRITE(IOCONS+1,*) ' File simres.hlp not found. Internal'//
c     *  ' help is used.'
C/// read of internal help is forbiden (I/99, J.S.)
        WRITE(IOCONS+1,*) 
     *  ' File simres.hlp not found in current directory.'
c        GOTO 14
        
        DO 941 I=1,53
941      WRITE(IOCONS+1,1010) HLP(I)
        GOTO 14        
37      READ(22,1000,ERR=39,END=39)STRDUM
1000    FORMAT(A)
1010    FORMAT(1X,A)
        i=LEN(STRDUM)
        do 1011 while ((STRDUM(I:I).eq.' ').and.(i.gt.0))
1011          i=i-1
        WRITE(IOCONS+1,1010) STRDUM(1:i)
        GOTO 37
        
C---------------------------------------------------------------------
C    NAME IS "DIR" COMMAND.
71      CONTINUE        
         WRITE(IOCONS+1,72)
72        FORMAT(' Path to the spectra files - ',$)
        READ(U,73) DIRNAME
73        FORMAT(A)       
        GO TO 14
C----------------------------------------------------------------
        
39        CONTINUE
        IF (IEROP.NE.0)WRITE (IOCONS+1,*) ' I/O error! ' 
        IF (IHUNIT.EQ.22) CLOSE(22)
38          GOTO 14
C---------------------------------------------------------------
        END
C
C
C
        SUBROUTINE FNDVAR(LINE,NUM,NPARS,IVNO)
C***********************************************************************
C     FINDS THE IDENTITY OF THE NEXT NAME ON THE LINE.
C     "IVNO" IS ITS POSITION IN THE ARRAY "NAM".
C***********************************************************************
        PARAMETER(IOCONS=5)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        REAL*8 CMD(8),NAM(100),IVAR,E1,E2
        DIMENSION LINE(120)
C
C        CHARACTER*120 LINE,NAM*4(100),CMD*5(7),E1,E2
C
        COMMON/CHNM/NAM
        DATA CMD/'LOOP','LIST','SAVE','FILE','EXIT','EXFF','HELP','DIR'/
        IVNO=0
        CALL GETNAM(LINE,NUM,IVAR)
C     "GETNAM" FINDS THE VALUE OF THE NEXT WORD ON THE LINE.
C     THIS VALUE IS STORED IN "IVAR"
        IF(IVAR.EQ.0)GO TO 4
C     IF "IVAR" = ZERO THE REST OF THE LINE WAS BLANK
C     AND "IVNO" KEEPS THE VALUE ZERO.
        DO 1 I=1,NPARS
C     CHECK AGAINST USER DEFINED NAMES.
        IVNO=I
        DECODE(5,1000,IVAR)E1
        DECODE(5,1000,NAM(IVNO))E2
        IF(E1.EQ.E2)GO TO 4
1        CONTINUE
        DO 2 I=1,8
C     CHECK AGAINST INTERNAL COMMANDS.
        IVNO=NPARS+I
        DECODE(5,1000,CMD(I))E2
        IF(E1.EQ.E2)GO TO 4
1000    FORMAT(A5)
2        CONTINUE
c        WRITE(IOCONS+1,3)IVAR
c3        FORMAT(' FNDVAR: Variable "',A5,'" Unknown. Line Ignored 
c     * From Here')
        IVNO=-5
C     IF "IVAR" WAS NOT A LEGAL NAME "IVNO" IS GIVEN A -VE VALUE.
4        RETURN
        END
C
        SUBROUTINE GETNAM(LINE,NUM,IVAR)
C***********************************************************************
C     GETS THE VALUE OF THE NEXT NAME ON THE LINE.
C     THIS IS PLACED IN "IVAR".
C***********************************************************************
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION LINE(120),KHAR(5)
        CHARACTER CH
        REAL*8 IVAR
        INTEGER DELIM(5)
        DATA DELIM/1H ,1H,,1H/,1H=,1H;/
        K=0
        IVAR=0
C---------------------------------------------------------------
C     SETTING UP "KHAR" WHICH WILL HOLD THE LETTERS, AS BLANKS.
        DO 1 I=1,5
1        KHAR(I)=DELIM(1)
C---------------------------------------------------------------
C     REMOVING PRECEDING BLANKS.
2        IF(LINE(NUM).NE.DELIM(1))GO TO 3
        NUM=NUM+1
        IF(NUM.GE.121)GO TO 8
C     IF "NUM" IS 121 THE LINE HAS ENDED WITHOUT A NAME FOUND.
        GO TO 2
C---------------------------------------------------------------
C     PUT LETTERS FROM LINE INTO "KHAR".
3        DO 4 I=1,5
        IF(LINE(NUM).EQ.DELIM(I))GO TO 5
C     CHECK FOR DELIMITER TO END NAME
4        CONTINUE
        K=K+1
        IF(K.GT.5)GO TO 5
C     ONLY 5 LETTERS MAY BE PACKED.
        KHAR(K)=LINE(NUM)
C     PUT LETTER IN CORRECT ELEMENT OF "KHAR".
        NUM=NUM+1
        IF(NUM.GE.121)GO TO 5
C     THE WORD MAY BE TERMINATED BY THE END OF THE LINE.
        GO TO 3
C---------------------------------------------------------------
C     PACK CHARACTERS INTO A SINGLE 36-BIT WORD.
5       CONTINUE 
998     FORMAT(A1)      
        DO 999 I=1,5
            WRITE(CH,998) KHAR(I)
            ICH=ICHAR(CH)
            IF (ICH.GT.92.AND.ICH.LT.123) ICH=ICH-32
            CH=CHAR(ICH)
            READ(CH,998) KHAR(I)
999     CONTINUE    
   
        ENCODE(5,1000,IVAR)KHAR
1000    FORMAT(5A1)
C---------------------------------------------------------------
8        RETURN
        END
C

C
C
        SUBROUTINE CHRNUM(LINE,NUM)
C***********************************************************************
C     ROUTINE FOR CONVERTING NUMBERS FROM CHARACTERS TO NUMERICS.
C***********************************************************************
        PArameter(iocons=5)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        INTEGER BLANK,COMMA,SLASH,EQUATE,SEMICO
        INTEGER ISIGN(2),LNUM,I,LINE(120),CHRNO,IEXP
        LOGICAL START,PT,USED(40)
        DIMENSION RET(40)
        COMMON/CHRT/NOS,RET
        COMMON/USED/ USED
        DATA BLANK/' '/,COMMA/','/,SLASH/'/'/,EQUATE/'='/,SEMICO/';'/
C---------------------------------------------------------------
C     INITIALISATION.
        NOS=0
        DO 1 I=1,40
        USED(I)=.FALSE.
1        RET(I)=0.
        USED(1)=.TRUE.
C===============================================================
C     REMOVE BLANKS AND DELIMITER
2        LNUM=NUM
        NUM=NUM+1
        IF(LNUM.GE.121)RETURN
        IF(LINE(LNUM).EQ.BLANK)GO TO 2
        IF(LINE(LNUM).EQ.COMMA)GO TO 3
        IF(LINE(LNUM).EQ.SLASH)GO TO 3
        IF(LINE(LNUM).EQ.EQUATE)GO TO 3
        IF(LINE(LNUM).EQ.SEMICO)GO TO 3
        NUM=LNUM
        GO TO 4
C---------------------------------------------------------------
C     REMOVE BLANKS AFTER DELIMITER
3        IF(NUM.GE.121)RETURN
        IF(LINE(NUM).NE.BLANK)GO TO 4
        NUM=NUM+1
        GO TO 3
C===============================================================
C     DECODE NUMBER
4        NUM=NUM-1
        NOS=NOS+1
        PT=.FALSE.
        START=.FALSE.
        ISIGN(1)=1
C---------------------------------------------------------------
C     GET NEXT CHARACTER
5        NUM=NUM+1
        I=CHRNO(LINE(NUM))
        IF(I.EQ.0)GO TO 15
        IF(I.LE.10)GO TO 8
        GO TO (6,7,7,10,10,14),I-10
C---------------------------------------------------------------
C     CHARACTER IS A DECIMAL POINT
6        IF(PT)GO TO 15
        START=.TRUE.
        PT=.TRUE.
        EXPON=1.
        GO TO 5
C---------------------------------------------------------------
C     CHARACTER IS SIGN
7        IF(START)GO TO 15
        IF(I.EQ.13)ISIGN(1)=-1
        START=.TRUE.
        GO TO 5
C---------------------------------------------------------------
C     CHARACTER IS A DIGIT
8        START=.TRUE.
        IF(PT)GO TO 9
        RET(NOS)=RET(NOS)*10+FLOAT(I-1)
        GO TO 5
9        EXPON=EXPON/10.
        RET(NOS)=RET(NOS)+FLOAT(I-1)*EXPON
        GO TO 5
C---------------------------------------------------------------
C     CHARACTER IS "E" OR "D"
10        IF(.NOT.START)GO TO 18
        IEXP=0
        NUM=NUM+1
C     GET NEXT CHARACTER
        ISIGN(2)=1
        I=CHRNO(LINE(NUM))
        IF(I.EQ.0)GO TO 16
        IF(I.NE.13.AND.I.NE.12)GO TO 11
C     THIS SECTION ONLY DEALS WITH A SIGN (+ OR -)
        IF(I.EQ.13)ISIGN(2)=-1
        NUM=NUM+1
C     IF CHARACTER WAS A SIGN NOW GET NEXT
        I=CHRNO(LINE(NUM))
        IF(I.EQ.0)GO TO 16
11        IF(I.GT.10)GO TO 16
C     THIS CHARACTER MUST BE NUMERIC
        IEXP=I-1
        NUM=NUM+1
C     GET NEXT CHARACTER
        I=CHRNO(LINE(NUM))
        IF(I.EQ.0)GO TO 16
        IF(I.GT.10)GO TO 12
C     IF THIS CHARACTER WAS NOT NUMERIC THE EXPONENT HAS ENDED
        IEXP=IEXP*10+I-1
12        RET(NOS)=RET(NOS)*10.**(ISIGN(2)*IEXP)
C     MULTIPLY NUMBER BY EXPONENT
        IF(I.GT.10)GO TO 13
        NUM=NUM+1
        I=CHRNO(LINE(NUM))
C     THE CHARACTER AFTER THE EXPONENT MUST BE A DELIMITER
13        IF(I.NE.16)GO TO 16
C---------------------------------------------------------------
C     END OF NUMBER (DELIMITER FOUND)
14        IF(.NOT.START)GO TO 2
        USED(NOS)=.TRUE.
        RET(NOS)=RET(NOS)*ISIGN(1)
        GO TO 2
C===============================================================
C     NON-NUMERIC CHARACTER FOUND
15        IF(START)GO TO 16
C     IF NUMBER HAS NOT STARTED CHARACTER IS LEGAL
18        NOS=NOS-1
        RETURN
16        WRITE(IOCONS+1,17)LINE(NUM),NUM
17        FORMAT(' ILLEGAL CHARACTER "',A1,'" AT POSITION ',I2/
     1        ' LINE IGNORED FROM HERE')
        NUM=121
        RETURN
        END

             INTEGER FUNCTION CHRNO(CHAR)
C***********************************************************************
C     FUNCTION TO GIVE POSITION OF ARGUMENT CHAR IN ARRAY KHAR
C
C***********************************************************************
        INTEGER CHAR,KHAR(20)
        DATA KHAR/1H0,1H1,1H2,1H3,1H4,1H5,1H6,1H7,1H8,1H9,
     1        1H.,1H+,1H-,1HD,1HE,1H ,1H,,1H/,1H=,1H;/
        DO 1 I=1,20
        IF(CHAR.EQ.KHAR(I))GO TO 2
1        CONTINUE
C     IF CHAR IS NOT FOUND IT IS ILLEGAL AND RETURN VALUE IS 0
        CHRNO=0
        RETURN
2        IF(I.GT.16)I=16
C     IF CHAR IS , / OR = IT IS A DELIMITER AND NO DISTINCTION
C      IS MADE BETWEEN THE THREE DELIMITERS OR BLANK IN THIS CASE
        CHRNO=I
        RETURN
        END

        SUBROUTINE UNITS(F,CUNIT,IER)
C***********************************************************************
C
C***********************************************************************
        PArameter(iocons=5)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        CHARACTER*5 CUNIT
        CHARACTER CHAR,CARCOUR
        COMMON/BATCH/ IUEXT,IUOUT
        IER=0 
2        WRITE(IOCONS+1,1)
1        FORMAT(' Energy In m[eV] Or T[Hz]? (meV) ',$)
        READ(IUEXT,5) CARCOUR
5        FORMAT(A)
        IASC=ICHAR(CARCOUR)
        IF (IASC.GT.92.AND. IASC.LT.123) IASC=IASC-32
        CARCOUR=CHAR(IASC)
        IF(CARCOUR.EQ.'T') THEN
        F=1.996854
        WRITE(CUNIT,503)
503        FORMAT('[THz]')
        ELSE
        F=.4826
        WRITE(CUNIT,501)
501        FORMAT('[meV]')
        ENDIF
        WRITE(IOCONS+1,6) CUNIT
6        FORMAT(' Units Are ',A)
        RETURN
        END
C
C
        SUBROUTINE BRAG(IOCONS)
C***********************************************************************
C
C***********************************************************************
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        PARAMETER(EPS1=1.D-5)
        REAL*8 KFIX
        DIMENSION A(4,4)
        CHARACTER*5 CUNIT
        COMMON DM,DA,ETAM,ETAA,ETAS,SM,SS,SA,KFIX,FX,ALF1
     1        ,ALF2,ALF3,ALF4,BET(4),DUMMA(3),DUMM(9),Q(3),EN,D(3),DE,
     2        G(3),GMOD,ROMH,ROMV,ROAH,ROAV,SDI,SHI

               COMMON/ARG/AGM(40)
        COMMON /EXTRA/ XXC,YYC,F,W,WI,WF
        COMMON /UNIT/CUNIT
        COMMON /MATRIX/A
        COMMON /S/S(3,3),SD(3,3),B(3),COSB(3)
        COMMON/TRSG/SX,SY,SXX,SYY,GX,GY,GXX,GYY
        LOGICAL LZER

        HUITLOG2=8.D0*LOG(2.D0)
        LZER=.FALSE.
        DO 11 I=1,4
        IF (ABS(A(I,I)).LT.EPS1) LZER=.TRUE.
11        CONTINUE

        IF (LZER) THEN
        WRITE (IOCONS+1,*) 
     1  ' BRAG: Sorry  Pb with matrix -> I am unable to divide by zero'
        RETURN
        ENDIF

        IF(AGM(1).GT.1.5)GO TO 4
C---------------------------------------------------------------------------
C   WIDTHS IN RECIPROCRAL ANGSTROMS
        DQX=SQRT(HUITLOG2/A(1,1))
        DQY=SQRT(HUITLOG2/A(2,2))
        DQZ=SQRT(HUITLOG2/A(3,3))
C   ENERGY AND VANADIUM WIDTHS.
        DEE=SQRT(HUITLOG2/A(4,4))
        DVN=(A(2,4)*A(1,1)-A(1,4)*A(1,2))**2/(A(1,1)*(A(2,2)*A(1,1)
     1        -A(1,2)**2))
        DVN=SQRT(HUITLOG2/(A(4,4)-A(1,4)**2/A(1,1)-DVN))
        WRITE(IOCONS+1,5)DQX,DQY,DQZ,CUNIT,DVN,DEE
5       FORMAT(' Bragg Widths (radial,tangential,vertical) [A-1]'/    
     1        ' DQR=',F9.5,' DQT=',F9.5,' DQV=',F9.5/
     1  ' '/
     2        ' Energy Widths (Vanadium, Bragg) ',A5/
     3        ' DVN=',F9.5,' DEE=',F9.5)
        RETURN
C---------------------------------------------------------------------------
C   WIDTH IN TERMS OF NUMBER OF STEPS.
4        SXY=A(1,1)*SXX*SXX+2.*A(1,2)*SXX*SYY+A(2,2)*SYY*SYY
        AB=SQRT(SXY)
        IF(AB.GT.EPS1) THEN
        W=2.3548/AB
        WRITE(IOCONS+1,13)W
13        FORMAT(' FWHM along (DH,DK,DL): ',F9.5,' [steps]')
        ELSE
        WRITE(IOCONS+1,601) D(1),D(2),D(3)
601        FORMAT(' BRAG: CHECK DH DK DL :',3(X,F9.3))
        ENDIF
        RETURN
        END

        SUBROUTINE RESOL(IOCONS)
C***********************************************************************
C
C***********************************************************************
C        PARAMETER(IOCONS=5)
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)
        DIMENSION A(4,4),ADA(4,4),B(4,4),AZ(4,4)
        CHARACTER*5 CUNIT
        COMMON/ARG/AGM(40)
        COMMON /UNIT/CUNIT
        COMMON /MATRIX/A
        COMMON /MCARLO/ ADA,B,AZ
        CHARACTER*2 IAX(4),FWHM*4
        CHARACTER*2 IAXP(4),STRVI*4,STRVF*4
        CHARACTER*54 MESOR(5)
        DATA IAX/'X','Y','Z','W'/,FWHM/'FWHM'/
        DATA IAXP/'X''','Y''','Z''','W'''/,STRVI/'VI ='/,STRVF/'VF ='/
        DATA MESOR/
     1        ' Resolution Volumes [A-3]: ',
     2        ' Resolution Matrix, x-axis along Q, [A-1]',
     3        ' Resolution Matrix, in recip. lattice coor., [r.l.u.]',
     4        ' Diagonalised Resolution Matrix, [r.l.u.]',
     5        ' Direction Cosines (w.r.t. reciprocal lattice)'/
        
        IND = 0
1       IND=IND+1
        IARG=IFIX(REAL(AGM(IND))+0.5)
        IF(IOCONS.NE.0) THEN
                GO TO (11,12,13,14,15),IARG
        ELSE
                CALL VIVF(VI,VF)
                CALL CONVRT(A,AZ)
                CALL DIAG(AZ,ADA,B)
                TMP=8.D0*LOG(2.D0)
                DO 331 I=1,4
331                ADA(I,I)= SQRT(TMP/ABS(ADA(I,I)))
        ENDIF
89      RETURN
C --------------------------------------------------------------------
C        CAS RESOL 1 ON CALCULE LES VOLUMES VI VF
C
11        CALL VIVF(VI,VF)
        WRITE(IOCONS+1,411) MESOR(1)(1:27),STRVI,VI,STRVF,VF
411        FORMAT(A,2(4X,A,G11.4))
        GO TO 1
CC
C --------------------------------------------------------------------
C        CAS RESOL 2 ON CALCULE LA MATRICE SUIVANT Q
C
12        WRITE(IOCONS+1,413) MESOR(2)(1:41),CUNIT,IAX
        WRITE(IOCONS+1,424) (IAX(I),(A(I,J),J=1,4),I=1,4)
        GO TO 1
C --------------------------------------------------------------------
C        CAS RESOL 3 ON CALCULE LA MATRICE SUIVANT THE Reciprocal Lattice

13        CALL CONVRT(A,AZ)
        WRITE(IOCONS+1,413)MESOR(3)(1:53),CUNIT,IAX
413        FORMAT(2A//5(10X,A))
412        FORMAT(A//5(10X,A))
        WRITE(IOCONS+1,424) (IAX(I),(AZ(J,I),J=1,4),I=1,4)
        GO TO 1
C --------------------------------------------------------------------
C        CAS RESOL 4 ON DIAGONALISE
C
14        CALL DIAG(AZ,ADA,B)
        IF (IOCONS.EQ.5) WRITE(IOCONS+1,413) MESOR(4)(1:41),CUNIT,IAXP
        DO 31 I=1,4
          DO 311 K=1,4
            IF (ADA(I,K).LT.1.D-20) ADA(I,K)=0.D0
311          CONTINUE
          IF (IOCONS.EQ.5) WRITE (IOCONS+1,424) IAXP(I),
     *    (ADA(I,K),K=1,4)
31        CONTINUE
424        FORMAT(2X,A,4E12.4)
C
C *** At Cooper&Nathans the res. function is exp(-.5*MijXiXj) !!!!!!
C 
500     FORMAT(1X,70('-'))
        TMP=8.D0*LOG(2.D0)
        DO 33 I=1,4
33        ADA(I,I)= SQRT(TMP/ABS(ADA(I,I)))
C       WRITE(IOCONS+1,500)
           WRITE(IOCONS+1,412) MESOR(5),IAX,FWHM
        DO 35 I=1,4
        WRITE (IOCONS+1,444) IAXP(I),(B(J,I),J=1,4),ADA(I,I)
35        CONTINUE
444        FORMAT (4(2X,A,5F12.5))
         GO TO 1
C --------------------------------------------------------------------
C        case RESOL 5 - covariance and resolution matrices 
C         by Monte Carlo simulation
15         CALL NESS(0,0.D0)
         GO TO 1
         
        END
C

C
C
        SUBROUTINE RESCON(IOCONS)
C***********************************************************************
C        RISO PROGRAM FOR FINDING CONTRIBUTIONS TO WIDTH OF A CONST-Q
C        SCAN OF TRIPLE AXIS SPECTROMETER.
C        A* AND B* - ORTHOGONAL RECIPROCAL LATTICE DIMENSIONS IN SCAT PLANE
C        TAUM & TAUA - 2*PI/(PLANE SPACING) FOR MONO. AND ANAL.
C        R- REACTOR. M - MONOCHROMATOR. S- SAMPLE. A - ANALYSER.
C        D - DETECTOR.
C        ** FOR CONVERSION NOTE 4.1377MEV =1THZ **
C***********************************************************************
        IMPLICIT DOUBLE PRECISION (A-H,O-Z)        
        PARAMETER (PI=3.141592653589793239D0)
        COMMON AA(5),SM,SS,SA,AAA(6)
C
        CHARACTER*5 CUNIT
C
        COMMON /UNIT/CUNIT
        COMMON /EXTRA/XXC,YYC,FF,W,WI,WF
        COMMON/TRSG/SS1,SS2,SXX,SYY,GG1,GG2,GXX,GYY
        DIMENSION A(6),B(4),C(6),D(7),E(4),F(4)
C
C
C        C2=.0002908882
C        C1=6.2831853
C
        C1=2.*PI
        C2= C1/(360.D0*60.D0)
        A(1)=XXC
        A(2)=YYC
        A(3)=WI
        A(4)=WF
        DO 120 I=1,2
        C(I)=AAA(2+I)*C2
  120        C(I+3)=AAA(7-I)*C2
        C(3)=AA(3)*C2
        C(6)=AA(4)*C2
        
        IF ((C(1)*C(2).LE.0).OR.(C(5)*C(6).LE.0)) GOTO 99
C        IF ((C(1)**2+C(2)**2.LE.0).OR.(C(4)**2+C(5)**2.LE.0)) GOTO 99
        
        
        
        
        T1=C1/AA(1)
        T2=C1/AA(2)
        IF(ABS(A(1))-.000001)163,163,170
163        A(1)=1E-06
170        P1=-SS*SM
        P2=SA*SM
        A(5)=0.695*SQRT(A(3))
        A(6)=0.695*SQRT(A(4))
        A7=SQRT(A(1)**2+A(2)**2)
        A9=ATAN(ABS(A(2)/A(1)))
        A9=SIGN(1.D0,A(1)*A(2))*A9+0.5*PI*(1.-SIGN(1.D0,A(1)))
        C1=(A(5)**2+A7**2-A(6)**2)/(2*A7*A(5))
        B1=ATAN(SQRT(1.-C1**2)/C1)
        IF(B1)300,310,310
  300        B1=B1+PI
  310        C1=(A(6)**2+A7**2-A(5)**2)/(2*A(6)*A7)
        B2=ATAN(SQRT(1.-C1**2)/C1)
        IF(B2)350,360,360
  350        B2=B2+PI
  360        S1=T1/(2*A(5))
        B(1)=ATAN(S1/SQRT(1.-S1**2))
        S1=T2/(2.*A(6))
        B(3)=ATAN(S1/SQRT(1.-S1**2))
        DO 630 J1=1,2
        J2=(J1-1)*3
        I1=1+J2
        I2=2+J2
        I3=3+J2
        I4=4+J2
        D(I1)=C(I1)*C(I2)/SQRT(C(I1)**2+C(I2)**2)
        C1=SQRT(4./(C(I1)**2+C(I2)**2)+1/C(I3)**2)
        D(I2)=(C(I2)**2-C(I1)**2)/(C(I2)**2+C(I1)**2)/C1
        D(I3)=2.*C(I2)**2/(C(I1)**2+C(I2)**2)/C1
        IF(C(I2)-C(I1))510,530,510
  510        D(I4)=2.*C(I2)**2/(C(I2)**2-C(I1)**2)
        GOTO 540
  530        D(I4)=1.E37
  540        J3=1+2*(J1-1)
        E(J3)=2.*A(2+J1)*D(I1)/TAN(B(J3))
        E(J3+1)=2.*A(2+J1)*D(I2)/TAN(B(J3))
        F(J3)=A(4+J1)*D(I1)/SIN(B(J3))
        B(J3+1)=ATAN(D(I4)*TAN(B(J3)))
        IF(B(J3+1))610,620,620
  610        B(J3+1)=B(J3+1)+3.14159265
  620        F(J3+1)=A(4+J1)*D(I3)/SIN(B(J3+1))
  630         CONTINUE
        S1=GG1*FF/.4826
        S2=GG2*FF/.4826
        S1P=S1+(1+SM)*(S2*COS(A9)-S1*SIN(A9))*SIN(A9)
        S2=S2+(1+SM)*(S1*SIN(A9)-S2*COS(A9))*COS(A9)
        S1=S1P
        X=A9+B1*P1+B(1)
        W1=-S1*F(1)*COS(X)-S2*F(1)*SIN(X)+E(1)
        X=A9+B1*P1+B(2)
        W2=-S1*F(2)*COS(X)-S2*F(2)*SIN(X)+E(2)
        X=A9+(3.14159265-B2*P1)-B(3)*P2
        W3=-S1*F(3)*COS(X)-S2*F(3)*SIN(X)+E(3)
        X=A9+(3.14159265-B2*P1)-B(4)*P2
        W4=-S1*F(4)*COS(X)-S2*F(4)*SIN(X)+E(4)
        FACT=W*.4826/FF
        W1=SQRT(W1*W1+W2*W2)*FACT
        W2=SQRT(W3*W3+W4*W4)*FACT
        WW=SQRT(W1*W1+W2*W2)
        WRITE(IOCONS+1,101)CUNIT
  101        FORMAT(' Phonon Width Contributions [A-1]',A)
        WRITE(IOCONS+1,100)W1,W2,WW
  100        FORMAT(' W( Mono.)= ',F10.5/
     1        ' W( Anal.)= ',F10.5/' W( Total)= ',F10.5)
        RETURN

   99   continue
  102        FORMAT(' ALF1..4 must be > 0 !')
        write(IOCONS+1,102)
        RETURN
        END
C
C-----------------------------------------------------------        
      SUBROUTINE GETROANAL(RO)
C return "optimal" monochromator and analyzer curvatures 
C calculated analytically
C *** J.S. 3/6/1997     
C-----------------------------------------------------------      

      IMPLICIT NONE
      REAL*4 RO(4)
      INTEGER*4 NGRADM,NGRADA
      REAL*8 TETAM,HIMON,ROHM,ROVM,ETM,PM,PM0,RKINM,RLAMM
      REAL*8 TETAA,HIANA,ROHA,ROVA,ETA,PA,PA0,RKINA,RLAMA
      REAL*8 VL0,VL1,VL2,VL3,VLMIN,VLMAX
      REAL*8 PI,TDR,TMR,R8LN2,HSQOVM,ZNORM
      REAL*8 VKI,SLAMDI,EI0,VKF,SLAMDF,EF0,DUM(17)
      REAL*8 CRYD(2)
      COMMON /MONOCH/TETAM,HIMON,ROHM,ROVM,ETM,NGRADM,PM,PM0,RKINM,RLAMM
      COMMON /ANALYS/TETAA,HIANA,ROHA,ROVA,ETA,NGRADA,PA,PA0,RKINA,RLAMA
      COMMON /SETP/VL0,VL1,VL2,VL3,VLMIN,VLMAX
      COMMON /CTANTS/PI,TDR,TMR,R8LN2,HSQOVM,ZNORM
      COMMON /MEAN/VKI,SLAMDI,EI0,VKF,SLAMDF,EF0,DUM
      COMMON /CRYD/CRYD
      REAL*8 THM,CHIM,THA,CHIA
      
      
                
      THM=PI/SQRT(EI0/HSQOVM*2.)/CRYD(1)
      THM=ABS(ASIN(THM))
      CHIM=HIMON*TDR      
      THA=PI/SQRT(EF0/HSQOVM*2.)/CRYD(2)
      THA=ABS(ASIN(THA))
      CHIA=HIANA*TDR
c      RO(1)=SIN(THM+CHIM)/2./VL1*100
c change to monochromatic focusing:
      RO(1)=SIN(THM+CHIM)/VL1*100
c      write(*,*) 'VL1, THM, CHIM: ',VL1,THM*180/PI,CHIM*180/PI
      RO(2)=1./VL1/(2.*SIN(THM)*COS(CHIM))*100
c      RO(3)=(VL2*SIN(THA+CHIA) + VL3*SIN(THA-CHIA))/2./VL2/VL3*100
      RO(3)=SIN(THA-CHIA)/VL2*100
      RO(4)=(1./VL2+1./VL3)/(2.*SIN(THA)*COS(CHIA))*100
      END

          
C-----------------------------------------------------------        
      SUBROUTINE GETRO(ICONS)
C generates "optimal" monochromator and analyzer curvatures      
C-----------------------------------------------------------                    
      IMPLICIT NONE
      REAL*4 RO(4)
      REAL*8 DUM(42),ROMH,ROMV,ROAH,ROAV,SDI,SHI  
      COMMON DUM,ROMH,ROMV,ROAH,ROAV,SDI,SHI       
      INTEGER*4 NOS,I,ICONS
      REAL*8 RET(40)
      COMMON/CHRT/NOS,RET     
      CHARACTER*10 REMARK(4)
            
1     FORMAT(' ROMH = ',F8.4,A20)
2     FORMAT(' ROMV = ',F8.4,A20)
3     FORMAT(' ROAH = ',F8.4,A20)
4     FORMAT(' ROAV = ',F8.4,A20)
                 
      CALL GETROANAL(RO)
      DO 10 I=1,4
10         REMARK(I)=' '
                  
      IF((NOS.EQ.0).OR.((NOS.GE.1).AND.(RET(1).EQ.1.))) THEN
         ROMH=RO(1)
         REMARK(1)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.2).AND.(RET(2).EQ.1.))) THEN
         ROMV=RO(2)
         REMARK(2)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.3).AND.(RET(3).EQ.1.))) THEN            
         ROAH=RO(3)
         REMARK(3)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.4).AND.(RET(4).EQ.1.))) THEN        
         ROAV=RO(4)
         REMARK(4)=' changed'             
      ENDIF
      IF(ICONS.GT.0) THEN      
        WRITE(ICONS+1,1) RO(1),' [m-1]  '//REMARK(1)   
        WRITE(ICONS+1,2) RO(2),' [m-1]  '//REMARK(2)        
        WRITE(ICONS+1,3) RO(3),' [m-1]  '//REMARK(3)        
        WRITE(ICONS+1,4) RO(4),' [m-1]  '//REMARK(4)                    
      ENDIF
      
      RETURN
      END


C-----------------------------------------------------------        
      SUBROUTINE GETROOPT(ICONS)
C generates "optimal" monochromator and analyzer curvatures 
C *** changed to optimize Vanad resolution and intensity
C *** J.S. 3/6/1997     
C-----------------------------------------------------------      
              
      IMPLICIT NONE
      REAL*4 RO(4),TOL ,DRO(4)
      REAL*8 DUM(42),ROMH,ROMV,ROAH,ROAV,SDI,SHI  
      COMMON DUM,ROMH,ROMV,ROAH,ROAV,SDI,SHI       
      INTEGER*4 NOS,IERR,I,ICONS
      REAL*8 RET(40)
      COMMON/CHRT/NOS,RET     
      CHARACTER*10 REMARK(4)
      COMMON /ERROR/IERR
      INTEGER*4 OPTPAR,OPTMERIT
      REAL*8 OPTEV
      COMMON /MCOPTIM/ OPTPAR,OPTMERIT,OPTEV
      REAL*4 OPTV 
      EXTERNAL OPTV 
      INTEGER*4 IUEXT,IUOUT
      COMMON/BATCH/ IUEXT,IUOUT
           
1     FORMAT(' ROMH = ',F8.4,A20)
2     FORMAT(' ROMV = ',F8.4,A20)
3     FORMAT(' ROAH = ',F8.4,A20)
4     FORMAT(' ROAV = ',F8.4,A20)
5     FORMAT(' Numerical optimization failed ! ') 
      
7     FORMAT('Optimization for Vanad peak:'/
     *  ' (1) Intensity (2) Resolution (3) I/R   (4) I/R^2 '/
     *  'Select figure of merit: ',$)


20    WRITE(*,7)
      READ(IUEXT,*) I
      WRITE(*,*)
      IF(I.LT.1.OR.I.GT.4) GOTO 20     
      OPTMERIT=I
      
                
      CALL GETROANAL(RO)
      TOL=0.01
      DO I=1,4
        DRO(I)=0.01 ! minimum increment for curvature = 0.01m^-1
      ENDDO  
      
      CALL LMOPT(OPTV,RO,4,TOL,DRO)
      
      IF (IERR.NE.0) WRITE(ICONS+1,5)
      
      DO 10 I=1,4
10         REMARK(I)=' '
                  
      IF((NOS.EQ.0).OR.((NOS.GE.1).AND.(RET(1).EQ.1.))) THEN
         ROMH=RO(1)
         REMARK(1)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.2).AND.(RET(2).EQ.1.))) THEN
         ROMV=RO(2)
         REMARK(2)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.3).AND.(RET(3).EQ.1.))) THEN          
         ROAH=RO(3)
         REMARK(3)=' changed'
      ENDIF
      IF((NOS.EQ.0).OR.((NOS.GE.4).AND.(RET(4).EQ.1.))) THEN        
         ROAV=RO(4)
         REMARK(4)=' changed'             
      ENDIF
      
      WRITE(ICONS+1,1) RO(1),' [m-1]  '//REMARK(1)   
      WRITE(ICONS+1,2) RO(2),' [m-1]  '//REMARK(2)   
      WRITE(ICONS+1,3) RO(3),' [m-1]  '//REMARK(3)   
      WRITE(ICONS+1,4) RO(4),' [m-1]  '//REMARK(4)   
      
      RETURN
      END
      
C

C-----------------------------------------------------------        
      SUBROUTINE GETROOPTMC(ICONS)
C Optimize curvature with M.C. simulation
C Only one of the curvatures can be optimized  
C *** J.S. 5/7/2001     
C-----------------------------------------------------------      
              
      IMPLICIT NONE
      REAL*4 RO(4),TOL 
      REAL*8 DUM(42),ROMH,ROMV,ROAH,ROAV,SDI,SHI  
      COMMON DUM,ROMH,ROMV,ROAH,ROAV,SDI,SHI       
      INTEGER*4 NOS,IERR,I,ICONS
      REAL*8 RET(40)
      COMMON/CHRT/NOS,RET    
      COMMON /ERROR/IERR
      CHARACTER*10 REMARK(4)
      INTEGER*4 OPTPAR,OPTMERIT
      REAL*8 OPTEV
      COMMON /MCOPTIM/ OPTPAR,OPTMERIT,OPTEV
      EXTERNAL OPTMC
      REAL*4 PAR(1),DPAR(1) 
      INTEGER*4 IUEXT,IUOUT
      COMMON/BATCH/ IUEXT,IUOUT
           
1     FORMAT(' ROMH = ',F8.4,A20)
2     FORMAT(' ROMV = ',F8.4,A20)
3     FORMAT(' ROAH = ',F8.4,A20)
4     FORMAT(' ROAV = ',F8.4,A20)
5     FORMAT(' Numerical optimization failed ! ') 
6     FORMAT(' Wrong syntax. Type> MRO n [e] ',/,
     *  'n=1 to 4 for ROMH, ROMV, ROAH, ROAV',/,
     *  'e   .. number of events in 1000 (default e=1)') 
                 
7     FORMAT(' (1) Incident flux',/,
     *       ' (2) flux/dE ',/,
     *       ' (3) flux/dE^2 ',/,
     *       ' (4) Powder peak',/,
     *       ' (5) Vanad peak ',/,
     * 'Select figure of merit: ',$)


20    WRITE(*,7)
      READ(IUEXT,*) I
      IF(I.LT.1.OR.I.GT.4) GOTO 20     
      OPTMERIT=I
c      WRITE(*,*) OPTMERIT
      
      CALL NESS_CONV(1)
      
      CALL GETROANAL(RO)  ! Analytical estimation    
      OPTPAR=NINT(RET(1))
      IF (NOS.GT.1) THEN
        OPTEV=RET(2)
      ELSE
        OPTEV=1.0        
      ENDIF        
      IF(OPTPAR.LT.0.OR.OPTPAR.GT.4.OR.NOS.LT.1) THEN
         WRITE(ICONS+1,6)
         RETURN
      ENDIF
      
      TOL=0.1
      DPAR(1)=0.05 ! minimum increment for vert. curvature = 0.05m^-1
      IF (OPTPAR.EQ.1.OR.OPTPAR.EQ.3) DPAR(1)=0.01 ! 0.01m^-1 for hor. curv.
      PAR(1)=RO(OPTPAR)
      CALL LMOPT(OPTMC,PAR,1,TOL,DPAR)
       
      IF (IERR.NE.0) THEN
         WRITE(ICONS+1,5)
         RETURN
      ENDIF         
      
      DO 10 I=1,4
10         REMARK(I)=' '
                  
      IF(OPTPAR.EQ.1) THEN
         ROMH=PAR(1)
         REMARK(1)=' changed'
         WRITE(ICONS+1,1) ROMH,' [m-1]  '//REMARK(1)   
      ENDIF
      IF(OPTPAR.EQ.2) THEN
         ROMV=PAR(1)
         REMARK(2)=' changed'
         WRITE(ICONS+1,2) ROMV,' [m-1]  '//REMARK(2)   
      ENDIF
      IF(OPTPAR.EQ.3) THEN
         ROAH=PAR(1)
         REMARK(3)=' changed'
         WRITE(ICONS+1,3) ROAH,' [m-1]  '//REMARK(3)   
      ENDIF
      IF(OPTPAR.EQ.4) THEN
         ROAV=PAR(1)
         REMARK(4)=' changed'             
         WRITE(ICONS+1,4) ROAV,' [m-1]  '//REMARK(4)   
      ENDIF
                      
      
      RETURN
      END
      
C
C
C
C***********************************************************************
      SUBROUTINE SAM_FLUX(ICOM)
C /// simulate flux at the sample (arg<>2) or at the detector (arg=2)
C /// by forward method (ICOM=1) or "from the sample" (ICOM=0)
C /// ICOM=2- monitor at position AGM(1)      
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'ness_common.inc'
        
      COMMON/ARG/AGM(40)              
      COMMON/CHRT/NOS,RET(40)
      
c      write(*,*) 'FLUX: ',ICOM
C// monitor 
      IF (ICOM.EQ.2) THEN
        CALL NESS_CONV(1)
        IMONIT=NINT(AGM(1))
        IF (IMONIT.LE.7) THEN 
           CALL NESS(7,0.D0)
        ELSE
           CALL NESS(6,0.D0)
        ENDIF   
        IMONIT=-1
        RETURN
      ENDIF
      
C// no monitor
      IMONIT=-1
      CALL NESS_CONV(1)
      IF (ICOM.EQ.3) THEN
         CALL NESS(8,0.D0)
      ELSE IF (ICOM.EQ.4) THEN
         CALL NESS(9,0.D0)
      ELSE IF (ICOM.EQ.5) THEN
         CALL NESS(10,0.D0)
      ELSE          
         IF (AGM(1).EQ.2) THEN      ! powder
           CALL NESS(4+ICOM,0.D0)
         ELSE IF (AGM(1).EQ.3) THEN !  TAS
           CALL NESS(1,0.D0) 
         ELSE IF (AGM(1).EQ.4) THEN !  TAS forward
           FTAS=1
           CALL NESS(1,0.D0)
           FTAS=0 
         ELSE IF (AGM(1).EQ.11) THEN !  double cryst.
           CALL NESS(11,0.D0) 
         ELSE                       ! flux at the sample
           CALL NESS(2+ICOM,0.D0)
         ENDIF
      ENDIF       
      RETURN
      END                   
        
C***********************************************************************
      SUBROUTINE SCAN_CHI(ICOM)
C /// simulate flux at the sample (arg<>2) or at the detector (arg=2)
C /// by forward method (ICOM=1) or "from the sample" (ICOM=0)
C /// ICOM=2- monitor at position AGM(1)      
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'ness_common.inc'
        
      INTEGER*4 IUEXT,IUOUT
      COMMON/BATCH/ IUEXT,IUOUT
      COMMON/ARG/AGM(40)              
      COMMON/CHRT/NOS,RET(40)
      COMMON/MONOCH/TETAM,HIMON,ROHM,ROVM,ETM,NGRADM,PM,PM0,RKINM,RLAMM
1     FORMAT('CHI = ',G12.6)      

      IF (NOS.GE.2) THEN
        CALL NESS_CONV(1)
        CHI0=HIMON
        DCHI=AGM(1)
        NCHI=NINT(AGM(2))
        EV=10.D0
        IF(NOS.GE.3) EV=AGM(3) ! number of events
        IF (NCHI.GT.100) NCHI=100
        IF (NCHI.LT.1) NCHI=1
        DO I=1,NCHI
          HIMON=CHI0+(I-(NCHI+1)/2)*DCHI
          CALL NESS_CONV(0)         
          WRITE(IUOUT,1) MON.CHI*180/PI
          CALL NESS(9,-ABS(EV))
        ENDDO
        HIMON=CHI0
        CALL NESS_CONV(0)         
      ENDIF  
      END

C***********************************************************************
      SUBROUTINE SCAN_TAS
C /// simulate standard TAS scan (DH,DK,DL,DE)
C//// using scattering cross-section defined by SQE_AMAG funciton
C accepts 4 arguments:
c a1 ... number of steps (obligatory) 
c a2 ... number of events (x1000) , default=10 
c a3 ... time (~monitor counts), default=100 
c a4 ... background (in cnts), default=0
C***********************************************************************
C
      IMPLICIT NONE

      INCLUDE 'ness_common.inc'
        
      INTEGER*4 I,J,NSTP
      REAL*8 EV,MONEF
      PARAMETER(MONEF=1D-8)
      
      REAL*8 DUM1(30),QHKL(4),STEP(4),DUM2(10) 
      COMMON DUM1,QHKL,STEP,DUM2
      INTEGER*4 IUEXT,IUOUT,NOS
      COMMON/BATCH/ IUEXT,IUOUT
      REAL*8 AGM(40),RET(40)
      COMMON/ARG/AGM              
      COMMON/CHRT/NOS,RET
      
      REAL*8 CNTS(128),CNTD(128),CNTE(128),KI(128),TIME,BCG
      REAL*8 QHKL0(4),K0     
      REAL*4 GASDEV
            
1     FORMAT(I3,2x,4(G10.4,1x),2(G12.4,2x))      
2     FORMAT('PNT   QH   QK   QL   EN   CNTS  MON')      
3     FORMAT('PNT   QH   QK   QL   EN   CNTS  TIME')      

      IF (NOS.GE.1) THEN
        K0=STP.KI
        CALL NESS_CONV(1)
        DO J=1,4
          QHKL0(J)=QHKL(J)
        ENDDO
c        WRITE(IUOUT,1) (QHKL0(J),J=1,4)
c        pause
        NSTP=NINT(AGM(1))
        EV=10.D0
        IF(NOS.GE.2) EV=AGM(2) ! number of events
        TIME=100.D0
        IF(NOS.GE.3) TIME=AGM(3) ! time
        BCG=0.D0
        IF(NOS.GE.4) BCG=AGM(4) ! background
        IF (NSTP.GT.100) NSTP=100
        IF (NSTP.LT.1) NSTP=1
        WRITE(IUOUT,2)
        DO I=1,NSTP
          DO J=1,4
            QHKL(J)=QHKL0(j)+(I-(NSTP+1)/2)*STEP(J)
          ENDDO            
          CALL NESS_CONV(1)         
          CALL NESS(8,ABS(EV))
          CNTS(I)=IINC
          CNTD(I)=TIME*I3AX
          CNTE(I)=TIME*DI3AX
          KI(I)=STP.KI
          WRITE(IUOUT,1) I,(QHKL(J),J=1,4),CNTD(I),CNTS(I) 
        ENDDO
        DO J=1,4
          QHKL(J)=QHKL0(J)
        ENDDO        
        CALL NESS_CONV(1)         
        DO I=1,NSTP
          SPCX(I)=QHKL0(4)+(I-(NSTP+1)/2)*STEP(4)
          SPCY(I)=CNTD(I)/CNTS(I)/MONEF ! normalize to monitor counts
          SPCD(I)=CNTE(I)/CNTS(I)/MONEF
        ENDDO
        WRITE(IUOUT,3)
        DO I=1,NSTP
          IF (BCG.GT.0) THEN ! add const. background and errors
            SPCY(I)=SPCY(I)+BCG
            SPCD(I)=SQRT(ABS(SPCY(I))+SPCD(I)**2)
            SPCY(I)=SPCY(I)+SQRT(ABS(SPCY(I)))*GASDEV()
          ENDIF  
          
          WRITE(IUOUT,1) I,(QHKL0(j)+(I-(NSTP+1)/2)*STEP(J),J=1,4),
     *                   SPCY(I),TIME/MONEF/CNTS(I)*KI(I)
        ENDDO
        SPCN=NSTP
      ENDIF  
      END

C***********************************************************************
      SUBROUTINE SCAN_THETA
C /// simulate standard TAS scan (A1,A2,A3,A4,A5,A6)
C accepts 4 arguments:
c a1 ... number of steps (obligatory) 
c a2 ... number of events (x1000) , default=10 
c a3 ... time (~monitor counts), default=100 
c a4 ... background (in cnts), default=0
C***********************************************************************
C
      IMPLICIT NONE

      INCLUDE 'ness_common.inc'
        
      INTEGER*4 MF
      PARAMETER (MF=65)
      INTEGER*4 NSTP
      REAL*8 EV,MONEF
      PARAMETER(MONEF=1D-8)
      
      INTEGER*4 IUEXT,IUOUT,NOS,I_IO
      COMMON/BATCH/ IUEXT,IUOUT
      REAL*8 AGM(40),RET(40)
      COMMON/ARG/AGM              
      COMMON/CHRT/NOS,RET
      CHARACTER*128 LINE
      character*50 filename
      REAL*8 AR(6)
      INTEGER*4 NA,IA(6),IWHAT,K,I,J
      REAL*4 FX(MF),FY(MF),DFY(MF),FY1(MF),DFY1(MF)
      REAL*8 CNTS(128),CNTD(128),CNTE(128),TIME
            
4     FORMAT('A',I1,'      ',$)      
5     FORMAT(G10.4,1x,$)
55    FORMAT(2(G10.4,1x))      
6     FORMAT('Axes [1..6]: ',$)
7     FORMAT(a)      
8     FORMAT('Steps [min]: ',$)
9     FORMAT('(1) Sample, (2) Powder, (3) Double-Crystal, (4) TAS :',$)
44     FORMAT(1x,4(2x,E13.5))

      IF (NOS.LT.1) THEN
        WRITE(IUOUT,*) 'Use number of points as the 1st argument'
        RETURN
      ENDIF  
      WRITE(IUOUT,9) 
      READ(IUEXT,*) IWHAT
      IF (IWHAT.EQ.4) THEN
        CALL SCAN_TAS
        RETURN
      ELSE IF (IWHAT.LT.1.OR.IWHAT.GT.4) THEN
        WRITE(*,*) 'UNDEFINED TASK: ',IWHAT 
        RETURN  
      ENDIF
C// initialize
        CALL NESS_CONV(1)
        DO J=1,6
          DTHAX(J)=0.
          IA(J)=J
          AR(J)=0.
        ENDDO
C// interpret arguments
        NSTP=NINT(AGM(1))
        EV=10.D0
        IF(NOS.GE.2) EV=AGM(2) ! number of events
        TIME=1.D0
        IF(NOS.GE.3) TIME=AGM(3) ! time
        IF (NSTP.GT.101) NSTP=101
        IF (NSTP.LT.1) NSTP=1
C// read angular steps from input        
        WRITE(IUOUT,6) 
        READ(IUEXT,7) LINE   ! read axes indexes
        CALL GETLINPARG(LINE,AR,6,NA)
        DO I=1,NA
          IA(I)=INT(AR(I))
          IF(IA(I).GT.6.OR.IA(I).LE.0) IA(I)=0
        ENDDO  
        WRITE(IUOUT,8) 
        READ(IUEXT,7) LINE  ! read axes steps
        
        CALL GETLINPARG(LINE,AR,6,J) 
        IF (J.NE.NA) THEN
          WRITE(*,*) 'EACH AXIS MUST HAVE A STEP DEFINED !!'
        ENDIF
          DO I=1,NA
             WRITE(IUOUT,4) IA(I) 
          ENDDO
          WRITE(IUOUT,*)
          DO I=1,NA
             WRITE(IUOUT,5) AR(I) 
          ENDDO
          WRITE(IUOUT,*)
C// get only valid steps
        K=0
        DO I=1,NA
          IF (IA(I).NE.0.AND.AR(I).NE.0) THEN
             K=K+1
             AR(K)=AR(I)
             IA(K)=IA(I)
             WRITE(IUOUT,4) IA(I)  ! write header
          ENDIF               
        ENDDO
        NA=K
        WRITE(IUOUT,*) 'CNTS    ERR' 
        DO J=1,MF
          FY1(J)=0
          DFY1(J)=0            
        ENDDO

C// Start scan
        DO I=1,NSTP
          DO J=1,NA
            DTHAX(IA(J))=(I-(NSTP+1)/2)*AR(J) 
          ENDDO
          DO J=1,NA
            WRITE(IUOUT,5) DTHAX(IA(J)) 
          ENDDO
          CALL NESS_CONV(0) 
          IF (IWHAT.EQ.1) THEN        
              CALL NESS(2,ABS(EV)) 
              CNTD(I)=TIME*IINC
              CNTE(I)=TIME*DIINC
          ELSE IF (IWHAT.EQ.2) THEN 
              CALL NESS(4,ABS(EV)) 
              CNTD(I)=TIME*IPWD
              CNTE(I)=TIME*DIPWD
          ELSE IF (IWHAT.EQ.3) THEN 
              CALL NESS(11,ABS(EV)) 
              CNTD(I)=TIME*I3AX
              CNTE(I)=TIME*DI3AX
          ENDIF              
          CNTS(I)=IINC
          WRITE(IUOUT,55) CNTD(I),CNTE(I)
          CALL PSD_ARRAY(FX,FY,DFY,MF)
          DO J=1,MF
            FY1(J)=FY1(J)+FY(J)
            DFY1(J)=DFY1(J)+DFY(J)**2            
          ENDDO
        ENDDO
C// End scan, reset configuration
        DO J=1,MF
            DFY1(J)=SQRT(DFY1(J))            
        ENDDO
        DO J=1,6
          DTHAX(J)=0.
        ENDDO
        CALL NESS_CONV(0)   
C// save integrated profile at the PSD        
      I_IO=22
      filename=' '
12    format(a50)
13    FORMAT(' PSD data output: ',$)
      WRITE(IUOUT,13)
      read(IUEXT,12) filename

      IF(filename(1:1).EQ.' '.OR.filename(1:1).EQ.CHAR(0)) then  ! generate automatic filename
        GOTO 200
      ELSE
        Open(Unit=i_IO,File=filename,err=999,Status='Unknown')
        write(i_io,*) 'X      INT       ERR    '
        do i=1,MF
          write(i_IO,44) FX(I),FY1(I),DFY1(I) 
        enddo  
        close(i_io)
      ENDIF
             
C// Fill arrays with results
200        DO I=1,NSTP
          SPCX(I)=(I-(NSTP+1)/2)*AR(1)
          SPCY(I)=CNTD(I)
          SPCD(I)=CNTE(I)
        ENDDO
C// List data
        DO J=1,NA
          WRITE(IUOUT,4) J 
        ENDDO
        WRITE(IUOUT,*) 'CNTS    ERR'         
        DO I=1,NSTP
          DO J=1,NA
            WRITE(IUOUT,5) (I-(NSTP+1)/2)*AR(J) 
          ENDDO
          WRITE(IUOUT,55) SPCY(I),SPCD(I)          
        ENDDO
        SPCN=NSTP
      return  
999   write(*,*) 'Cannot open file as unit ',i_IO      
      return
      END
C
C***********************************************************************
      SUBROUTINE BENCH
C /// simulate flux at the sample (arg<>2) or at the detector (arg=2)
C /// by forward method (ICOM=1) or "from the sample" (ICOM=0)      
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'ness_common.inc'
        
      COMMON/ARG/AGM(40)              
      COMMON/CHRT/NOS,RET(40)

      NB=1000
      IMONIT=-1
      IF(NOS.GT.0) NB=NINT(1000*AGM(1))
      CALL NESS_CONV(1)
      CALL NESS(2,0.D0)
      CALL NESS_BENCH(NB)
      
      RETURN
      END                   
C
C
C***********************************************************************
      SUBROUTINE ROCK(ICR)
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'ness_common.inc'
        
      COMMON/ARG/AGM(40)              
      COMMON/CHRT/NOS,RET(40)
      REAL*8 RTH(129),DTH,DIVH,DIVV

      NC=1000
      NTH=65
      DTH=2./180./60.*PI
      DIVH=0.
      DIVV=0.
      IF(NOS.GT.0) NC=NINT(1000*AGM(1))
      IF(NOS.GT.1) NTH=NINT(AGM(2))
      IF(NTH.LT.11) NTH=11
      IF(NTH.GT.129) NTH=129
      IF(NOS.GT.2) DTH=AGM(3)/180/60*PI
      IF(NOS.GT.3) DIVH=AGM(4)/180/60*PI
      IF(NOS.GT.4) DIVV=AGM(5)/180/60*PI
      CALL NESS_CONV(1)
      CALL SPEC_INI(0,3)
      if (NOS.GT.0.AND.AGM(1).EQ.0.) THEN
        CALL TEST_SYMMETRY(NC,NTH,DTH,DIVH,DIVV)      
      ELSE
        CALL NESS_ROCK(ICR,NC,NTH,DTH,RTH,DIVH,DIVV)
        OPEN(22,FILE='rcurve.dat',STATUS='unknown',ERR=100)
1       format(a)
2       format(2(E11.5,4x))
        write(22,1) 'theta[min]   r(theta)'
        DO I=1,NTH
          write(22,2) (-(NTH-1)/2.+I*1.)*DTH*180*60/PI,RTH(I)
        ENDDO      
100     close(22)
      ENDIF  
      
      RETURN
      END                   
C
C***********************************************************************
      SUBROUTINE TYPECFG
C /// print complete configuration of all components
C***********************************************************************
C
      IMPLICIT DOUBLE PRECISION (A-H,O-Z)

      INCLUDE 'ness_common.inc'
        
      COMMON/ARG/AGM(40)              
      COMMON/CHRT/NOS,RET(40)
      COMMON/BATCH/ IUEXT,IUOUT        

      CALL NESS_CONV(1)
      CALL SPEC_INI(0,3)
      IF(NOS.GT.0) THEN
         N=NINT(AGM(1))
         IF(N.EQ.1) CALL SLIT_WRITE(IUOUT,Source)
         IF(N.EQ.2) CALL BENDER_WRITE(IUOUT,gdea)
         IF(N.EQ.3) CALL BENDER_WRITE(IUOUT,guide)
         IF(N.EQ.4) CALL BENDER_WRITE(IUOUT,sol1)
         IF(N.EQ.5) CALL CRYST_WRITE(IUOUT,mon)
         IF(N.EQ.6) CALL BENDER_WRITE(IUOUT,sol2a)
         IF(N.EQ.7) CALL BENDER_WRITE(IUOUT,sol2)
         IF(N.EQ.8) CALL SLIT_WRITE(IUOUT,sam)
         IF(N.EQ.9) CALL BENDER_WRITE(IUOUT,sol3)
         IF(N.EQ.11) CALL SLIT_WRITE(IUOUT,det)
      ELSE
          CALL WRITE_SETUP(IUOUT,8)
      ENDIF
      
      RETURN
      END                   

C-----------------------------------------------------------        
      SUBROUTINE SET_DEVICE
C set graphics device string for PGPLOT
C *** J.S. 13/6/1998     
C-----------------------------------------------------------      
      CHARACTER*30 DEVSTR
      COMMON /GRAPH/ DEVSTR
      COMMON/BATCH/ IUEXT,IUOUT        
1     FORMAT(' Graphics device (/? for help) : ',$)
2     format(a)
      write(IUOUT,1)
      read(IUEXT,2) DEVSTR   
      return
      end
      
C-----------------------------------------------------------        
      SUBROUTINE SETVAR(IVAR)
C    
C-----------------------------------------------------------      
      IMPLICIT NONE
      INCLUDE 'const.inc'
      INTEGER*4 IVAR,NOS,IUEXT,IUOUT
      REAL*8 RET(40)
      COMMON/BATCH/ IUEXT,IUOUT        
      COMMON/CHRT/NOS,RET
1     FORMAT(' Source flux [1e14 n/s/cm^2] : ',F10.4)
2     FORMAT(' Source temperature [K] : ',F6.0)
      if (IVAR.EQ.1) then         
         IF(NOS.NE.0) SFLUX=RET(1)
         write(IUOUT,1) SFLUX
      endif   
      if (IVAR.EQ.2) then         
         IF(NOS.NE.0) STEMP=RET(1)
         write(IUOUT,2) STEMP
      endif   
      return
      end


C--------------------------------
      SUBROUTINE SETCFG(ICOM)
c--------------------------------
      CHARACTER*60 CFGNAME
      CHARACTER*128 CFGPATH
      COMMON/CONFNAME/ CFGNAME,CFGPATH

      COMMON/BATCH/ IUEXT,IUOUT        
1     FORMAT(' Configuration filename : ',$)
2     format(a)
      if(ICOM.EQ.0) THEN
         cfgname='simres00.cfg'
      else
        write(IUOUT,1)
        read(IUEXT,2) CFGNAME
        if(CFGNAME(1:1).eq.' '.or.cfgname(1:1).eq.char(0)) THEN
           cfgname='simres00.cfg'
        endif        
        write(IUOUT,*) CFGNAME(1:50)
      endif  
      return
      end
      
C-----------------------------------------------------------------------
      SUBROUTINE SETRESPATH(SARG)
C select search path for configuration files      
C-----------------------------------------------------------------------
      IMPLICIT NONE
      
      CHARACTER*(*) SARG
      CHARACTER*128 MYPATH
      INTEGER*4 IS,IL
      
      INTEGER*4 IUEXT,IUOUT
      COMMON/BATCH/ IUEXT,IUOUT        
      
      CHARACTER*60 CFGNAME
      CHARACTER*128 CFGPATH
      COMMON/CONFNAME/ CFGNAME,CFGPATH
      CHARACTER*1 PATHDEL
      
      DATA PATHDEL /'/'/      
      
1     FORMAT(' Additional search path for configuration files [',
     *         a,'] : ',$)
2     FORMAT(A)       
3     FORMAT(' Configurations will be searched in ',a)
       
      CALL BOUNDS(CFGPATH,IS,IL)
C Get pathname from dialog or from the argument SARG 
      IF (SARG.EQ.' ') THEN
        IF (IL.LE.0) THEN
          write(IUOUT,1) 'current folder'
          read(IUEXT,2) MYPATH
        ELSE
          write(IUOUT,1) CFGPATH(IS:IS+IL-1)
          read(IUEXT,2) MYPATH 
          IF (MYPATH(1:1).EQ.' '.OR.MYPATH(1:1).EQ.CHAR(0)) THEN
             MYPATH=CFGPATH(IS:IS+IL-1)
          ENDIF
        ENDIF  
      ELSE
        MYPATH=SARG
      ENDIF
C Interpret MYPATH, ensure that ending / is present       
      CALL BOUNDS(MYPATH,IS,IL)
      IF ((IL.LE.0).OR.
     *    (IL.EQ.1.AND.MYPATH(IS:IS+IL-1).EQ.'.').OR.
     *    (IL.EQ.2.AND.MYPATH(IS:IS+IL-1).EQ.'.'//PATHDEL)) THEN
         CFGPATH=' '
         write(IUOUT,3) 'current folder'
         RETURN
      ENDIF
      
c      write(*,*) '<'//MYPATH(IS:IS+IL-1)//'>'
      IF(MYPATH(IS+IL-1:IS+IL-1).NE.PATHDEL) THEN
         CFGPATH=MYPATH(IS:IS+IL-1)//PATHDEL
         IL=IL+1
      ELSE   
         CFGPATH=MYPATH(IS:IS+IL-1)
      ENDIF
      write(IUOUT,3) CFGPATH(1:IL)
      END  

CC--------------------------------
      SUBROUTINE DOSHELL(COMM)
c--------------------------------
      CHARACTER*72 COMM,COMM1
      COMMON/BATCH/ IUEXT,IUOUT        
1     FORMAT(' Command : ',$)
2     format(a)
      if((COMM.EQ.' ').OR.(COMM(2:2).EQ.CHAR(0))) THEN 
        write(IUOUT,1)
        read(IUEXT,2) COMM1
      else
        COMM1=COMM
      endif
      write(*,*) COMM1
      CALL SYSTEM(COMM1)
      return
      end
      

      BLOCK DATA
C***********************************************************************

      REAL*8 NAM(100)
      COMMON/CHNM/NAM
      CHARACTER*63 HLP(53),HLPOPT(18)      
      COMMON/HELP/ HLP,HLPOPT
      INCLUDE 'const.inc'
      REAL*8  SPINT
      COMMON /SPIN/ SPINT
        
      DATA NAM/
     1  'DM','DA','ETAM','ETAA','ETAS','SM','SS','SA','KFIX','FX',
     2  'ALF1','ALF2','ALF3','ALF4','BET1','BET2','BET3','BET4',
     3  'AS','BS','CS','AA','BB','CC','AX','AY','AZ','BX','BY','BZ',
     4  'QH','QK','QL','EN','DH','DK','DL','DE','GH','GK',
     5  'GL','GMOD','ROMH','ROMV','ROAH','ROAV','SDI','SHI','DTH',
     6  'BRAG','MPROF','AMOD','EMOD ','PLOT','PRINT','RO','SPOS',
     7  'WRITE','PWDS','FLUX','CFG','SETUP','MONO','THETA','FLIP','MAG'
     8  ,'SPIN','GRFDEV','FLX','TEMP','POLAR','CRYST','SHELL','NFLUX'   !25
     9  ,'BENCH','OSC','SVOL','LPROF','DPROF','MONIT','ROCK','NRO'
     *  ,'MRO','TAS','PWD','TPROF','SCHI','SRCX','SRCY','SCAN','AROCK'
     &  ,'ANAL',8*' '/

      DATA HLP/
Cxxxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccxxx
     1'LIST    - list of parameters                                  ',
     2'SAVE    - saves parameters in the current parameter file      ',
     3'SAVE as - saves parameters, asks for filename                 ',     
     4'EXIT    - exit with warning if data file not saved            ',
     5'EXFF    - exit without saving data                            ', 
     6'CFG     - asks for configuration file name                    ', 
     7'FILE    - load data file or RESCAL parameter file             ',
     8'CRYST 1 - set monochromator name, select from "crystal.lib"   ',
     9'CRYST 2 - set analyzer name, select from "crystal.lib"        ',
     *'                                                              ', 
     1'GRFDEV  - set output device nape for PGPLOT                   ',  
     2'PLOT    - plot the beam at the sample in given projection     ',  
     3'PRINT   - print the last plotted image                        ',  
     4'                                                              ', 
     5'RO      - sets "optimal" crystal curvatures                   ', 
     6'TEMP    - Temperature of the moderator [K]                    ',
     6'FLX     - Integral flux [n/s/cm2] at the channel nose         ',
     7'OSC     - switch the oscilation of collimators on/off         ',
     8'SPOS x y z  - changes sample position (for Monte Carlo only)  ',
     9'              (x // ki, z is vertical coordinate)             ',
     *'MONO a b c  set optional parameters to the monochromator      ',
     *'ANAL a b c  set optional parameters to the analyzer           ',
     1'      a - d-gradient in [0.001/cm]                            ',  
     2'      b - orientation of the gradient [deg],(=0 along z)      ',  
     3'      c - lamella thickness [um] to calculate primary ext.    ',  
     4'THETA a b - set rocking angles [min] for monoch. and analyzer ', 
     4'SETUP n - prints out the parameters of the n-th component     ', 
     4'                                                              ', 
     5'FLUX    - simulate flux at the sample, in forward direction   ', 
     6'FLUX  2 - flux at the detector, supposing ideal powder sample ',
     7'NFLUX   - as FLUX, but starts simulation at the sample        ',
     8'NFLUX 2 - as FLUX 2, but starts simulation at the sample      ',
     8'TAS     - R(Q,E), TAS resolution funciton                     ',
     8'SCAN np [nev mon  bcg]  .. TAS scan simulation                ',
     8'PWD     - R(Q), diffraction resolution funciton               ',
     8'PWDS    - scan powder peak by multidetector with step DTH     ',
     9'MONIT n - run simulation from source to monitor n             ', 
     9'SVOL    - plots the image of sample gauge volume              ', 
     *'LPROF   - plots the lambda distribution across the sample     ', 
     1'DPROF   - plots the beam profile along the detector           ', 
     1'TPROF   - plots the powder peak profile                       ', 
     1'MPROF   - plots multidetector scan                            ', 
     2'BENCH n - start NFLUX command, then run n*1000 events through ', 
     3'          a benchmark (returns performance of components)     ',
     4'ROCK  [nev np dth] - cal. rocking curve of the monochromator  ', 
     4'AROCK [nev np dth] - cal. rocking curve of the analyzer       ', 
     5'SHELL   - permits to perform a single shell command           ',
     6'SRCX [A B X0] - source inhomogeneity along x                  ',
     6'SRCY [A B Y0] - source inhomogeneity along y                  ',
     6'AMOD [0|1] - turn off|on analyzer flat-cone mode              ',
     6'EMOD [0|1] - turn off|on E=const. mode                        ',
     6'                                                              ',
     7'       (read simres.hlp and manuals for other details)       '/ 
  
      DATA HLPOPT/
Cxxxxxccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccxxx
     1'RESTRAX options:                                              ',
     2'                                                              ',
     2'-d        debug mode (no sampling optimization)               ',
     3'-sx       randomize with x, Example: -s10001                  ',
     4'-tx       test random generator with dim=x, e.g. -t5          ', 
     5'-flxfile  read flux distribution dPhi/dLambda from file       ',
     6'            Example: -flxdist.dat                             ',
     7'-flhx     set horizontal source divergence limit to x [deg]   ',
     8'            Example: -flh1.5                                  ',
     9'-flvx     set vertical source divergence limit to x [deg]     ',
     *'            Example: -flv2.3                                  ',
     1'-Voigt    use pseudo-Voigt mosaic distribution                ',
     2'-ran1     use Numerical Recipes RAN1 generator                ',
     3'-rand     use system random number generator (RAN)            ',
     3'-noopt    no automatic sampling optimization                  ',
     3'-nmon     weight by monitor efficiency ~ 1/ki                 ',
     3'-help     show this help                                      ',
     4'                                                              '/

      DATA SPINT /-1/      
      DATA PI,HOVM,HSQOV2M,GammaNi/3.1415926535897932,6.296E-1,2.0721,
     * 0.00173/
      DATA SFLUX,STEMP /1.0,310./ 
      DATA FLXH,FLXV,FLXN / 0, 0, 0/
      END
                       
c dPhi/dLambda is in [1e12/s/cm^2/Ang]      


