C////////////////////////////////////////////////////////////////////////
C $Id: simres_main.f,v 1.135 2019/06/12 17:55:13 saroun Exp $
C
C                 ************************************
C                 ***                              ***
C                 ***        S I M R E S           ***
C                 ***                              ***
C                 ***   (C) J.Saroun & J.Kulda     ***
C                 ***        ILL Grenoble          ***
C                 ***      ray-tracing version     ***
C                 ************************************
C
C A clone of RESTRAX: http://neutron.ujf.cas.cz/restrax
C
C Provides more flexible (and more realistic) ray-tracing code useful for simulation
C of newly designed or upgraded instruments and optimisation of their configuration.
C This version permits to simulate intensity and distribution of neutron beam in both
C real and momentum subspaces at different positions along the TAS beamline.
C Arrangements derived from TAS setup can also be simulated - they involve e.g. powder
C diffractometers equipped with multidetectors, neutron guides or different configurations
C of primary spectormeter (i.e. crystal monochromator with series of collimator or
C guide segments).
C
C****************************************************************************
C *** For all additional information contact the authors:                 ***
C ***                                                                     ***
C ***                 kulda@ill.fr        saroun@ujf.cas.cz               ***
C ***                                                                     ***
C****************************************************************************

C---------------------------------------------
      PROGRAM RESTRAX
c      write(*,*) "RESTRAX MAIN"
      CALL RESTRAX_MAIN
      END

C-----------------------------
      SUBROUTINE LOGO
C-----------------------------
      USE IO
      IMPLICIT NONE
      INCLUDE 'config.inc'

1     FORMAT(2x,'-----------------------------------------------',/,
     &       2x,'S I M R E S  - Monte Carlo ray-tracing ',/,
     &       2x,'Version: ',a40,/,
     &       2x,'Build:   ',a40,/,
     &       2x,'-----------------------------------------------',/,
     &       2x,'(C) J.Saroun & J.Kulda',/,
     &       2x,'ILL Grenoble, NPI Rez near Prague',/,
     &       2x,' ',/,
     &       2x,'-----------------------------------------------',/,
     &       2x,'type ? for command list',/,/)


      write(sout,1) PACKAGE_VERSION,PACKAGE_DATE
      end SUBROUTINE LOGO

C-----------------------------
      SUBROUTINE VERSIONINFO
C-----------------------------
      USE XMLINFO
      IMPLICIT NONE
      INCLUDE 'config.inc'
      if (XMLOUT.gt.0) then
        call XML_RSXDUMP(sout,'VERSIONINFO',1)
        call XML_TAG(sout,'VERSION',trim(PACKAGE_VERSION),2)
        call XML_TAG(sout,'BUILD',trim(PACKAGE_DATE),2)
        call XML_RSXDUMP(sout,'VERSIONINFO',0)
      endif
      end SUBROUTINE VERSIONINFO

C-------------------------------------------------------------
      SUBROUTINE SETXML(ARG)
C Read instrument configuration in XML format from console.
C Only for GUI.
c-------------------------------------------------------------
      use XMLINFO
      USE XMLSIMRES
      use XMLWRITER
      use MESSAGES
1     format('SETXML [',a,']')
      CHARACTER*(*) ARG
      integer :: ierr=0
      ! write(*,1) trim(ARG)
      if (XMLOUT.gt.0) then
      ! no argument: load from console and then send back updated cfg.
        if (len_trim(ARG).le.0) then
          call load_config('STDIN',ierr)
          if (ierr==0) call save_config('STDIN')
          call FLUSH_MESSAGES
      ! UPDATE: send config. to GUI
        else if(trim(ARG).eq.'UPDATE') then
          call save_config('STDIN')
      ! other, assumed component name, send it to GUI
        else
          call save_component('STDIN',trim(ARG))
        endif
      else
        call MSG_WARN('Can''t update configuration from console: not in XML mode',1)
      ENDIF
      END SUBROUTINE SETXML

C--------------------------------
      SUBROUTINE DOSHELL(COMM)
c execute shell command
c--------------------------------
!      USE CONSTANTS
      use IO
      IMPLICIT NONE
      CHARACTER*(*) COMM
      CHARACTER*256 LINE
      INTEGER*4 IS,L

1     FORMAT(' Command : ',$)
2     format(A)
      LINE=' '
      if((COMM(1:1).EQ.' ').OR.(COMM(2:2).EQ.CHAR(0))) THEN
        write(SOUT,1)
        call RSXREAD(LINE)
      else
        L=MIN(256,LEN(COMM))
        LINE=COMM(1:L)
      endif
      CALL BOUNDS(LINE,IS,L)

      IF (L.GT.0) THEN
c       write(*,*) 'calling shell: ',LINE(IS:IS+L-1)
        CALL SYSTEM(LINE(IS:IS+L-1))
      ENDIF
      end

C-----------------------------------------------------
      SUBROUTINE READINIFILE(JOBNAME)
C   read initialization file
C CFGNAME = configuration file
C DATAPATH = path to the data files
C OPENFILE = data or RESCAL file to open
C return JOBNAME .. filename of a job file to be executed at the startup
c-----------------------------------------------------
      USE CONSTANTS
      use FILETOOLS
      use IO
      use SOURCES_TABLE
      IMPLICIT NONE
      CHARACTER*128 LINE,S
      character(MAX_FNAME_LENGTH) :: FRES,FRESPATH
      CHARACTER*(*) JOBNAME
      INTEGER*4 IRES,IERR,IU

1     FORMAT(A)
      JOBNAME=' '
      S=' '
      call OPENRESFILE('restrax.ini','ini',0,FRES,FRESPATH,IU,IERR)
!      CALL OPENRESFILE('restrax.ini',22,IRES,0)
      IF(IERR.EQ.0) THEN
         IRES=0
         DO WHILE(IRES.EQ.0)
           READ(IU,1,END=100,iostat=IRES) LINE
           IF(LINE(1:1).NE.'#') THEN
!           WRITE(*,*) 'INI: ',trim(LINE)
!             read(*,*)
             CALL READ_STR('CFGNAME',LINE,CFGNAME,IERR)
             CALL READ_STR('XMLCONFIG',LINE,XCONFIG_NAME,IERR)
          !   CALL READ_STR('DATAPATH',LINE,DATPATH,IERR)
             CALL READ_STR('JOB',LINE,JOBNAME,IERR)
             !CALL READ_STR('OPENFILE',LINE,RESCAL_NAME,IERR)
          !   CALL READ_STR('FLUX',LINE,S,IERR)
          !   IF (IERR.EQ.0) CALL READ_FLUX(S)
             !CALL READ_STR('PULSE',LINE,S,IERR)
             !IF (IERR.EQ.0) CALL READ_PULSE(S)
c             WRITE(*,*) 'READ: |',trim(S),'|',IERR
c             read(*,*)
           ENDIF
         ENDDO
100      CLOSE(IU)
      ENDIF
      END


C-----------------------------------------------------------------------------
      SUBROUTINE PROCARG
C// Process command line arguments for SIMRES
C-----------------------------------------------------------------------------
      USE COMPONENTS
      USE TRACINGDATA
      USE CONSTANTS
      use FILETOOLS
      use IO
      USE GENERATOR
      use RNDGEN
      implicit none
    !  INCLUDE 'ness_common.inc'

      INTEGER*4 I,J,M,IS,IL
      character*128 S
      INTEGER*4 IARGC
      REAL(KIND(1.D0)) :: AUX

C Handle command-line options
      M=IARGC()

C Derive path to default setup files from executable pathname

      call GETARG(0,S)
      call ConvertToSlash(S)
      i=INDEX(S,'/'//'bin/')
      IF (i.gt.1) then
        CFGPATH=S(1:i)//'setup/'
      ELSE
        CFGPATH='setup/'
      ENDIF
      M=IARGC()
      DO I=1,M
        call GETARG(I,S)
! disable console input (for Windows wrapper)
        if (S(1:7).eq.'-nocons') THEN
          write(SOUT,*) '-nocons, console disabled'
          ConsoleEnabled=.FALSE.
        endif
        if (S(1:5).eq.'-dir=') THEN
           CALL BOUNDS(S,IS,IL)
           IS=IS+5
           IL=IL-5
           IF(IL.GT.0) THEN
             RESPATH=S(IS:IS+IL-1)
             call ConvertToSlash(RESPATH)
             IF(RESPATH(IL:IL).NE.'/') THEN ! add path delimiter
                RESPATH=RESPATH(1:IL)//'/'
                IL=IL+1
             ENDIF
             write(SOUT,*) 'dir='//RESPATH(1:IL)
           ENDIF
        else if (S(1:5).eq.'-dll=') then
           write(SOUT,*) 'dll mode'
C// generate XML output for GUI
        else if (S(1:7).eq.'-xmlout') THEN
            write(SOUT,*) 'XML mode'
            XMLOUT=1
        endif
10      if (S(1:2).eq.'-s') THEN
           read(S(3:30),*) J
           if(J.NE.0) THEN
               ISEED=ABS(J)
               write(SOUT,*) 'SEED=',ISEED
           ENDIF
        endif
        if (S(1:2).eq.'-t') THEN
           read(S(3:30),*,err=20) J
           CALL RAN1SEED(ISEED)
           write(*,*) 'Test of the random number generator:'
           call RAN1TEST(J,1000000*J)
20         call RSXSTOP
           RETURN
        endif
      !  if (S(1:4).eq.'-flx') THEN
      !    CALL READ_FLUX(S(5:))
      !  endif
        if (S(1:6).eq.'-dlam=') THEN
          read(S(7:30),*,err=21) AUX
          GEN_DLAM=AUX
        endif
        if (S(1:8).eq.'-rprobe=') THEN
          READ(S(9:),*,ERR=21) AUX
          call SETPROBE_RANGE(AUX)
        endif
21      if (S(1:8).eq.'-nprobe=') THEN
          READ(S(9:),*,ERR=22) AUX
          call SETPROBE_FREQ(NINT(AUX))
        endif
22      if (S(1:8).eq.'-limits=') THEN
          READ(S(9:),*,ERR=23) AUX
          call SET_LIMITS(AUX)
        endif
23      if (S(1:4).eq.'-pulse=') THEN
          ! CALL READ_PULSE(S(8:))
          write(*,*) '-pulse option is no more accepted'
        endif
        if (S(1:4).eq.'-flh') THEN
           write(*,*) '-flh option is no more accepted'
           !read(S(5:30),*,err=30) FLXH
           !FLXH=FLXH*PI/180
        endif
30      if (S(1:4).eq.'-flv') THEN
          write(*,*) '-flv option is no more accepted'
          !read(S(5:30),*,err=50) FLXV
          !FLXV=FLXV*PI/180
        endif
50      if (S(1:6).eq.'-Voigt') THEN
           MDIST=cmos_voigt
           write(*,*) 'pseudo-Voigt mosaic distribution'
        else if (S(1:7).eq.'-Lorenz') THEN
           MDIST=cmos_lorentz
           write(*,*) 'Lorenzian mosaic distribution'
        else if (S(1:4).eq.'-Uni') THEN
           MDIST=cmos_box
           write(*,*) 'Uniform mosaic distribution'
        else
           MDIST=cmos_gauss
          ! write(*,*) 'Gaussian mosaic distribution'
        endif
        if (S(1:4).eq.'-sil') THEN
           read(S(5:30),*,err=60) SILENT
           write(SOUT,*) 'SILENT=',SILENT
        endif
60      if (S(1:5).eq.'-ran1') THEN
           IGEN=1
           write(*,*) 'Numerical Recipes RAN1 generator'
        endif
        if (S(1:5).eq.'-rand') THEN
           IGEN=2
           write(*,*) 'System random number generator'
        endif
        if (S(1:5).eq.'-help'.or.S(1:1).eq.'?') THEN
          call HLPOUT
        endif
        ENDDO
      END


C
C***********************************************************************
      SUBROUTINE BENCH
C /// simulate flux at the sample (arg<>2) or at the detector (arg=2)
C /// by forward method (ICOM=1) or "from the sample" (ICOM=0)
C***********************************************************************
      USE TRACING
      use XMLSIMRES
      USE CONSTANTS
      IMPLICIT NONE
!     integer :: ierr
!      character(32) CNUM
!      real(KIND(1.D0)) :: Z

c      Z=-PI/1.D10
c      do i=1,20
c        call FLOAT2STR(Z,CNUM)
c        write(*,*) trim(CNUM)
c        Z=Z*1.D1
c      enddo
c      read(*,*)

    !  call NESS_CONV(IERR)

! call TRACING_PROC
!      call WRITE_BEAMLINE(20)

!      call save_config("simres.xml",ierr)
!      call report_classes
      write(*,*) 'BENCH: no action defined'
      END SUBROUTINE BENCH
C

C-----------------------------------------------------------------------------
        SUBROUTINE RESINIT
C// initialize RESTRAX
C// include all actions necessary to allocate memory,
C// initialize variables, print LOGO etc..
C-----------------------------------------------------------------------------
      USE COMPONENTS
      USE TRACINGDATA
      USE CONSTANTS
      use XMLSIMRES
      use XMLWRITER
      use FILETOOLS
      USE XMLINFO
      use GRFEXEC
      use TABLES
      USE CMDHANDLER
      use VCTABLE_COMPONENTS
      use SOURCES_TABLE
      use MIRROR_TABLE
      use MESSAGES
      use CRYSTALS
      IMPLICIT NONE
      INCLUDE 'config.inc'

      INTEGER*4 IRES,IERR
      character(MAX_FNAME_LENGTH) extname,outname,FNAME,JOBNAME
      REAL*8 HMOS_DIST
      EXTERNAL HMOS_DIST

C OS info
      write(*,*) 'SYS=',trim(SYSNAME)

! initialize command buffer
      call INIT_IO
! initialize messages
      call INIT_MESSAGES


C initialize LINP
      call INIT_COMMANDS
      call CLASSES_INIT

! clear lookup tables
      call MIRROR_CLEARALL
      call STRAIN_CLEARALL
      CALL CLEAR_FLUX_TABLE

C Set default parameters (some are changed by command options):
      SILENT=0         ! default silence level:
      IGEN=0           ! default rnd generator: Mersenne Twister
      MDIST=cmos_gauss ! default gaussian mosaic
      XMLOUT=0         ! no XML output by default

C Handle command-line options
! IMPORTANT: PROCARG also determines paths to setup files
      CALL PROCARG

C initialize error function
      select case(MDIST)
      case(cmos_gauss)
        RDIST=4.D0
      case(cmos_box)
        RDIST=1.D0
      case default
        RDIST=6.D0
      end select
      CALL ERF_INIT(HMOS_DIST,-RDIST,RDIST)

      CALL RAN1SEED(ISEED)       ! Initialize random number generator
      if (XMLOUT.eq.0) CALL READINIFILE(JOBNAME)  ! read restrax.ini file (only in console mode)
      call SETRESPATH(RESPATH) ! set default path for configuration
      ! call SETDATAPATH(DATPATH)  ! set data path to current dir.
      call LISTGRFDEV   ! fill the list with available grahics devices
      CALL GETENV('PGPLOT_DEV',FNAME) ! select device according to the PGPLOT_DEV variable
      IF(FNAME(1:1).NE.' ') THEN
        call SELGRFDEV(trim(FNAME))
      else
        call SELGRFDEV('/NULL')
      ENDIF
! loading of classes definition and defualt instrument setup is obligatoy
! Stop the program in case of failure
      call parse_xmlclasses(trim(CFGPATH)//trim(XCONFIG_CLASSES),ierr)
      !write(*,*) 'RESINIT parsed classes: ',ierr
      if (ierr.ne.0) then
        write(*,*) 'RESINIT: Fatal error in parse_xmlclasses.'
        call RESEND
      endif
! read crystal.lib AFTER classes and BEFORE configuration !
      call READ_ATOMS
      !write(*,*) 'RESINIT READ_ATOMS OK'
      !write(*,*) 'going to read ',trim(CFGPATH)//trim(XCONFIG_CRYSTALS)
      call parse_xmlcrystals(trim(CFGPATH)//trim(XCONFIG_CRYSTALS),ierr)
      !write(*,*) 'RESINIT parse_xmlcrystals OK'
      call READCRYST
      !write(*,*) 'RESINIT READCRYST OK'
      if (XMLOUT.le.0) then
! print logo
      CALL LOGO
! read configuration, try (1) XCONFIG_NAME and (2) XCONFIG_DEFAULT
! Stop the program in case of failure
        call load_config(trim(XCONFIG_NAME),ierr)
        if (ierr.ne.0) call load_config(trim(CFGPATH)//trim(XCONFIG_DEFAULT),ierr)
        !write(*,*) 'RESINIT parsed config: ',ierr
        if (ierr.ne.0) then
          write(*,*) 'RESINIT: Fatal error in load_config.'
          call RESEND
        endif
      ELSE
! Send version info to GUI.
! Also inform GUI that RESTRAX is ready to accept commands
! GUI then streams the initial configuration to RESTRAX
        !write(*,*) 'RESINIT calling VERSIONINFO'
        call VERSIONINFO
        !write(*,*) 'VERSIONINFO done'
      endif

! execute initial job file
      if (XMLOUT.le.0) then
    ! run job file required by restrax.ini
        IF (JOBNAME(1:1).NE.' ') THEN
          CALL REINP(JOBNAME)
          CALL LINPSETIO(SINP,SOUT,SMES)
          RETURN
        ELSE
    ! ask for a job file: console mode only
          write(*,*) 'Input batch file (press ENTER to skip)'
          call DLG_FILEOPEN(' ','|'//trim(RESPATH)//'|'//trim(CFGPATH),'inp',1,0,IRES,extname)
          if (ires.gt.0) then
            write(*,*) 'Output from batch execution (press ENTER for console output)'
            call DLG_FILEOPEN(' ',' ','out',1,1,IRES,outname)
            if (IRES.gt.0) then
              CALL REOUT(trim(outname))
            endif
            CALL REINP(trim(extname))
            write(SOUT,*) 'RESTRAX - batch job '//trim(extname)
            CALL LINPSETIO(SINP,SOUT,SMES)
          endif
        ENDIF
      endif
      call FLUSH_MESSAGES
      write(*,*) 'RESINIT OK'
      END SUBROUTINE RESINIT

!-----------------------------------------------------------------------------
      SUBROUTINE RESEND
! end of RESTRAX
! include all actions necessary to deallocate memory etc...
!-----------------------------------------------------------------------------
      use FILETOOLS
      USE DATASETS
      use CLASSES
      use CLASSDEF
      use TABLE_CRYSTALS
      use XTALS_REF
      USE NSTORE
      USE VCTABLE_COMPONENTS
      use MESSAGES
      IMPLICIT NONE
      WRITE(SMES,*) 'RESEND called'
      ! release dynamically loaded libraries
      call RELEASEMCPL
      ! restore IO
      CALL REINP(' ')
      CALL REOUT(' ')
    ! deallocate dynamically allocated data
      !call NSTORE_DESTROY
      call DSET_CLEAR ! also destroys associated NSTORE data
      call NSTORE_UNREGISTER_ALL ! destroy the remaining NSTORE data
      call FREE_IO    ! command buffer

      call FREE_MESSAGES
      call CLEAR_CLASSES
      call CLEAR_TABLE_CRYSTALS
      call XTALS_REF_DISPOSE
      call PCOM_DISPOSE_ALL
      WRITE(SMES,*) ' -> End of ResTrax'
      call XML_RSXDUMP(SMES,'EXIT',1)
      call XML_RSXDUMP(SMES,'EXIT',0)
      STOP
      END SUBROUTINE RESEND


C-------------------------------------
      SUBROUTINE RESTRAX_MAIN
C Main unit for console application
C Should be called by the main procedure
C-------------------------------------
      USE IO
      IMPLICIT NONE
      INCLUDE 'linp.inc'
      CHARACTER*128 LINE
      character*64 PRMPT
      INTEGER*4 IRES,RESTRAX_HANDLE,IS,IL,PROCESS_BUFFER

2     FORMAT(a,$)
C initialization
      CALL RESINIT
C run
      IRES=0
      DO WHILE (IRES.EQ.0)
        call RESTRAX_PROMPT(PRMPT)
        call BOUNDS(PRMPT,IS,IL)
c// PRMPT has char(0) at the end
        if (IL.GT.2) WRITE(linp_out,2) PRMPT(IS:IS+IL-2)
        IRES=PROCESS_BUFFER() ! process buffered input
        IF (IRES.EQ.0) THEN
          call RSXREAD(LINE)
          IRES=RESTRAX_HANDLE(LINE)
          IF (SOUT.NE.6) call flush(SOUT)
        ENDIF
      ENDDO
C finalization
      write(*,*) 'RESTRAX_MAIN finished'
      CALL RESEND
      END SUBROUTINE RESTRAX_MAIN

