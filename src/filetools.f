!//////////////////////////////////////////////////////////////////////
!////  $Id: filetools.f,v 1.33 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.33 $
!////     $Date: 2019/08/16 17:16:25 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Input/Output operations on files
!////  - ensure search order for config. files
!////  - save info about current config. files
!////
!//////////////////////////////////////////////////////////////////////
      MODULE FILETOOLS
      use XMLINFO
      implicit none


      integer,parameter :: MAX_FNAME_LENGTH=256
      character(MAX_FNAME_LENGTH),parameter :: XCONFIG_DEFAULT='default.xml'

! search paths
    !  CHARACTER(MAX_FNAME_LENGTH) :: DATPATH='./'   ! exp. data
      CHARACTER(MAX_FNAME_LENGTH) :: RESPATH   ! configurations, lookup tables, etc., user-defined
      CHARACTER(MAX_FNAME_LENGTH) :: CFGPATH   ! configurations, lookup tables, etc., system-defined (./setup)
      CHARACTER(MAX_FNAME_LENGTH) :: OUTPATH   ! output data

! config. file names
      character(MAX_FNAME_LENGTH) :: XCONFIG_NAME=XCONFIG_DEFAULT      ! XML configuration file name
      character(MAX_FNAME_LENGTH) :: XCONFIG_CLASSES='classes.xml'     ! XML file with class definitions
      character(MAX_FNAME_LENGTH) :: XCONFIG_COMMANDS='commands.xml'   ! XML commands definition
      character(MAX_FNAME_LENGTH) :: XCONFIG_CRYSTALS='crystals.xml'   ! XML file with cdefinitions of crystals
      character(MAX_FNAME_LENGTH) :: RESCAL_NAME='default.res'       ! Rescal file name (old format)
      character(MAX_FNAME_LENGTH) :: CFGNAME='default.cfg'           ! CFG file name (old format)

      CHARACTER(128) :: CFGTITLE='default configuration'

      logical :: checkOverwrite=.true.

      ! temporary info for 'file' read rediected from GUI
      integer :: file_iread=0
      integer :: file_id=0
      character(512) :: file_name='default'

! default file IO flags
      integer,parameter :: i_OK=0
      integer,parameter :: i_EOF=1
      integer,parameter :: i_VALID=2
      integer,parameter :: i_ERR=3


      contains



!---------------------------------------------------------------
      SUBROUTINE READ_TABLE_ROW(IU,iline,LINE,ierr)
! Read one line from IU, ignore comments and emty lines
! IU ... input unit
! LINE ... line text (empty for comments or errors)
! ierr = 0 OK
! ierr = 1 EOF
! ierr = other error
!---------------------------------------------------------------
      integer,intent(in) :: IU
      character(*),intent(out) :: LINE
      integer,intent(out) :: iline,ierr
      character(1024) :: S
      integer :: IS,IL,IC
1     format(a)

      IL=0
      do while(IL<=0)
        iline=iline+1
        S=' '
        Read(IU,1,iostat=IERR,end=30,err=40) S
        CALL BOUNDS(S,IS,IL)
        IC=index(S,'#')
        if (IC>0) IL=IC-IS
      enddo
      LINE=S(IS:IS+IL-1)
      ierr=i_OK
      Return
30    ierr=i_EOF
      LINE=' '
      return
40    ierr=i_ERR
      LINE=' '
      return
      end SUBROUTINE READ_TABLE_ROW

!--------------------------------------------------
      integer function FILE_READF(IU,A,MA,IMIN,IMAX)
! Read float numbers into a double array in given range of indexes
! return
! <0 error or EOF
! 0 empty line
! >0 number of elements actually read
!-----------------------------------------------------
      integer,intent(in) :: IU,IMIN,IMAX,MA
      real(kind(1.D0)),intent(out) :: A(MA)
      character(256) :: CLINE
      integer :: i,ird,ierr,nr
      if ((IMIN<=0).or.(imax>MA)) then
        FILE_READF=-2
        return
      endif
      nr=0
      do i=IMIN,IMAX
        ird=FILE_READ(IU,CLINE)
        if (ird<=0) then ! error
          ierr=ird
          write(*,*) 'FILE_READF error ',ird,trim(CLINE)
          goto 98
        else ! valid input
          READ(CLINE,*,ERR=98,END=10,IOSTAT=IERR) A(i)
10        nr=nr+1
        endif
      enddo
      FILE_READF=nr
      return
98    FILE_READF=ierr
      end function FILE_READF


!--------------------------------------------------
      integer function FILE_READ(IU,LINE)
! Read one line from buffer using RSXREAD
! NOTE: file IO in restrax is redirected via GUI
! IU unit number is not necessary, but we use it to check the
! consistency of input sequence (OPEN, READ, CLOSE ...)
! RETURN:
! >=0   ... number of characters read (trimmed)
! -1    ... end of file
! -2    ... file already closed or unit number mismatch
! -3    ... not in buffer mode (buffer empty)
!-----------------------------------------------------
      !INCLUDE 'restrax_cmd.inc'
      CHARACTER*(*) :: LINE
      integer,intent(in) :: IU
      integer :: ires,id,ierr
      LINE=' '
      if (IU==file_iread) then
        if (InputMode.EQ.imBUFFER) then
          call RSXREAD(LINE)
          ires=len_trim(LINE)
          if (ires>0) then
            call READ_I4('CLOSEFILE',LINE,id,ierr)
            if (ierr==0) then
              if (id==file_id) then
                file_iread=0
                file_name=' '
                ires=-1
              endif
            endif
          endif
        else
          ires=-3
        endif
      else
        ires=-2 ! read unit closed or unit number mismatch
      endif
      !write(*,*) 'FILE_READ ',trim(LINE),' ',ires,' ',InputMode.EQ.imBUFFER
      ! close unit automaticaly when error or end of file encountered
      if (ires<0) call CLOSEFILE_READ(IU)
      FILE_READ=ires
      end function FILE_READ

!--------------------------------------------------
      integer function OPENFILE_READ(FNAME,IU)
! Open file for read
! NOTE: file IO in restrax is redirected via GUI
! Only register file name and unit number
!--------------------------------------------------
      CHARACTER(*),intent(in) :: FNAME
      integer,intent(in) :: IU
      character(512) :: MSG,LINE
      integer :: ierr,id
1     format('File unit is already open: iu=',I5,', name=',a)
2     format('Wrong file input format, missing OPENFILE=ID line. Check sending GUI procedure.')
      if (file_iread==0) then
        file_iread=IU
        if (len_trim(FNAME)>0) then
          file_name=trim(FNAME)
        else
          file_name='default'
        endif
        ierr=FILE_READ(file_iread,LINE)
        if (ierr>0) call READ_I4('OPENFILE',LINE,id,ierr)
        if (ierr.ne.0) then
          write(MSG,2)
          call MSG_ERROR('OPENFILE_READ',trim(MSG),1,1)
          call CLOSEFILE_READ(IU)
        endif
        file_id=id
        !write(*,*) 'OPENFILE_READ ',IU,' ',trim(file_name)
      else
        write(MSG,1) file_iread,trim(file_name)
        call MSG_ERROR('OPENFILE_READ',trim(MSG),1,1)
        call CLOSEFILE_READ(IU)
      endif
      OPENFILE_READ=file_iread
      end function OPENFILE_READ


!--------------------------------------------------
      subroutine CLOSEFILE_READ(IU)
! Open file for read.
! NOTE: file IO in restrax is redirected via GUI
! Only reset file unit and name
!--------------------------------------------------
      !include 'restrax_cmd.inc'
      integer,intent(in) :: IU
      character(512) :: LINE
      logical :: iop
      integer :: id,ierr
      !write(*,*) 'CLOSEFILE_READ ',IU,' ',trim(file_name),' ',file_iread
      if (file_iread==0) return
      file_iread=0
      file_name=' '
      iop=(InputMode.EQ.imBUFFER)
      ! clean buffer from the rest of file data
      ! may be not empty in case of a read error
      do while (iop)
        call RSXREAD(LINE)
        if (len_trim(LINE)>0) then
          call READ_I4('CLOSEFILE',LINE,id,ierr)
          !write(*,*) 'CLOSEFILE_READ ',ierr,id,' ',trim(LINE)
          iop=((ierr.ne.0).or.(id.ne.file_id))
        endif
        iop=(iop.and.(InputMode.EQ.imBUFFER))
      enddo
      end subroutine CLOSEFILE_READ


!--------------------------------------------------
      integer function OPENFILE_WRITE(FNAME)
! Open file for write and return unit number.
! NOTE: file IO in restrax is redirected via GUI
! Only send start XML tags to STDOUT. File save is handled by GUI !!
!--------------------------------------------------
      CHARACTER(*),intent(in) :: FNAME
      call XML_PROGRESS(smes,0,10,5,'')
      call XML_OPENWSTART(fname)
      OPENFILE_WRITE=smes
      end function OPENFILE_WRITE


!--------------------------------------------------
      SUBROUTINE CLOSEFILE_WRITE(IU)
! Close file for write and return unit number
! NOTE: file IO in restrax is redirected via GUI
! Only send closing XML tags to STDOUT. File save is handled by GUI !!
!--------------------------------------------------
      integer,intent(in) :: IU
      call XML_OPENWEND
      call XML_PROGRESS(smes,2,0,0,'')
      end SUBROUTINE CLOSEFILE_WRITE


!--------------------------------------------------
      integer function OPENFILEUNIT(FNAME,mustexist)
! Open file for input or output
! return IO unit>0 if successful
! return 6 if FNAME=STDOUT
!--------------------------------------------------
      IMPLICIT NONE
      CHARACTER(*),intent(in) :: FNAME
      logical,intent(in) :: mustexist
      integer :: lun,i,ierr
      logical :: isopen,exists
1     format('OPEN [',a,'] len=',G11.5,' char=',G11.5)
      lun=-1

      if (index(FNAME,char(0))>0) then
        write(*,*) 'OPENFILEUNIT WARINIG: file name with null characters'
        cALL DELNUL(FNAME)
      endif

      if (trim(FNAME).eq.'STDOUT') then
        OPENFILEUNIT=6
        return
      else if (len_trim(FNAME).eq.0) then
        OPENFILEUNIT=-1
        return
      endif

      ! write(*,1) trim(FNAME), len_trim(FNAME), ichar(FNAME(1:1))
      do i = 10,99
        inquire(unit=i,opened=isopen)
        if (.not.isopen) then
          if (mustexist) then
            inquire(file=trim(adjustl(fname)),exist=exists)
            if (.not.exists) then
              call MSG_ERROR('OPENFILEUNIT','file does not exist: '//trim(fname),0,1)
            else
              lun = i
            endif
          else
            lun = i
          endif
          if (lun>0) then
            open(unit=lun,file =trim(adjustl(fname)),IOSTAT=ierr)
            if (ierr.ne.0) lun=-1
          endif
          EXIT
        endif
      enddo
      OPENFILEUNIT=lun
      end function OPENFILEUNIT

!-------------------------------------------------------------------------
      integer function OPENFILEOUT_SFX(FNAME,SFX)
! Open file for output with specifed name suffix.
! If there is no path in FNAME, use OUTPATH
! return IO unit>0 if successful
! return 6 if FNAME=STDOUT
!-------------------------------------------------------------------------
      character*(*) :: FNAME,SFX
      character*256 :: fout,fdir,fn
      character*512 :: S
      character*32 :: fe
      integer :: is,il,IO, LSTR
call bounds(fname,is,il)
      if (trim(FNAME).eq.'STDOUT') then
        IO=6
      else
        call bounds(fname,is,il)
        call SplitPathName(fname(is:is+il-1),fdir,fn,fe)
        if (len_trim(fdir).eq.0) fdir=trim(OUTPATH)
        S=trim(fdir)//'/'//trim(fn)
        LSTR=min(256,len_trim(S))
        fout=S(1:LSTR)
        if (len_trim(sfx).gt.0) fout=trim(fout)//'_'//trim(sfx)
        if (len_trim(fe).gt.0) fout=trim(fout)//'.'//trim(fe)
        IO=OPENFILEUNIT(trim(fout),.false.)
      endif
      OPENFILEOUT_SFX=IO
      end function OPENFILEOUT_SFX


C-------------------------------------------------------------------------
      SUBROUTINE OPENRESFILE(FNAME,FEXT,IPROMPT,FRES,FRESPATH,IUNIT,IERR)
C Open file for input in RESTRAX, searching in following directories:
C current:RESPATH:CFGPATH
C INPUT:
C   FNAME ... filename
C   FEXT  ... default extension
C   IUIT  ... file unit number
C   IPROMPT   ... force user input if ird>0, even if fname is not empty
C OUTPUT:
C   FRES  ... resulting filename without path
C   FRESPATH  ... path where the file was found
C   IERR  ... <>0 if cannot open file
C-------------------------------------------------------------------------
      character(*), intent(in) :: FNAME,FEXT
      integer,intent(in) :: IPROMPT
      character(*), intent(out) :: FRES,FRESPATH
      integer,intent(out) :: IUNIT,IERR
      INTEGER :: IRES,IL1,IL2
      character(2*MAX_FNAME_LENGTH) :: pathlist
      character(MAX_FNAME_LENGTH) :: ffn,fn,fe


      IERR=-1
      FRES=' '
      FRESPATH=' '
      pathlist='|'//trim(RESPATH)//'|'//trim(CFGPATH)
      CALL DLG_FILEOPEN(FNAME,pathlist,FEXT,IPROMPT,0,IRES,FFN)
      IF (ires.gt.0) THEN
        IUNIT=OPENFILEUNIT(trim(FFN),.true.)
    !    OPEN(UNIT=IUNIT,FILE=trim(FFN),STATUS='OLD',action='READ',IOSTAT=IERR)
        if (IUNIT>0) IERR=0
        if (IERR == 0) then
          call SplitPathName(trim(FFN),FRESPATH,FN,FE)
          IL1=len_trim(FN)
          IL2=len_trim(FE)
          IF ((IL1.GT.0).AND.(IL2.GT.0)) THEN
            FRES=trim(FN)//'.'//trim(FE)
          ELSE IF (IL1.GT.0) THEN
            FRES=trim(FN)
          ENDIF
        endif
      ENDIF
      END SUBROUTINE OPENRESFILE


!*********************  General file operations and dialogs *******************

C-------------------------------------------------------------------
      SUBROUTINE DLG_FILEOPEN(FNAME,FPATH,FEXT,IPROMPT,IWRITE,IRES,FRES)
C Get a fully qualified filename, test existence, etc.
C Does not actually open any file !!
C INPUT:
C   FNAME  ... filename (accepts also full pathname)
C   FPATH  ... | delimited list of search directories
C   FEXT   ... default extension
C   IPROMPT  ... if>0, read the filename interactively
C   IWRITE ... if=0, open for read, else open for write
C RETURN:
C   FRES   ... resulting filename (incl. path)
C   IRES>0 ... ord. number of the path string from fpath
C   IRES=0 ... not found
C-------------------------------------------------------------------
      CHARACTER(*),intent(in) :: FNAME,FPATH,FEXT
      INTEGER,intent(in) :: IPROMPT,IWRITE
      INTEGER,intent(out) :: IRES
      integer, parameter :: MAXSTR=2*MAX_FNAME_LENGTH+1
      CHARACTER(*),intent(out) :: FRES
      CHARACTER(MAX_FNAME_LENGTH) :: FD,FN,FE,FSEARCH,AUX
      CHARACTER(MAXSTR) :: S
      CHARACTER(1024) :: LINE
      INTEGER :: ILD,ILN,ILE,LEXT,LS,LP,LL,LRES,LSTR
      LOGICAL :: APEXT, ASKNAME, LOG1
1     FORMAT(' Open file [',a,'] : ',$)
2     FORMAT(' Open file : ',$)
3     FORMAT(' Save to file [',a,'] : ',$)
4     FORMAT(' Save to file : ',$)
10    format('DLG_FILEOPEN: ',1x,a,'=[',a,'] len=',G13.5)
12    format('DLG_FILEOPEN length:',6(1x,G12.5))


      !write(*,10) 'FNAME',trim(FNAME),len_trim(FNAME)
      !write(*,10) 'FPATH',trim(FPATH),len_trim(FPATH)
      !write(*,10) 'FEXT',trim(FEXT),len_trim(FEXT)
      if (trim(fname).eq.'STDIN') then
          FRES='STDIN'
          IRES=1
          return
      endif
      IRES=0
      LRES=LEN(FRES)
      FSEARCH=trim(adjustl(FPATH))
      call ConvertToSlash(FSEARCH)
      !write(*,10) 'FSEARCH',trim(FSEARCH)
c split to path, name, extension
      call SplitPathName(FNAME,FD,FN,FE)
      !write(*,10) 'FD',trim(FD),len_trim(FD)
      !write(*,10) 'FN',trim(FN),len_trim(FN)
      !write(*,10) 'FE',trim(FE),len_trim(FE)
c check if FNAME contains full path
c if yes, use it as the only search path, otherwise use FPATH
c if not, use FD/FN as the file name appended to a search path
      if (isFullPathName(FD)) then
        FSEARCH=trim(FD)
      else
        if (len_trim(FD).gt.0) then
          S=trim(FD)//'/'//trim(FN)
          LSTR=min(MAX_FNAME_LENGTH,len_trim(S))
          FN=S(1:LSTR)
        endif
      endif
      LEXT=len_trim(FEXT)
      ILN=len_trim(FN)
      ILE=len_trim(FE)

      ILD=len_trim(FD)

C append extension ? (yes if there is extension suggested AND extension is missing)
      APEXT=((LEXT.GT.0).AND.(ILE.LE.0))
C ask for filename interactively ? (yes if explicitly required OR no filename is given)
      ASKNAME=((ILN.LE.0).OR.(IPROMPT.GT.0))

      ! write(*,12) LEXT,ILN,ILE,ILD,APEXT,ASKNAME

C format prompt
      S=' '
      IF (APEXT) THEN
        IF (ILN.GT.0) THEN
          S=trim(FN)//'.'//trim(FEXT)
        ELSE
          S='*.'//trim(FEXT)
        ENDIF
      ELSE IF (ILN.GT.0) THEN
        if (ILE.GT.0) then
          S=trim(FN)//'.'//trim(FE)
        else
          S=trim(FN)
        endif
      ENDIF
      LS=len_trim(S)
C ask for filename
      IF (ASKNAME) THEN ! get filename interactively
        select case(IWRITE)
        case(0)
          IF (LS.gt.0) THEN
            WRITE(*,1) trim(S)
          ELSE
            WRITE(*,2)
          ENDIF
        case default
          IF (LS.gt.0) THEN
            WRITE(*,3) trim(S)
          ELSE
            WRITE(*,4)
          ENDIF
        end select
        call RSXREAD(AUX)
        !call ConvertToSlash(AUX)
        IF (len_trim(AUX).LE.0) THEN ! use default
          IF (INDEX(S,'*').le.0) THEN
            LINE=trim(S)
          ELSE
            RETURN
          ENDIF
        else
          LINE=trim(adjustl(AUX))
        ENDIF
        call SplitPathName(LINE,FD,FN,FE)
        if (isFullPathName(FD)) then
          FSEARCH=trim(FD)
        else
          if (len_trim(FD).gt.0) then
            S=trim(FD)//'/'//trim(FN)
            LSTR=min(MAX_FNAME_LENGTH,len_trim(S))
            FN=S(1:LSTR)
          endif
        endif
        ILN=len_trim(FN)
        ILE=len_trim(FE)
        ILD=len_trim(FD)
        APEXT=((LEXT.GT.0).AND.(ILE.LE.0)) ! append extension?
      ENDIF
c construct filename with extension (excl. full path)
c append extension when required
      if (APEXT) then
        S=trim(FN)//'.'//trim(FEXT)
      else if (ILE.GT.0) THEN
        S=trim(FN)//'.'//trim(FE)
      else
        S=trim(FN)
      endif
      LSTR=min(MAX_FNAME_LENGTH,len_trim(AUX))
      AUX=S(1:LSTR)

      select case(IWRITE)
  ! open for read
      case(0)
    ! find the first file that exists in a directory listed in FDIR
        CALL CHECKRESFILE(AUX,FSEARCH,IRES,FRES)
        if (IRES.le.0) call MSG_WARN('File not found: '//trim(AUX),1)
  ! open for write
      case default
    ! prepend the path name
        IRES=1
        LP=len_trim(FSEARCH)
        IF ((LP.GT.0).AND.(FSEARCH(LP:LP).NE.'/')) THEN
          S=trim(FSEARCH)//'/'//trim(AUX)
        ELSE IF (LP.GT.0) THEN
          S=trim(FSEARCH)//trim(AUX)
        ELSE
          S=trim(AUX)
        ENDIF
        call STRIM(S,FN)

        LL=min(LRES,len_trim(FN))
        FRES=FN(1:LL)
C Overwrite confirmation only for IOVER>0
    !    write(*,*) 'DLG_FILEOPEN ',checkOverwrite,IRES
        IF (checkOverwrite) THEN
          S='y'
          INQUIRE(FILE=FN(1:LL),EXIST=LOG1)
          IF (LOG1) THEN  ! ask before overwrite
            write(*,'("Overwrite ",a,"? [yes/no]")') FN(1:LL)
            call RSXREAD(S)
          ENDIF
          select case(S(1:1))
          case('y','Y')
            IRES=1
          case default
            IRES=0
          end select
        ENDIF
      end select

      END SUBROUTINE DLG_FILEOPEN

!-----------------------------------------------------------------------
      SUBROUTINE DLG_GETPATH(SARG,SPROMPT,ANSWER,PNAME)
! Get directory name from dialog or from the argument
! Existence is not checked
! INPUT:
!   SARG      ... input string with the path name
!   PROMPT    ... input prompt text, results in: > prompt [default] : _
!   PNAME     ... default name (used on plain ENTER)
! RETURN:
!   PNAME     ... resulting pathname
!-----------------------------------------------------------------------
      CHARACTER(*),intent(in) :: SARG,SPROMPT,ANSWER
      CHARACTER(*),intent(inout) :: PNAME
      CHARACTER(MAX_FNAME_LENGTH) :: MYPATH
      INTEGER :: LL
      logical :: ASK

1     FORMAT(a,' [',a,'] : ',$)
2     FORMAT(a,' [current directory] : ',$)
3     FORMAT(a,' ',a)

      LL=len_trim(PNAME)
      ASK=(len_trim(SARG).le.0)

! run dialog if no argument is given
      IF (ASK) THEN
        IF (LL.GT.0) THEN
          write(*,1) trim(SPROMPT),trim(PNAME)
        ELSE
          write(*,2) trim(SPROMPT)
        endif
        call RSXREAD(MYPATH)
        IF (len_trim(MYPATH).eq.0) THEN
          MYPATH=trim(adjustl(PNAME))
        ENDIF
      ELSE
        MYPATH=trim(adjustl(SARG))
      ENDIF
      call ConvertToSlash(MYPATH)
      PNAME=trim(MYPATH)
      LL=len_trim(MYPATH)
C Interpret MYPATH, ensure that ending / is present
  ! empty DIRNAME = current directory
      IF ((LL.LE.0).OR.(trim(MYPATH).EQ.'.').OR.(trim(MYPATH).EQ.'./')) THEN
        PNAME=' '
        write(*,3) trim(ANSWER),'current folder'
      else
        IF(MYPATH(LL:LL).NE.'/') PNAME=trim(MYPATH)//'/'
        write(*,3) trim(ANSWER),trim(PNAME)
      endif
      END SUBROUTINE DLG_GETPATH

C-----------------------------------------------------------------------
      SUBROUTINE SETRESPATH(SARG)
C select search path for configuration files
C-----------------------------------------------------------------------
      CHARACTER*(*) SARG
      CHARACTER*128 PROMPT,ANSWER

      PROMPT=' Search path for configuration files'
      ANSWER=' Configurations will be searched in'
      CALL DLG_GETPATH(SARG,PROMPT,ANSWER,RESPATH)
      if (XMLOUT.LE.0) then
        call XML_RSXDUMP(sout,' ',1)
        call XML_TAG(sout,'CPATH',trim(RESPATH),2)
        call XML_RSXDUMP(sout,' ',0)
      endif
      END SUBROUTINE SETRESPATH

C-----------------------------------------------------------------------
      SUBROUTINE SETOUTPATH(SARG)
C select search path for output
C-----------------------------------------------------------------------
      CHARACTER*(*) SARG
      CHARACTER*128 PROMPT,ANSWER

      PROMPT=' Search path for output files'
      ANSWER=' Output will be saved in'
      CALL DLG_GETPATH(SARG,PROMPT,ANSWER,OUTPATH)
      if (XMLOUT.LE.0) then
        call XML_RSXDUMP(sout,' ',1)
        call XML_TAG(sout,'OPATH',trim(OUTPATH),2)
        call XML_RSXDUMP(sout,' ',0)
      endif
      END SUBROUTINE SETOUTPATH


C--------------------------------------------------------
      SUBROUTINE CHECKRESFILE(FNAME,FPATH,IRES,FRES)
C Test existence of a file
C INPUT:
C   FNAME  ... filename
C   FPATH  ... | delimited list of search directories
C RETURN:
C   FRES   ... resulting filename (incl. path)
C   IRES>0 ... ord. number of the path string from fpath
C   IRES=0 ... not found
C--------------------------------------------------------
      CHARACTER(*),intent(in) :: FNAME,FPATH
      CHARACTER(*), intent(out) ::FRES
      integer,intent(out) :: IRES
      INTEGER*4 IS,IL,ISF,ILF,LL,LRES,J,IP
      LOGICAL*4 LOG1
      CHARACTER(MAX_FNAME_LENGTH) :: FFN

!      write(*,*) 'CHECKRESFILE file=['//trim(FNAME)//']'
!      write(*,*) 'CHECKRESFILE path=['//trim(FPATH)//']'
      IRES=-1
      CALL BOUNDS(FNAME,ISF,ILF)
      LRES=LEN(FRES)
      LOG1=.FALSE.
      IL=1
      FFN=' '
      J=1
      IP=0 ! index of a parameter in FPATH
      LL=1
      DO WHILE (.NOT.LOG1.AND.IL.GE.0)
        IP=IP+1
        CALL FINDSTRPAR(fpath,'|',IP,IS,IL)
        IF (IL.GE.0) THEN
          j=is+il-1 ! last character of the IP-th path sring
          IF (IL.EQ.0) THEN ! no path string
            FFN=fname(ISF:ISF+ILF-1)
            LL=ILF
          ELSE IF ((J.GT.0).AND.(fpath(j:j).NE.'/'))  THEN
            FFN=fpath(is:j)//'/'//fname(ISF:ISF+ILF-1)
            LL=j-is+2+ILF
          ELSE ! last character is PATH delimiter
            FFN=fpath(is:j)//fname(ISF:ISF+ILF-1)
            LL=j-is+1+ILF
          ENDIF
          INQUIRE(FILE=FFN,EXIST=LOG1)
!          write(*,*) 'CHECKRESFILE inquire=['//trim(FFN)//']',LOG1
        ENDIF
      ENDDO
      IF (LL.GT.LRES) LL=LRES  ! clip result size to the output buffer length
      FRES=FFN(1:LL)
      IF (LOG1) IRES=IP
      END SUBROUTINE CHECKRESFILE




!-------------------------------------------
! Return true, if fname is the full path name
! Does not check existence
!-------------------------------------------
      logical function isFullPathName(fname)
      character(*) fname
      logical :: LOG1
        LOG1=.false. ! (len_trim(fname).gt.0)
        LOG1=(LOG1.or.(fname(1:1).eq.'/'))
        LOG1=(LOG1.or.(fname(1:2).eq.'//'))
        LOG1=(LOG1.or.(fname(2:2).eq.':'))
        isFullPathName=LOG1
      end function isFullPathName


!--------------------------------------------------------
! Convert backslash to slash
! NOTE: windows should recognize slash as a delimiter
!--------------------------------------------------------
      SUBROUTINE ConvertToSlash(PNAME)
      CHARACTER(*) PNAME
      integer :: i
      i=index(PNAME,'\')
      do while (i.gt.0)
        if (i.gt.0) PNAME(i:i)='/'
        i=index(PNAME,'\')
      enddo
      end SUBROUTINE ConvertToSlash

C-----------------------------------------------------------------
C Split the full filename (LINE) to path, name and extension
C path includes terminal delimiter, extension includes initial dot
C INPUT:
C line   ... full filename
C dlm    ... path delimiter
C OUTPUT:
C fpath   ... path
C fname   ... name
C fext    ... extension
C-----------------------------------------------------------------
      SUBROUTINE SplitPathName(PNAME,FPATH,FNAME,FEXT)
      CHARACTER(*) PNAME,FPATH,FNAME,FEXT
      INTEGER*4 L,ip,ie,i
      character(MAX_FNAME_LENGTH) :: AUX
      AUX=trim(adjustl(PNAME))
      call ConvertToSlash(AUX)
      L=len_trim(AUX)
! ip .. position of the last path delimiter
      ip=0
      i=INDEX(AUX(1:L),'/')
      if (i.gt.0) then
        DO WHILE (i.GT.0)
          ip=ip+i
          i=INDEX(AUX(ip:L),'/')
        ENDDo
        ip=ip-1
      endif
! ie .. position of the last dot (extension delimiter)
      ie=0
      i=INDEX(AUX(1:L),'.')
      if (i.gt.0) then
        DO WHILE (i.GT.0)
          ie=ie+i
          i=INDEX(AUX(ie:L),'.')
        ENDDo
        ie=ie-1
      endif
! split name
      if (ie.le.0) ie=L+1
      FPATH=' '
      FNAME=' '
      FEXT=' '
      if (ip.gt.0) FPATH=AUX(1:ip-1)
      FNAME=AUX(ip+1:ie-1)
      if (ie.lt.L) FEXT=AUX(ie+1:L)
      end SUBROUTINE SplitPathName


!-------------------------------------------------------------
      SUBROUTINE REINP(SARG)
! redirection of input
!-------------------------------------------------------------
      character(MAX_FNAME_LENGTH) :: FRES,FRESPATH
      character*(*) SARG
      INTEGER*4 IU,IERR
      if((SARG(1:1).ne.' ').and.(SARG(1:1).ne.char(0))) THEN
        call OPENRESFILE(SARG,'inp',0,FRES,FRESPATH,IU,IERR)
        IF (IERR.NE.0) then
          write(smes,*) 'Cannot open input file '//SARG
          IU=5
        endif
      else
        IU=5
      endif
      call CHANGE_INPUT(IU)
      if (IU.eq.5) then
        call XML_SWITCH(SOUT,'BATCHMODE',0)
        call MSG_INFO('Leaving batch mode',1)
      else
        call XML_SWITCH(SOUT,'BATCHMODE',1)
        call MSG_INFO('Entering batch mode: '//trim(SARG),1)
      endif
      END SUBROUTINE REINP

!-------------------------------------------------------------
      SUBROUTINE REOUT(SARG)
! redirection of input
!-------------------------------------------------------------
      character*(*) SARG
      integer :: IU
      if((SARG(1:1).ne.' ').and.(SARG(1:1).ne.char(0))) THEN
        IU=OPENFILEUNIT(trim(SARG),.false.)
        if (IU.le.0) then
          write(smes,*) 'Cannot open output file  '//SARG
          IU=-1
        endif
      else
        IU=6
      endif
      call CHANGE_OUTPUT(IU)
      END SUBROUTINE REOUT

      end module FILETOOLS
