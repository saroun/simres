C//////////////////////////////////////////////////////////////////////
C////
C////  R E S T R A X   4.1
C////
C////  Some subroutines for I/O operations:
C////
C////
C//////////////////////////////////////////////////////////////////////

C-------------------------------------------
      SUBROUTINE READ_RESCAL(IU,IERR)
C Read RESCAL parameters from unit IU (assume RESCAL format)
C IERR<>0 ... error while reading file
C-------------------------------------------
      use XMLINFO
      IMPLICIT NONE
      INCLUDE 'const.inc'
      INCLUDE 'rescal.inc'

      INTEGER*4 IU,IERR,IE
      CHARACTER*30 LINE
      character*128 MSG
      character*256 CLINE
      INTEGER*4 I,J,IS,IL
      REAL*8 VER,Z
102   FORMAT('Error ',I5,' in RESCAL file, line=',I5)
1      FORMAT(a)
      IERR=0
      I=1
      VER=0.0
      READ(IU,1,ERR=98,END=98,IOSTAT=IERR) LINE
      CALL READ_R8('version',LINE,VER,IERR)
      if ((IERR.NE.0).OR.(VER.LT.4.77)) then ! till  version 4.77
        READ(LINE,*,ERR=98,END=98,IOSTAT=IERR)  RES_DAT(1)
        DO I=2,i_A3-1
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
        RES_DAT(i_A3)=0.0
        RES_DAT(i_A4)=0.0
        DO I=i_A4+1,i_DA3-1
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
        RES_DAT(i_DA3)=0.0
        RES_DAT(i_DA4)=0.0
        DO I=i_DA4+1,i_CHAN-1
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
        RES_DAT(i_CHAN)=1  ! fc_ic0
        DO I=i_CHAN+1,RES_NVAR
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
      else IF (VER.LE.5.05) THEN    ! till version 5.05
        DO I=1,i_A3-1
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
        RES_DAT(i_A3)=0.0
        RES_DAT(i_A4)=0.0
        DO I=i_A4+1,i_CHAN-1
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
        RES_DAT(i_CHAN)=1 ! fc_ic0
        DO I=i_CHAN+1,RES_NVAR
          READ(IU,*,ERR=98,END=98,IOSTAT=IERR) RES_DAT(I)
        ENDDO
      ELSE   ! current version, read format NAME=number
        J=0
        DO WHILE (IERR.EQ.0)
          READ(IU,1,ERR=98,END=90,IOSTAT=IERR) CLINE
          CALL BOUNDS(RES_NAM(i),IS,IL)
          DO I=1,RES_NVAR
            CALL READ_R8(RES_NAM(i)(IS:IS+IL-1),CLINE,Z,IE)
            if (IE.EQ.0) then
              RES_DAT(I)=Z
              J=J+1
            endif
          enddo
        ENDDO
90      if (J.LT.RES_NVAR) then
          call MSG_WARN('Input file contains incomplete set of parmeters',1)
        endif
      ENDIF
      IERR=0
      return
c// handle errors
98    if (IERR.NE.0) then
        write(MSG,102) IERR,I
        call MSG_ERROR('READ_RESCAL',MSG,0,1)
      endif
      END SUBROUTINE READ_RESCAL

C-------------------------------------------
      SUBROUTINE OPEN_RESCAL(FNAME,IERR)
C Read RESCAL parameters from a *.res file
C-------------------------------------------
      USE CONSTANTS
      USE FILETOOLS
      IMPLICIT NONE
      INCLUDE 'const.inc'
      INCLUDE 'rescal.inc'
      character(*),intent(in) :: FNAME
      integer,intent(out) :: IERR
      character(MAX_FNAME_LENGTH) :: FN,FRES,FRESPATH
      INTEGER :: IU,IPROMPT
      REAL(KIND(1.D0)) :: DAT(MRES)

C make a copy of RESCAL parameters
      DAT(1:RES_NVAR)=RES_DAT(1:RES_NVAR)
C open file
      if (len_trim(FNAME).le.0) then
        FN=adjustl(RESCAL_NAME)
        IPROMPT=1
      else
        FN=adjustl(FNAME)
        IPROMPT=0
      endif
      call OPENRESFILE(trim(FN),'res',IPROMPT,FRES,FRESPATH,IU,IERR)
      if (IERR.eq.0) then
        call READ_RESCAL(IU,IERR)
        CLOSE(IU)
      else
        call MSG_ERROR('OPEN_RESCAL',' Can''t open RESCAL file '//trim(FN),0,1)
      endif
      if (ierr.ne.0) then
        RES_DAT(1:RES_NVAR)=DAT(1:RES_NVAR)
      else
        call RESCAL_IMPORT
        RESCAL_NAME=trim(FN)
        call MSG_INFO(' Rescal parameters imported from '//trim(FRES),1)
      endif
      END SUBROUTINE OPEN_RESCAL


