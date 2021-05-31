!//////////////////////////////////////////////////////////////////////
!////  $Id: sample_poly_table.f,v 1.5 2019/06/22 20:48:04 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.5 $
!////     $Date: 2019/06/22 20:48:04 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Lookup tables forpolycystal reflections
!////
!////////////////////////////////////////////////////////////////////////
      MODULE SAMPLE_POLY_TABLE
      use SIMATH
      use FILETOOLS
      implicit none
      save
      private

! reflections table
      integer,parameter :: PCRYST_MAXREF=512
      character(27),parameter :: DEFREFPAR='SIGA:SIGI:SIGF:SIGS:BT:BMPH'
      ! reflection parameters:
      ! dhkl,jhkl,Fhkl = dhkl[A], multiplicity, FHKL/cell_volume [fm/A^3]
      ! reflection cross-section in [1/mm]
      ! s_dhkl = 1.D-3*0.5*jhkl*dhkl*(Fhkl*lambda)**2*DW
      REAL(KIND(1.D0)) :: REFTAB(3,PCRYST_MAXREF) ! dhkl,jhkl,F2
      ! parameters needed for wavelength dependent total cross-section:
      ! capture: SIGA*lambda
      ! single-phonon: SIGS*lambda
      ! multi-phonon: SIGF*(1 - exp(-BMPH/lambda**2)
      ! incoherent: SIGI*DW_AVE
      ! Debye-Waller factor: DW = exp(-2*BT/dhkl**2)
      ! averaged DW factor: DW_AVE = (1 - exp(-8*BT/lambda**2))/(8*BT/lambda**2)
      REAL(KIND(1.D0)) :: SIGA, SIGI, SIGF, SIGS, BT, BMPH
      INTEGER :: PCRYST_NREF ! number of points, log-scale
      character*(128) :: REF_TABLE_NAME

      public READ_REF_TABLE,REF_TABLE_NAME,PCRYST_MAXREF,PCRYST_NREF,REFTAB
      public CLEAR_REF_TABLE, REF_TABLE_SSCATT, REF_TABLE_SABS, REPORT_SIGTOT

      contains

!---------------------------------------------------------
      subroutine REPORT_SIGTOT(IU, lmin, lmax, nlam)
!---------------------------------------------------------
      REAL(KIND(1.D0)), intent(in) :: lmin,lmax
      integer,intent(in) :: IU,nlam
      REAL(KIND(1.D0)) :: dlam,lam
      REAL(KIND(1.D0)) :: stot, sdif, sinc, sabs, ssph, smph
      REAL(KIND(1.D0)) :: row(7)
      integer:: i,j
1     format(G13.6,2x,G13.6)
2     format(a)
      dlam = (lmax-lmin)/(nlam-1)
      write(IU,2) '# Interaction cross-sections, 1/cm'
      write(IU,2) '# lambda, stot, sdif, sinc, sabs, ssph, smph'
      do i=1, nlam
        lam = lmin + dlam*i
        stot = 0.D0
        sdif = 0.D0
        sabs = 0.D0
        sinc = 0.D0
        ssph = 0.D0
        smph = 0.D0
        do j=1,PCRYST_NREF
          sdif = sdif + REF_TABLE_SSCATT(j,lam)
        enddo
        call SABS_COMPONENTS(lam, sinc, sabs, ssph, smph)
        stot = sdif+sinc+sabs+ssph+smph
        row = 10.0*(/0.1*lam, stot, sdif, sinc, sabs, ssph, smph/)
        write(IU,"(6(G13.6,a),G13.6)") (row(j),achar(9),j=1,6),row(7)
      enddo
      end subroutine REPORT_SIGTOT

!---------------------------------------------------------
      REAL(KIND(1.D0)) function REF_TABLE_SSCATT(iref,lambda)
! scattering cross section [1/mm] for given reflection
!---------------------------------------------------------
        integer, intent(in) :: iref
        REAL(KIND(1.D0)), intent(in) :: lambda
        REAL(KIND(1.D0)) :: SIG
        SIG = 0.D0
        if ((iref>0).and.(iref<=PCRYST_NREF)) then
          if (lambda<2.D0*REFTAB(1,iref)) then
            SIG=0.5D-3*REFTAB(2,iref)*REFTAB(1,iref)*(lambda*REFTAB(3,iref))**2
          endif
        endif
        REF_TABLE_SSCATT=SIG
      end function REF_TABLE_SSCATT


!---------------------------------------------------------
      subroutine SABS_COMPONENTS(lambda, SINC, SABS, SSPH, SMPH)
! total removal cross-section except diffraction (capture, incoherent, TDS)
!---------------------------------------------------------
        REAL(KIND(1.D0)), intent(in) :: lambda
        REAL(KIND(1.D0)), intent(out) :: SINC, SABS, SSPH, SMPH
        REAL(KIND(1.D0)) :: DW_AVE, Z
        Z = 8*BT/lambda**2
        DW_AVE = (1.D0 - exp(-Z))/Z
        SINC = SIGI*DW_AVE
        SABS = SIGA*lambda
        SSPH = SIGS*lambda
        SMPH = SIGF*(1.D0 - exp(-BMPH/lambda**2))
      end subroutine SABS_COMPONENTS

!---------------------------------------------------------
      REAL(KIND(1.D0)) function REF_TABLE_SABS(lambda)
! total removal cross-section except diffraction (capture, incoherent, TDS)
!---------------------------------------------------------
        REAL(KIND(1.D0)), intent(in) :: lambda
        REAL(KIND(1.D0)) :: SIG, DW_AVE, Z
        ! capture + single phonon
        SIG = (SIGA+SIGS)*lambda
        ! incoherent
        Z = 8*BT/lambda**2
        DW_AVE = (1.D0 - exp(-Z))/Z
        SIG = SIG + SIGI*DW_AVE
        ! multi-phonon
        SIG = SIG + SIGF*(1.D0 - exp(-BMPH/lambda**2))
        REF_TABLE_SABS=SIG
      end function REF_TABLE_SABS


!---------------------------------------------------------
      SUBROUTINE CLEAR_REF_TABLE
!---------------------------------------------------------
        PCRYST_NREF=0
        REFTAB=0.D0
        REF_TABLE_NAME='none'
        ! default parameters for iron:
        SIGA=0.0121
        SIGI=0.0034
        SIGF=0.0954
        SIGS=5.799D-4
        BT=0.088 ! [AA**2]
        BMPH=0.3076 ! [AA**2]

      end SUBROUTINE CLEAR_REF_TABLE


!---------------------------------------------------------------
      subroutine READ_REFS_PARAM(LINE,IERR)
! Read parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! IC      table index
! return: ID and index of the parameter from DEFREFPAR
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      character(*),intent(in) :: LINE
      integer, intent(out) :: IERR
      character*(16) :: ID
      real(kind(1.D0)) :: Z(1)
      integer :: IE,NA
      integer :: idx
1     format(a,' : ',3(1x,G12.6))
      ID=''
      IE=1
      idx=0
      call READ_ID_R8(LINE,ID,Z,1,NA,IE)
      if ((IE==0).and.(NA>0)) then
        select case(trim(ID))
        case('SIGA')
          SIGA=Z(1)
        case('SIGI')
          SIGI=Z(1)
        case('SIGF')
          SIGF=Z(1)
        case('SIGS')
          SIGS=Z(1)
        case('BT')
          BT=Z(1)
        case('BMPH')
          BMPH=Z(1)
        case default
          IE=1
        end select
        ! write(*,1) trim(ID),Z(1)
      endif
      IERR=IE
      end subroutine READ_REFS_PARAM



!---------------------------------------------------------------
      SUBROUTINE READ_REFS(IU,ilin,ierr)
! read reflections as a table with columns:
! dhkl [A], multiplicity, FHKL/CellVolume [fm A^3]
!---------------------------------------------------------------
      INTEGER,intent(in) :: IU
      INTEGER,intent(out) :: ilin,ierr
      CHARACTER*(1024) S
      integer :: j, ie
      real(kind(1.D0)) :: Z(3)
      PCRYST_NREF=0
      S=' '
      ierr=0
      do while (IERR.eq.0)
        call READ_TABLE_ROW(IU,ilin,S,ierr)
        J=LEN_TRIM(S)
        if ((ierr.eq.i_OK).and.(PCRYST_NREF.ge.PCRYST_MAXREF)) ierr=i_VALID
        if (ierr.eq.i_OK) then
          if (INDEX(S(1:J),'=')>0) then
            call READ_REFS_PARAM(S(1:J),ie)
            if (ie>=0) ierr=i_OK
          else
            Read(S(1:J),*,iostat=IERR,ERR=50) Z(1), Z(2), Z(3)
50          if (ierr.ne.0) then
              ierr=i_VALID
            else
              PCRYST_NREF=PCRYST_NREF+1
              REFTAB(1,PCRYST_NREF)=Z(1)
              REFTAB(2,PCRYST_NREF)=Z(2)
              REFTAB(3,PCRYST_NREF)=Z(3)*exp(-BT/Z(1)**2)
            endif
          endif
        endif
      enddo
      if (PCRYST_NREF.le.1) ierr=i_ERR
      end SUBROUTINE READ_REFS


!---------------------------------------------------------------
      SUBROUTINE READ_REF_TABLE(FNAME,ilin)
! Read flux table from a file
!---------------------------------------------------------------
      CHARACTER*(*),intent(in) :: FNAME
      integer,intent(out) :: ilin
      character(MAX_FNAME_LENGTH) :: FRES,FRESPATH
      INTEGER :: IU,IDCOMPARE
      INTEGER :: IERR
      CHARACTER*128 :: S,MSG
      character*32 :: CNUM1

      IERR=i_OK
      ilin=0
! empty string => clear table and exit
      if (LEN_TRIM(FNAME).eq.0) then
        call CLEAR_REF_TABLE
        Return
      endif
! nothing changed
      if (IDCOMPARE(REF_TABLE_NAME,FNAME).eq.0) then
        return
      endif
      ilin=1 ! to raise report mesage later
! from now, assume that a new table will be read
      call CLEAR_REF_TABLE

! no table => exit with clear table
      if (IDCOMPARE('none',trim(FNAME)).eq.0) then
        Return
      endif
! open file
      S=trim(FNAME)
      call OPENRESFILE(trim(S),' ',0,FRES,FRESPATH,IU,IERR)
      IF(IERR.NE.0) GOTO 100
      call READ_REFS(IU,ilin,ierr)
      select case (ierr)
      case(i_OK,i_EOF,i_VALID)
        goto 30
      case(i_ERR)
        goto 50
      case default
        goto 40
      end select
! close file, all is OK
30    CLOSE(IU)
      REF_TABLE_NAME=trim(FNAME)
      call INT2STR(PCRYST_NREF,CNUM1)
2     FORMAT('Loaded reflections table (',a,' reflections)')
      write(S,2) trim(CNUM1)
      call MSG_INFO(S,1)
      return

! error while reading
40    call CLEAR_REF_TABLE
4     FORMAT('Error ',I5,' while reading reflections table, line ',I5)
      write(S,4) IERR,ilin
      call MSG_ERROR('READ_REF_TABLE',S,0,1)
      return

! structure error, e.g. non-equidistant data
50    CLOSE(IU)
      call CLEAR_REF_TABLE
5     FORMAT('Error in table structure, ',a,': ',a)
      write(S,5) trim(FNAME),trim(MSG)
      call MSG_ERROR('READ_REF_TABLE',S,0,1)
      return

! error on file open
100   call MSG_ERROR('READ_REF_TABLE','Cannot open table with reflections: ['//trim(FNAME)//'].',0,1)
      call CLEAR_REF_TABLE
      return

      end SUBROUTINE READ_REF_TABLE

      end MODULE SAMPLE_POLY_TABLE
