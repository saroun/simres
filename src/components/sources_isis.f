!//////////////////////////////////////////////////////////////////////
!////  $Id: sources_isis.f,v 1.2 2017/09/07 18:52:28 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2016, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.2 $
!////     $Date: 2017/09/07 18:52:28 $
!////////////////////////////////////////////////////////////////////////
!////
!////  Handles moderator description using ISIS brightness tables for McStas
!////
!////////////////////////////////////////////////////////////////////////////////////

      MODULE SOURCES_ISIS
      use SIMATH
      use FILETOOLS
      implicit none
      save
      private


      TYPE SRCISIS; SEQUENCE
  ! Maxwell spectrum
        character(128) :: FNAME
        real(kind(1.D0)) :: LMIN,LMAX ! wavelength [A]
        real(kind(1.D0)) :: TMAX ! time range [us]
        real(kind(1.D0)) :: POW  ! power [kW]
        logical :: VALID
      end TYPE SRCISIS
      TYPE(SRCISIS) :: ISIS
! temporary bin data
      real(kind(1.D0)) :: EB(2),TB(2),TB0

! line read result, returned by PARSE_ISIS_LINE
      integer,parameter :: ires_none=0 ! ignore line
      integer,parameter :: ires_ebin=1 ! e-nbin read
      integer,parameter :: ires_tstart=2 ! starting flag for time section
      integer,parameter :: ires_tstop=3  ! stop flag for time section
      integer,parameter :: ires_tbin=4   ! t-bin read
      integer,parameter :: ires_err=-1   ! error
! read status flags, used by PARSE_ISIS_LINE
      integer,parameter :: sta_none=0  ! outside read section
      integer,parameter :: sta_ebin=1  ! new E-bin is defined
      integer,parameter :: sta_tbin=2  ! reading time bins
      real(kind(1.D0)),parameter  :: PCURR60=3.744905847e4*1.1879451/50 ! proton current / 1e10 at 60 kW


      logical :: dbg
      integer :: idbg


      public PARSE_ISIS_TABLE,ISIS_DEFAULT,SRCISIS,ISIS,READ_PARAM_ISIS

      contains
!---------------------------------------------------------------
      subroutine ISIS_DEFAULT
! set default values in ISIS data structure
!---------------------------------------------------------------
      ISIS%LMIN=3.D-1
      ISIS%LMAX=2.D1
      ISIS%TMAX=5.D2
      ISIS%POW=160.D0
      ISIS%FNAME=' '
      ISIS%VALID=.false.
      EB=0.D0
      TB=0.D0
      TB0=0.D0
      dbg=.false.
      idbg=0
      end subroutine ISIS_DEFAULT

!---------------------------------------------------------------
      subroutine READ_PARAM_ISIS(LINE,ID,IERR)
! read parameters for ISIS source
! Read source parameter value from a text line, assumed format is ID=value
! LINE    input text
! ID      parameter name
! return: ID and index of the parameter from NAMDEF123
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      character(*),intent(in) :: LINE
      character(*),intent(out) :: ID
      integer, intent(out) :: IERR
      character*(16) :: S2
      character*(128) :: FN
      real(kind(1.D0)) :: Z(3)
      integer :: IE,NA
      ID=''
      if (INDEX(LINE,'FNAME=').gt.0) then
        call READ_STR('FNAME',LINE,FN,IE)
        ID='FNAME'
        if ((IE.eq.0).and.(len_trim(FN).gt.0)) ISIS%FNAME=trim(FN)
      else
        call READ_ID_R8(LINE,S2,Z,3,NA,IE)
        if (IE==0) then
          call MKUPCASE(S2)
          select case(trim(S2))
              case('LMIN')
                ISIS%LMIN=Z(1)
              case('LMAX')
                ISIS%LMAX=Z(1)
              case('TMAX')
                ISIS%TMAX=Z(1)
              case('POW')
                ISIS%POW=Z(1)
              case default
                IE=1 ! unkown ID
          end select
        endif
        ID=trim(S2)
      endif
      IERR=IE
      end subroutine READ_PARAM_ISIS

!---------------------------------------------------------------
      SUBROUTINE PARSE_ISIS_LINE(LINE,STATE,Z,IRES)
! interpret one input line of the ISIS moderator file (McStas format)
!---------------------------------------------------------------
      character(*) :: LINE
      real(kind(1.D0)),intent(out) :: Z(2)
      INTEGER,intent(in) :: STATE
      INTEGER,intent(out) :: IRES
      integer :: LL,IS,IS1,LL1,NA
      character(128) :: S2,row
      IRES=ires_none
      LL=len_trim(LINE)
      select case(STATE)
      case(sta_none)
    ! search for "energy bin:"
        IS=index(LINE,'energy bin:')
        if (IS.gt.0) then
          S2=trim(LINE(IS+11:))
          IS1=index(S2,' to ')
          LL1=len_trim(S2)
          if (IS1>0) then
            row=S2(1:IS1)//' '//S2(IS1+3:LL1)
            call READ_ARRAY8(row,2,Z,NA)
            if (NA.gt.1) then
              IRES=ires_ebin
            endif
          else
            IRES=ires_err
          endif
          return
        endif
      case(sta_ebin)
    ! search for "time" = start of time bins
        IS=index(trim(LINE),'time')
        if (IS.eq.1) IRES=ires_tstart
      case(sta_tbin)
    ! search for "total" = end of time bins
        IS=index(trim(LINE),'total')
        if (IS.eq.1) then
          IRES=ires_tstop
          return
        else
    ! else read time bins
          call READ_ARRAY8(LINE,2,Z,NA)
          if (NA.gt.1) then
            IRES=ires_tbin
          else
            IRES=ires_err
          endif
        endif
      case default
        IRES=ires_none
      end select
      end SUBROUTINE PARSE_ISIS_LINE

!---------------------------------------------------------------
      SUBROUTINE ADD_BIN(FLUX,L0,DLAM,DT,NLAM,NT)
! fill flux table with the values of energy and time bins EB and TB
!---------------------------------------------------------------
      INTEGER,intent(in) :: NLAM,NT
      real(kind(1.D0)),intent(inout) :: FLUX(NLAM,NT)
      real(kind(1.D0)),intent(in) :: L0,DLAM,DT
      integer :: i,j,i1,i2,j1,j2
      real(kind(1.D0)) :: z,L1,L2,LL1,LL2,T1,T2,LMAX,Tfrac,Lfrac,tot,dpix
1     format(a,': ',6(G12.5))

      LMAX=L0*EXP((NLAM-1)*DLAM)
! calculate wavelength bin limits, convert from MeV to A
      L1=TWOPI*SQRT(HSQOV2M/EB(2)/1.D9)
      if (EB(1)>0.D0) then
        L2=TWOPI*SQRT(HSQOV2M/EB(1)/1.D9)
      else
        L2=L1+LMAX
      endif
! time bin limits, convert to us
      T1=TB0*1.D-2
      T2=TB(1)*1.D-2

      !dbg=((EB(1)>1e-8).and.(EB(1)<2.e-8))
      !dbg=(dbg.and.(T1>1.0))
      if (dbg) idbg=idbg+1
      !dbg=((idbg>0).and.(idbg<20))
      ! dbg=((L1>1).and.(L2<2))
      if (dbg) write(*,1) 'bin of energy,time',EB(1:2),TB(1:2)
      if (dbg) write(*,1) 'limits',L1,L2,T1,T2,DT

      if ((L2.le.L0).or.(L1.GE.LMAX).or.(T2.LE.0.D0).or.(T1.GE.(NT-1)*DT)) then
        return
      endif


! range of flux table cell indexes
! wavelength
      z=LOG(L1/L0)/DLAM
      i1=max(0,INT(z))
      z=LOG(L2/L0)/DLAM
      i2=min(NLAM-1,INT(z)+1)
! time
      j1=max(0,INT(T1/DT))
      j2=min(NT-1,INT(T2/DT)+1)
      if (dbg) write(*,1) 'bin range',i1,i2,j1,j2
      tot=0.D0
      dpix=(L2-L1)*(T2-T1)
      LL1=L0*EXP(i1*DLAM)
      LL2=L0*EXP((i1+1)*DLAM)
      if (dbg) write(*,1) 'rel. pix size',(LL2-LL1)/(L2-L1),DT/(T2-T1)
! scan table cells and add corresponding fractions of the source bin
      do i=i1,i2-1
        LL1=L0*EXP(i*DLAM)
        LL2=L0*EXP((i+1)*DLAM)
        Lfrac=min(L2,LL2)-max(L1,LL1)
        if (dbg) write(*,1) 'L-frac',i,Lfrac,Lfrac/(L2-L1)
        do j=j1,j2-1
          Tfrac=min(T2,(j+1)*DT)-max(T1,j*DT)
          if (dbg) write(*,1) 'T-frac',i,j,max(T1,j*DT),min(T2,(j+1)*DT),Tfrac/(T2-T1)
          ! add flux, convert to 1e10/s/cm2/sr/A
          FLUX(i+1,j+1)=FLUX(i+1,j+1)+TB(2)*Tfrac*Lfrac/dpix*PCURR60*ISIS%POW
          tot=tot+Tfrac*Lfrac/dpix
        enddo
      enddo
      if (dbg) write(*,1) 'fraction added',tot, TB(2)*tot*PCURR60*ISIS%POW
      end SUBROUTINE ADD_BIN

!---------------------------------------------------------------
      SUBROUTINE PARSE_ISIS_TABLE(FLUX,AVE,L0,DLAM,DT,NLAM,NT,IERR)
! parsing of a ISIS h.xxx moderator file
! OUTPUT:
! FLUX ... instant brightness, 2D lookup table
! AVE ... time averaged brightness, 1D lookup table
! L0 ... minimum wavelength
! DLAM ... logaritmic bin width for wavelength histogram [A]
! DT ... linear bin width for time histogram [us]
! NLAM ... number of wavelenth bins
! NT ... number of time bins
! leading index (rows) corresponds to wavelength
!---------------------------------------------------------------
      INTEGER,intent(in) :: NLAM,NT
      INTEGER,intent(out) :: IERR
      real(kind(1.D0)),intent(out) :: FLUX(NLAM,NT)
      real(kind(1.D0)),intent(out) :: AVE(NLAM)
      real(kind(1.D0)),intent(out) :: L0,DLAM,DT
      integer :: i,j,ilin,neb,ntb,STATE,ie,IU
      CHARACTER*512 :: S,MSG
      CHARACTER*1024 :: LINE
      character(MAX_FNAME_LENGTH) :: FRES,FRESPATH
      real(kind(1.D0)) :: Z(2),DUM,DL
      real(kind(1.D0)) :: y(NT)
1     format(a,': ',6(G12.5))
      DLAM=LOG(ISIS%LMAX/ISIS%LMIN)/(NLAM-1)
      L0=ISIS%LMIN
      DT=ISIS%TMAX/(NT-1)
      IERR=i_OK
      ilin=0
      ISIS%VALID=.false.
      do i=1,NLAM
        do j=1,NT
          FLUX(i,j)=0.D0
        enddo
      enddo
! empty string => clear table and exit
      if (LEN_TRIM(ISIS%FNAME).eq.0) then
        Return
      endif
! open file
      S=trim(ISIS%FNAME)
      call OPENRESFILE(trim(S),' ',0,FRES,FRESPATH,IU,IERR)
      IF(IERR.NE.i_OK) GOTO 100
! assume 1-line header
      neb=0
      ntb=0
      STATE=sta_none
      ! write(*,*) 'PARSE_ISIS_TABLE '//trim(S),' DT=',DT
      do while ((ierr.eq.i_VALID).or.(ierr.eq.i_OK))
        call READ_TABLE_ROW(IU,ilin,LINE,ierr)
        if (ierr.eq.i_OK) then
          call PARSE_ISIS_LINE(LINE,STATE,Z,IE)
          select case(IE)
          case(ires_ebin)
            ! new valid energy bin
            neb=neb+1
            EB=Z(1:2)
            TB(1)=0.D0
            STATE=sta_ebin
            !dbg=((EB(1)>1.8e-06).and.(EB(2)<2.2e-06))
            !if (dbg) write(*,1) 'EBIN ',EB
            ! dbg=(EB(1)>1.D-8)
            ! dbg=((EB(1)>1.8e-06).and.(EB(2)<2.2e-06))
          case(ires_tstart)
            ! start new time array
            if (STATE.eq.sta_ebin) STATE=sta_tbin
            ntb=0
          case(ires_tstop)
            ! end of time array
            STATE=sta_none
          case(ires_tbin)
            ! valid time bin read
            ntb=ntb+1
            TB0=TB(1)
            TB=Z(1:2)
            !dbg=((EB(1)>1.8e-06).and.(EB(2)<2.2e-06))
            !dbg=(dbg.and.(TB(2).gt.0.D0))
            call ADD_BIN(FLUX,L0,DLAM,DT,NLAM,NT)
            !dbg=.false.
          end select
        endif
      enddo

      select case (ierr)
      case(i_OK,i_EOF)
        goto 30
      case(i_ERR)
        goto 50
      case default
        goto 40
      end select

30    CLOSE(IU)


!     normalize to bin size, get integrals
      do i=1,NLAM
        DUM=0.D0
        DL=L0*(EXP((i+1)*DLAM)-EXP(i*DLAM))
        do j=1,NT
          FLUX(i,j) = FLUX(i,j)/DL/DT
          y(j)=FLUX(i,j)
          DUM=DUM+FLUX(i,j)
        enddo
        AVE(i)=INTEG1D(y,DT,NT)
        !write(*,1) 'AVE: ',L0*exp((i-1)*DLAM),L0*exp(i*DLAM),AVE(i),DUM*DT
      enddo
      !write(*,1) 'AVE integral:', INTEG1DLOG(AVE,L0,DLAM,NLAM)
      !write(*,*) 'PARSE_ISIS_TABLE OK, bins: ',neb,ntb,' lines: ',ilin
      ISIS%VALID=.true.
      return

! error while reading
4     FORMAT('Error ',I5,' while reading flux table, line ',I5)
40    write(S,4) IERR,ilin
      call MSG_ERROR('PARSE_ISIS_TABLE',S,0,1)
      return

! structure error, e.g. non-equidistant data
50    CLOSE(IU)
5     FORMAT('Error in table structure, ',a,': ',a)
      write(S,5) trim(ISIS%FNAME),trim(MSG)
      call MSG_ERROR('PARSE_ISIS_TABLE',S,0,1)
      return

! error on file open
100   call MSG_ERROR('PARSE_ISIS_TABLE','Cannot open table: ['//trim(ISIS%FNAME)//'].',0,1)
      return

      end SUBROUTINE PARSE_ISIS_TABLE


      end module SOURCES_ISIS