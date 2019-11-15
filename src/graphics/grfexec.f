!//////////////////////////////////////////////////////////////////////
!////  $Id: grfexec.f,v 1.59 2019/08/16 17:16:26 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.59 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Subroutines executing top-level graphics
!//////////////////////////////////////////////////////////////////////
      MODULE GRFEXEC
      use ENUMDEF
      use GRFPLOT
      use XMLINFO
      use FILETOOLS
      use COMMANDS
      use RESULTS
      USE RESOL1D
      use BEAM1D
      use BEAM2D
      use RESOL2D
      use COMPONENTS
      use EVENTMONITOR

      implicit NONE
! plot types
      integer,parameter :: pl_replot=-1
      integer,parameter :: pl_BEAM1D=0
      integer,parameter :: pl_BEAM2D=1
      integer,parameter :: pl_SCAN1D=2
      integer,parameter :: pl_SCAN2D=3
      integer,parameter :: pl_RES1D=4
      integer,parameter :: pl_RES2D=5
      integer,parameter :: pl_DET1D=6
      integer,parameter :: pl_DET2D=7
      integer,parameter :: pl_GAUGE2D=8

      contains

!----------------------------------------------------------------------
      SUBROUTINE DOPLOT(REPLOT)
! Execue PLOT command
!----------------------------------------------------------------------
      logical,intent(in) :: REPLOT
      INTEGER :: ICOM,TOPRINT,ENABLED
      ENABLED=NINT(getCmdParam('PLOT','ENABLED'))
      TOPRINT=NINT(getCmdParam('PLOT','PRINT'))
      if (REPLOT) TOPRINT=0
      if (.not.REPLOT) then
        ICOM=NINT(getCmdParam('PLOT','TYPE'))
        !write(*,*) 'DOPLOT ICOM=',ICOM
! define new graphic output
        if (TOPRINT.LE.0) then
          call ClrAllGraphs
          select case(ICOM)
          CASE(pl_BEAM1D)
            CALL VIEWBEAM1D('BEAM1D')
          CASE(pl_DET1D)
            CALL VIEWBEAM1D('DET1D')
          CASE(pl_BEAM2D)
            CALL VIEWBEAM2D('BEAM2D')
          CASE(pl_DET2D)
            CALL VIEWBEAM2D('DET2D')
          CASE(pl_SCAN1D)
            CALL VIEWSCAN1D
          CASE(pl_SCAN2D)
            CALL VIEWSCAN2D
          CASE(pl_RES1D)
            CALL VIEWRES1D
          CASE(pl_RES2D)
            CALL VIEWRES2D
          CASE(pl_GAUGE2D)
            CALL VIEWBEAM2D('GAUGE2D')
          end select
        endif
      endif
      IF (ENABLED.le.0) THEN
        WRITE(SOUT,*) 'graphics output is disabled'
        RETURN
      ENDIF
! open PGPLOT device
      if (INITGRF(TOPRINT)) then
    ! replot all
        CALL PlotAllGraphs
    ! close PGPLOT device
        CALL CLOSEGRF(TOPRINT)
      else
        call SELGRFDEV('/NULL')
      endif
      END SUBROUTINE DOPLOT


!----------------------------------------------------------------------
      SUBROUTINE VIEWRES1D
! Prepare plot of 1D resolution function in given direction (integrate over others)
! Plot parameters are taken fom command object RES1D
!----------------------------------------------------------------------
      character*20 cstat(4),cfit(3)
      real :: par(3),dw,YFIT(SPCM)
      REAL(KIND(1.D0)) :: X0,DX,suma,center,fwhm,wspread,ymax
      integer :: ip,ig,i,IVAR,ER
      logical :: auto,gfit

14    FORMAT('fwhm: ',G10.3,' spread: ',G10.3,' center: ',G10.3)
21    FORMAT('Imax: ',G10.4)
22    FORMAT('sum: ',G10.4)
23    FORMAT('fwhm: ',G10.4)
24    FORMAT('spread: ',G10.4)
25    FORMAT('center: ',G10.4)
26    FORMAT('mean: ',G10.4)

      ip=1  ! port index
  ! set viewport screen
      GRFPORT(ip)%DC=(/0.1,0.9,0.3,0.8/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
! GRAPH DATA
      ig=getNextGraph(grtyp_1D)  ! graph index
  ! number of points
      GRF1D(ig)%NP=MIN(SPCM,NINT(getCmdParam('RES1D','NP')))
  ! variable index
      IVAR=MAX(1,NINT(getCmdParam('RES1D','X'))+1)
      IVAR=MIN(IVAR,QE_MAX)
  ! center
      X0=getCmdParam('RES1D','X0')
  ! range
      DX=getCmdParam('RES1D','DX')
  ! autoscale
      auto=(getCmdParam('RES1D','XAUTO').gt.0.D0)
      if (DX.le.0) auto=.true.
  ! logscale
      if (getCmdParam('RES1D','YLOG').gt.0.D0) GRFPORT(ip)%LOGSCALE(2)=10.0
  ! fit gauss
      gfit=(getCmdParam('RES1D','GFIT').gt.0.D0)
  ! errors
      ER=NINT(getCmdParam('RES1D','ERR'))
  ! captions
      !call FINDSTRPAR(QE_NAMES,':',IVAR,IS,IL)
      !call FINDSTRPAR(QE_UNITS,':',IVAR,IS1,IL1)
      !if (IL1.gt.0) then
      !  GRF1D(ig)%X_CAP=QE_NAMES(IS:IS+IL-1)//' ['//QE_UNITS(IS1:IS1+IL1-1)//']'
      !else
      !  GRF1D(ig)%X_CAP=QE_NAMES(IS:IS+IL-1)
      !endif
      call QVAR_GETNAME(IVAR,.true.,GRF1D(ig)%X_CAP)
      GRF1D(ig)%Y_CAP='Intensity'
      GRFPORT(ip)%HEAD='Resolution function'

  ! apply autoscale
      if (auto) call GETQRANGE(IVAR,X0,DX)
  ! get curve data
      call DORESOL1D(IVAR,X0,DX,GRF1D(ig)%NP,RES1D)
      GRF1D(ig)%NP=MIN(RES1D%N,MAXPNTS)
      ymax=0.D0
      do I=1,GRF1D(ig)%NP
        GRF1D(ig)%Y(i)=RES1D%Y(I)
        GRF1D(ig)%DY(i)=RES1D%DY(I)
        GRF1D(ig)%X(i)=RES1D%X(I)
        if (RES1D%Y(I)>ymax) ymax=RES1D%Y(I)
      enddo
      GRF1D(ig)%ER=ER
  ! get peak parameters (optional)
      call PEAKPARAM(RES1D,suma,center,fwhm,wspread)
  ! format text with pek parameters
      WRITE(cstat(1),22) suma
      WRITE(cstat(2),23) fwhm
      WRITE(cstat(3),24) wspread
      WRITE(cstat(4),26) center
  ! set styles
      GRF1D(ig)%COLOR=3
      GRF1D(ig)%POINT=4
      GRF1D(ig)%LINE=1
      GRF1D(ig)%PSIZE=MAX(0.3,1.5*dw)
  ! assign viewport
      call AssignGraph(ig,ip,grtyp_1D)
  ! set y-scale
      call SetYscale(ip)

  ! fit Gaussian (optional)
      gfit=(gfit.and.(GRFPORT(ip)%WC(4).gt.0.D0))
      if (gfit) then
        ig=getNextGraph(grtyp_1D)  ! new graph index
        ! PAR(1)=GRFPORT(ip)%WC(4)
        PAR(1)=ymax
        PAR(2)=center
        !PAR(3)=wspread
        PAR(3)=fwhm
        CALL FITGAUSS(PAR,RES1D%X,YFIT,RES1D%N)
        RES1D%Y(1:RES1D%N)=YFIT(1:RES1D%N)
        WRITE(cfit(1),21) par(1)
        WRITE(cfit(2),25) par(2)
        WRITE(cfit(3),24) par(3)
        GRF1D(ig)=GRF1D(ig-1)      ! copy of previous
        GRF1D(ig)%Y(1:RES1D%N)=YFIT(1:RES1D%N)
        GRF1D(ig)%COLOR=2
        GRF1D(ig)%POINT=0
        GRF1D(ig)%LINE=1
        GRF1D(ig)%ER=0
        call AssignGraph(ig,ip,grtyp_1D)
      endif

! VIEWPORT
  ! set physical x-scale
      GRFPORT(ip)%WC(1:2)=(/x0-dx/2,x0+dx/2/)
  ! reference points
      GRFPORT(ip)%REF=0.0
  ! viewport titles etc..
      GRFPORT(ip)%XTIT=GRF1D(ig)%X_CAP
      GRFPORT(ip)%YTIT=GRF1D(ig)%Y_CAP
      call getCmdParamS('RES1D',1,'COM',GRFPORT(ip)%comment)
      GRFPORT(ip)%legend=' '
      WRITE(GRFPORT(ip)%legend,14) fwhm,wspread,center

  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))
! LEGENDS
  ! peak parameters
      ig=getNextGraph(grtyp_txt)
      GRFTXT(ig)%COLOR=1
      GRFTXT(ig)%DIST=10.0
      GRFTXT(ig)%SHIFT=-0.05
      GRFTXT(ig)%SIDE='B'
      GRFTXT(ig)%LWIDTH=2
      GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
      GRFTXT(ig)%TEXT=' '
      call APPENDPAR(GRFTXT(ig)%TEXT,'|','statistics:')
      do i=1,4
        call APPENDPAR(GRFTXT(ig)%TEXT,'|',cstat(i))
      enddo
      call AssignGraph(ig,ip,grtyp_txt)
  ! fit result
      if (gfit) then
        ig=getNextGraph(grtyp_txt)
        GRFTXT(ig)=GRFTXT(ig-1)
        GRFTXT(ig)%COLOR=2
        GRFTXT(ig)%SHIFT=0.5
        GRFTXT(ig)%TEXT=' '
        call APPENDPAR(GRFTXT(ig)%TEXT,'|','gaussian fit:')
        do i=1,3
          call APPENDPAR(GRFTXT(ig)%TEXT,'|',cfit(i))
        enddo
        call AssignGraph(ig,ip,grtyp_txt)
      endif

! print additional info on console
  !    WRITE(sout,14) fwhm,wspread,center
      END subroutine VIEWRES1D

!----------------------------------------------------------------------
      SUBROUTINE VIEWBEAM1D(TYP)
! Prepare plot of 1D beam projection in given direction (integrate over others)
! Plot parameters are taken fom command object BEAM1D
!----------------------------------------------------------------------
      use OPTIMIZATION
      character(*) :: TYP
      character*20 cstat(4),cfit(3)
      real :: par(3),dw,YFIT(SPCM),MPAR(3*MAXPEAKS),dMPAR(3*MAXPEAKS)
      REAL(KIND(1.D0)) :: X0,DX,ADX,suma,center,fwhm,wspread,CHISQR
      integer :: ip,ig,i,IVAR,ER,what,np,filt
      logical :: auto,autoctr,gfit,app
      real :: ARR2(2)
      character*128 :: COMM
      TYPE(PDETECTOR) :: DET
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
1     format(a,6(G13.5,1x))
14    FORMAT('fwhm: ',G10.3,' spread: ',G10.3,' center: ',G10.3)
21    FORMAT('Imax: ',G10.4)
22    FORMAT('sum: ',G10.4)
23    FORMAT('fwhm: ',G10.4)
24    FORMAT('spread: ',G10.4)
25    FORMAT('center: ',G10.4)
26    FORMAT('mean: ',G10.4)
27    FORMAT('chi2: ',G10.4)
28    FORMAT('peaks: ',I3)
29    FORMAT('test version')

      select case(trim(TYP))
      case('BEAM1D')
        what=show_monitor
      case('DET1D')
        what=show_detector
        if (.not.GetDetector(DET)) return
      end select
      filt=-1 ! default event filter
      ip=1  ! port index
! set viewport screen
      GRFPORT(ip)%DC=(/0.1,0.9,0.3,0.8/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
! GRAPH DATA
      ig=getNextGraph(grtyp_1D)  ! graph index
  ! number of points
      GRF1D(ig)%NP=MIN(SPCM,NINT(getCmdParam(TYP,'NP')))
  ! variable index
      IVAR=MAX(1,NINT(getCmdParam(TYP,'X'))+1)
      select case(what)
      case(show_monitor)
        IVAR=MIN(IVAR,RK_MAX)
      case(show_detector)
        IVAR=MIN(IVAR,DT_MAX)
      end select
  ! center
      X0=getCmdParam(TYP,'X0')
  ! range
      DX=getCmdParam(TYP,'DX')
  ! autoscale
      auto=(getCmdParam(TYP,'XAUTO').gt.0.D0)
      autoctr=(getCmdParam(TYP,'CAUTO').gt.0.D0)
      if (DX.le.0) auto=.true.
  ! logscale
      if (getCmdParam(TYP,'YLOG').gt.0.D0) GRFPORT(ip)%LOGSCALE(2)=10.0
  ! fit gauss
      gfit=(getCmdParam(TYP,'GFIT').gt.0.D0)
  ! errors
      ER=NINT(getCmdParam(TYP,'ERR'))
  ! append
      APP=(getCmdParam(TYP,'APPEND').gt.0.D0)
  ! coordinates
      if (what==show_monitor) then
        i=NINT(getCmdParam(TYP,'COORD'))
      else
        i=coord_axis
      endif
  ! filter
      filt=NINT(getCmdParam(TYP,'FILT'))
      call BMONITORS_SETCOORD(i)
      modulation_dt=0.D0
      modulation_d0=0.D0
      call getCmdParamS(TYP,1,'COM',COMM)
      GRFPORT(ip)%comment=trim(COMM)
      if (INDEX(COMM,"mod=")==1) then
        call STR2ARRAY4(COMM(5:),':',ARR2,2,i)
        if (i>1) then
          modulation_dt=ARR2(1)
          modulation_d0=ARR2(2)
        endif
        write(*,1) 'modulation dt,d0: ',modulation_dt,modulation_d0
      endif

  ! captions
      !call FINDSTRPAR(RK_NAMES,':',IVAR,IS,IL)
      !call FINDSTRPAR(RK_UNITS,':',IVAR,IS1,IL1)
      !if (IL1.gt.0) then
      !  GRF1D(ig)%X_CAP=RK_NAMES(IS:IS+IL-1)//' ['//RK_UNITS(IS1:IS1+IL1-1)//']'
      !else
      !  GRF1D(ig)%X_CAP=RK_NAMES(IS:IS+IL-1)
      !endif
      if (PULSED_SRC) then
        GRF1D(ig)%Y_CAP='Intensity [n/pulse]'
      else
        GRF1D(ig)%Y_CAP='Intensity [n/s]'
      endif
      select case(what)
      case(show_monitor)
        call VAR_GETNAME(IVAR,.true.,GRF1D(ig)%X_CAP)
        GRFPORT(ip)%HEAD='Beam profile'
        if (auto.or.autoctr) then
          call GETRANGE1D(IVAR,filt,X0,ADX)
          if (auto) DX=ADX
        endif
      case(show_detector)
        call DVAR_GETNAME(IVAR,.true.,GRF1D(ig)%X_CAP)
        GRFPORT(ip)%HEAD='Detector data'
        if (auto.or.autoctr) then
          call DETRANGE1D(IVAR,DET%X,X0,ADX)
          if (auto) DX=ADX
        endif
      end select

  ! get curve data
      RES1D%APPEND=app
      RES1D%FILT=filt
      !write(*,1) 'VIEWBEAM1D before DOBEAM1D'
      call DOBEAM1D(IVAR,what,(/X0-0.5*DX, X0+0.5*DX/),GRF1D(ig)%NP,RES1D)
      !write(*,1) 'VIEWBEAM1D after DOBEAM1D'
      GRF1D(ig)%NP=MIN(RES1D%N,MAXPNTS)
      !write(*,1) 'VIEWBEAM1D ',GRF1D(ig)%NP,RES1D%N
      do I=1,GRF1D(ig)%NP
        GRF1D(ig)%Y(i)=RES1D%Y(I)
        GRF1D(ig)%DY(i)=RES1D%DY(I)
        GRF1D(ig)%X(i)=RES1D%X(I)
        !write(*,1) 'VIEWBEAM1D ',GRF1D(ig)%X(i),GRF1D(ig)%Y(i)
      enddo
       !write(*,1) 'VIEWBEAM1D range ',X0,DX
      GRF1D(ig)%ER=ER

  ! get peak parameters (optional)
      if ((RES1D%NPK==0).or.(RES1D%PKMAX==RES1D%PKMIN)) then
        call PEAKPARAM(RES1D,suma,center,fwhm,wspread)
  ! format text with pek parameters
        WRITE(cstat(1),22) suma
        WRITE(cstat(2),23) fwhm
        WRITE(cstat(3),24) wspread
        WRITE(cstat(4),26) center
      endif
  ! set styles
      GRF1D(ig)%COLOR=3
      GRF1D(ig)%POINT=4
      GRF1D(ig)%LINE=1
      GRF1D(ig)%PSIZE=MAX(0.3,1.5*dw)
  ! assign viewport
      call AssignGraph(ig,ip,grtyp_1D)
  ! set y-scale
      call SetYscale(ip)
  ! fit Gaussian (optional)
      gfit=(gfit.and.(GRFPORT(ip)%WC(4).gt.0.D0))
      if (gfit) then
        !write(*,*) 'VIEWBEAM1D NPK=',RES1D%NPK
        ig=getNextGraph(grtyp_1D)  ! new graph index
        GRF1D(ig)=GRF1D(ig-1)      ! copy of previous
        if ((RES1D%NPK==0).or.(RES1D%PKMAX==RES1D%PKMIN)) then
          PAR(1)=GRFPORT(ip)%WC(4)
          PAR(2)=center
          PAR(3)=wspread
          CALL FITGAUSS(PAR,RES1D%X,YFIT,RES1D%N)
          !RES1D%Y(1:RES1D%N)=YFIT(1:RES1D%N)
          WRITE(cfit(1),21) par(1)
          WRITE(cfit(2),25) par(2)
          WRITE(cfit(3),24) par(3)
          GRF1D(ig)%NPEAKS=1
          GRF1D(ig)%PKPAR(1:3)=PAR(1:3)
        else
          CHISQR=FITMGAUSS2(IVAR,MPAR,dMPAR,RES1D%X,YFIT,RES1D%N,NP)
          !RES1D%Y(1:RES1D%N)=YFIT(1:RES1D%N)
          WRITE(cfit(1),27)  CHISQR
          WRITE(cfit(2),28)  RES1D%PKMAX-RES1D%PKMIN+1
          WRITE(cfit(3),29)
          GRF1D(ig)%NPEAKS=RES1D%PKMAX-RES1D%PKMIN+1
          GRF1D(ig)%PKPAR(1:3*GRF1D(ig)%NPEAKS)=MPAR(1:3*GRF1D(ig)%NPEAKS)
          GRF1D(ig)%PKERR(1:3*GRF1D(ig)%NPEAKS)=dMPAR(1:3*GRF1D(ig)%NPEAKS)
          !write(*,*) 'GRFDATA:  ',ig,GRF1D(ig)%NPEAKS,GRF1D(ig)%PKPAR(1:3)
        endif
        GRF1D(ig)%Y(1:RES1D%N)=YFIT(1:RES1D%N)
        GRF1D(ig)%COLOR=2
        GRF1D(ig)%POINT=0
        GRF1D(ig)%LINE=1
        GRF1D(ig)%ER=0
        call AssignGraph(ig,ip,grtyp_1D)
      endif

! VIEWPORT
  ! set physical x-scale
      GRFPORT(ip)%WC(1:2)=(/x0-dx/2,x0+dx/2/)
  ! reference points
      GRFPORT(ip)%REF=0.0
  ! viewport titles etc..
      GRFPORT(ip)%XTIT=GRF1D(ig)%X_CAP
      GRFPORT(ip)%YTIT=GRF1D(ig)%Y_CAP
      GRFPORT(ip)%legend=' '
      WRITE(GRFPORT(ip)%legend,14) fwhm,wspread,center

  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))
! LEGENDS
  ! peak parameters
      if ((RES1D%NPK==0).or.(RES1D%PKMAX==RES1D%PKMIN)) then
        ig=getNextGraph(grtyp_txt)
        GRFTXT(ig)%COLOR=1
        GRFTXT(ig)%DIST=8.0
        GRFTXT(ig)%SHIFT=-0.05
        GRFTXT(ig)%SIDE='B'
        GRFTXT(ig)%LWIDTH=2
        GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
        GRFTXT(ig)%TEXT=' '
        call APPENDPAR(GRFTXT(ig)%TEXT,'|','statistics:')
        do i=1,4
          if (len_trim(cstat(i))>0) call APPENDPAR(GRFTXT(ig)%TEXT,'|',cstat(i))
        enddo
        call AssignGraph(ig,ip,grtyp_txt)
      endif
  ! fit result
      if (gfit) then
        ig=getNextGraph(grtyp_txt)
        GRFTXT(ig)%DIST=8.0
        GRFTXT(ig)%SIDE='B'
        GRFTXT(ig)%LWIDTH=2
        GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
        GRFTXT(ig)%COLOR=2
        GRFTXT(ig)%SHIFT=0.5
        GRFTXT(ig)%TEXT=' '
        call APPENDPAR(GRFTXT(ig)%TEXT,'|','gaussian fit:')
        do i=1,3
          if (len_trim(cfit(i))>0) call APPENDPAR(GRFTXT(ig)%TEXT,'|',cfit(i))
        enddo
        call AssignGraph(ig,ip,grtyp_txt)
      endif
! TOF
      if (TOF_MODE) then
        call LABEL_TOF(ip)
      endif
! monitor ID
      if (what==show_monitor) then
        call LABEL_MONITOR(ip)
      endif
! print additional info on console
  !    WRITE(sout,14) fwhm,wspread,center
      END subroutine VIEWBEAM1D

!----------------------------------------------------------------------
      SUBROUTINE VIEWBEAM2D(TYP)
! Prepare plot of 2D beam projection in given direction (integrate over others)
! Plot parameters are taken from command object BEAM2D
!----------------------------------------------------------------------
      character(*) :: TYP
      REAL(KIND(1.D0)) :: X0,DX,Y0,DY
      integer :: ip,ig,i,j,IX,IY,IZ,what,IU,filt
      real :: dw
      logical :: auto,APP
      TYPE(PDETECTOR) :: DET
      character*128 :: COMM
      character*32 :: SS
      character*512 :: FN
1     format(a,' ',3(G12.5,2x))

      IU=-1
      FN=' '
      filt=-1
      IZ=0
      select case(trim(TYP))
      case('BEAM2D')
        what=show_monitor
      case('DET2D')
        what=show_detector
        if (.not.GetDetector(DET)) then
          call MSG_WARN('Instrument has no detector.',1)
          return
        endif
      case('GAUGE2D')
        what=show_gauge
        if (.not.GetDetector(DET)) then
          call MSG_WARN('Instrument has no detector.',1)
          return
        endif
      case default
        return
      end select
      ! write(*,*) 'VIEWBEAM2D ',trim(TYP),what

      ip=1  ! port index
! GRAPH DATA
      ig=getNextGraph(grtyp_2D)  ! graph index
  ! number of pixels
      i=MIN(MAXPIX,NINT(getCmdParam(TYP,'NX')))
      GRF2D(ig)%NX=max(16,i)
      i=MIN(MAXPIX,NINT(getCmdParam(TYP,'NY')))
      GRF2D(ig)%NY=max(16,i)
  ! variable index
      IX=MAX(1,NINT(getCmdParam(TYP,'X'))+1)
      IY=MAX(1,NINT(getCmdParam(TYP,'Y'))+1)
      select case(what)
      case(show_monitor)
        IX=MIN(IX,RK_MAX)
        IY=MIN(IY,RK_MAX)
        IZ=0
      case(show_detector)
        IX=MIN(IX,DT_MAX)
        IY=MIN(IY,DT_MAX)
        IZ=0
      case(show_gauge)
      ! NOTE: list of gauge variables = intensity, followed by the list of detector variables
        IX=MIN(IX,RK_MAX)
        IY=MIN(IY,RK_MAX)
        IZ=MAX(1,NINT(getCmdParam(TYP,'Z'))+1)
        IZ=MIN(IZ,DT_MAX+1)
      end select
  ! center
      X0=getCmdParam(TYP,'X0')
      Y0=getCmdParam(TYP,'Y0')
  ! range
      DX=getCmdParam(TYP,'DX')
      DY=getCmdParam(TYP,'DY')
  ! autoscale
      auto=(getCmdParam(TYP,'AUTO').gt.0.D0)
      if (DX.le.0) auto=.true.
      if (DY.le.0) auto=.true.
  ! logscale
      if (getCmdParam(TYP,'ZLOG').gt.0.D0) GRFPORT(ip)%LOGSCALE(3)=10.0
  ! minimum-z method
      GRFPORT(ip)%MINZ=NINT(getCmdParam(TYP,'MINZ'))
  ! append
      APP=(getCmdParam(TYP,'APPEND').gt.0.D0)
  ! coordinates
      select case(what)
      case(show_monitor,show_gauge)
        i=NINT(getCmdParam(TYP,'COORD'))
      case(show_detector)
        i=coord_axis
      end select
      call BMONITORS_SETCOORD(i)
  ! filter
      filt=NINT(getCmdParam(TYP,'FILT'))
  ! captions
      !call FINDSTRPAR(RK_NAMES,':',IX,IS,IL)
      !call FINDSTRPAR(RK_UNITS,':',IX,IS1,IL1)
      !if (IL1.gt.0) then
      !  GRFPORT(ip)%XTIT=RK_NAMES(IS:IS+IL-1)//' ['//RK_UNITS(IS1:IS1+IL1-1)//']'
      !else
      !  GRFPORT(ip)%XTIT=RK_NAMES(IS:IS+IL-1)
      !endif
  ! comment
      call getCmdParamS(TYP,1,'COM',COMM)
  ! debug option: dump file name is given as a comment
      if (INDEX(COMM,"dump=")>0) then
        FN=trim(OUTPATH)//'/'//trim(COMM(6:))
        IU=OPENFILEUNIT(trim(FN),.false.)
        if (IU>0) then
          write(*,*) 'Dump events to '//trim(FN)
          write(IU,1) '# Dump of events from 2D plot'
        endif
      endif

      select case(what)
      case(show_monitor)
        call VAR_GETNAME(IX,.true.,GRFPORT(ip)%XTIT)
        call VAR_GETNAME(IY,.true.,GRFPORT(ip)%YTIT)
        GRFPORT(ip)%HEAD='Beam profile'
        if (auto) call GETRANGE2D(BMONITOR_REF,IX,IY,filt,X0,DX,Y0,DY)
      case(show_detector)
        call DVAR_GETNAME(IX,.true.,GRFPORT(ip)%XTIT)
        call DVAR_GETNAME(IY,.true.,GRFPORT(ip)%YTIT)
        GRFPORT(ip)%HEAD='Detector data map'
        if (auto) call DETRANGE2D(IX,IY,DET%X,X0,DX,Y0,DY)
      case(show_gauge)
        call VAR_GETNAME(IX,.true.,GRFPORT(ip)%XTIT)
        call VAR_GETNAME(IY,.true.,GRFPORT(ip)%YTIT)
        if (IZ==1) then
          SS='Intensity'
        else
          call DVAR_GETNAME(IZ-1,.true.,SS)
        endif
        GRFPORT(ip)%HEAD='Gauge map, '//trim(SS)
        if (auto) call GETRANGE2D(BMONITOR_SEC,IX,IY,filt,X0,DX,Y0,DY)
      end select
      if (IU>0) then
        write(IU,1) '# "'//trim(GRFPORT(ip)%XTIT)//'" "'//trim(GRFPORT(ip)%YTIT)//'" "weight"'
        write(IU,1) '# center:',X0,Y0
        write(IU,1) '# range:',DX,DY
      endif

  ! reference points
      GRFPORT(ip)%REF=0.0
  ! get curve data
      if (.not.APP) call CLEAR_RES2D(RES2D)
      RES2D%APPEND=APP
      RES2D%FILT=filt
      if (what==show_gauge) then
        call DOGAUGE2D(IX,IY,IZ-1,X0,Y0,DX,DY,GRF2D(ig)%NX,GRF2D(ig)%NY,RES2D,IU)
      else
        call DOBEAM2D(IX,IY,what,X0,Y0,DX,DY,GRF2D(ig)%NX,GRF2D(ig)%NY,RES2D,IU)
      endif
      if (IU>0) close(IU)
      GRF2D(ig)%NX=MIN(RES2D%NX,MAXPIX)
      GRF2D(ig)%NY=MIN(RES2D%NY,MAXPIX)
      do J=1,GRF2D(ig)%NY
        do I=1,GRF2D(ig)%NX
          GRF2D(ig)%Z(i,j)=RES2D%Z(i,j)
        enddo
      enddo
  ! assign viewport
      call AssignGraph(ig,ip,grtyp_2D)
! VIEWPORT
  ! set viewport screen
      GRFPORT(ip)%DC=(/0.15,0.93,0.31,0.89/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
  ! set physical scale
      GRFPORT(ip)%WC=(/X0-DX/2,X0+DX/2,Y0-DY/2,Y0+DY/2/)
  ! viewport titles etc..
      GRFPORT(ip)%comment=trim(COMM)
      !call getCmdParamS(TYP,1,'COM',GRFPORT(ip)%comment)
      GRFPORT(ip)%legend=' '
  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))


  ! TOF
      if (TOF_MODE) then
        call LABEL_TOF(ip)
      endif

! monitor ID
      if (what==show_monitor) then
        call LABEL_MONITOR(ip)
      endif
      END subroutine VIEWBEAM2D

!----------------------------------------------------------------------
      SUBROUTINE LABEL_MONITOR(ip)
! add monitor name label
!----------------------------------------------------------------------
      integer,intent(in) :: ip
      integer :: ig
      character*32 monid
        ig=getNextGraph(grtyp_txt)
        GRFTXT(ig)%COLOR=1
        GRFTXT(ig)%DIST=1
        GRFTXT(ig)%SHIFT=0.7
        GRFTXT(ig)%SIDE='T'
        GRFTXT(ig)%LWIDTH=2
        GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
        call NSTORE_GETID(BMONITOR_REF,monid)
        GRFTXT(ig)%TEXT='monitor: '//trim(monid)
        call AssignGraph(ig,ip,grtyp_txt)
      end SUBROUTINE LABEL_MONITOR


!----------------------------------------------------------------------
      SUBROUTINE LABEL_TOF(ip)
! add monitor name label
!----------------------------------------------------------------------
      integer,intent(in) :: ip
      integer :: ig
      character*32 monid
        ig=getNextGraph(grtyp_txt)
        GRFTXT(ig)%COLOR=1
        GRFTXT(ig)%DIST=1
        GRFTXT(ig)%SHIFT=0.0
        GRFTXT(ig)%SIDE='T'
        GRFTXT(ig)%LWIDTH=2
        GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
        call FLOAT2STR(TOF_ZERO/HOVM,monid)
        GRFTXT(ig)%TEXT='TOF_0='//trim(monid)
        call AssignGraph(ig,ip,grtyp_txt)
      end SUBROUTINE LABEL_TOF


!----------------------------------------------------------------------
      SUBROUTINE VIEWRES2D
! Prepare plot of 2D beam projection in given direction (integrate over others)
! Plot parameters are taken from command object BEAM2D
!----------------------------------------------------------------------
      REAL(KIND(1.D0)) :: X0,DX,Y0,DY
      integer :: ip,ig,i,j,IX,IY
      real :: dw
      logical :: auto
      ip=1  ! port index
! GRAPH DATA
      ig=getNextGraph(grtyp_2D)  ! graph index
  ! number of pixels
      i=MIN(MAXPIX,NINT(getCmdParam('RES2D','NX')))
      GRF2D(ig)%NX=max(16,i)
      i=MIN(MAXPIX,NINT(getCmdParam('RES2D','NY')))
      GRF2D(ig)%NY=max(16,i)
  ! variable index
      IX=MAX(1,NINT(getCmdParam('RES2D','X'))+1)
      IX=MIN(IX,QE_MAX)
      IY=MAX(1,NINT(getCmdParam('RES2D','Y'))+1)
      IY=MIN(IY,QE_MAX)
  ! center
      X0=getCmdParam('RES2D','X0')
      Y0=getCmdParam('RES2D','Y0')
  ! range
      DX=getCmdParam('RES2D','DX')
      DY=getCmdParam('RES2D','DY')
  ! autoscale
      auto=(getCmdParam('RES2D','AUTO').gt.0.D0)
      if (DX.le.0) auto=.true.
      if (DY.le.0) auto=.true.
  ! logscale
      if (getCmdParam('RES2D','ZLOG').gt.0.D0) GRFPORT(ip)%LOGSCALE(3)=10.0
  ! captions
      !call FINDSTRPAR(QE_NAMES,':',IX,IS,IL)
      !call FINDSTRPAR(QE_UNITS,':',IX,IS1,IL1)
      !if (IL1.gt.0) then
      !  GRFPORT(ip)%XTIT=QE_NAMES(IS:IS+IL-1)//' ['//QE_UNITS(IS1:IS1+IL1-1)//']'
      !else
      !  GRFPORT(ip)%XTIT=QE_NAMES(IS:IS+IL-1)
      !endif
      !call FINDSTRPAR(QE_NAMES,':',IY,IS,IL)
      !call FINDSTRPAR(QE_UNITS,':',IY,IS1,IL1)
      !if (IL1.gt.0) then
      !  GRFPORT(ip)%YTIT=QE_NAMES(IS:IS+IL-1)//' ['//QE_UNITS(IS1:IS1+IL1-1)//']'
      !else
      !  GRFPORT(ip)%YTIT=QE_NAMES(IS:IS+IL-1)
      !endif
      call QVAR_GETNAME(IX,.true.,GRFPORT(ip)%XTIT)
      call QVAR_GETNAME(IY,.true.,GRFPORT(ip)%YTIT)
      GRFPORT(ip)%HEAD='Resolution funtion'
  ! reference points
      GRFPORT(ip)%REF=0.0
  ! apply autoscale
      if (auto) call GETQRANGE2D(IX,IY,X0,DX,Y0,DY)
  ! get map data
      call DORESOL2D(.false.,IX,IY,X0,Y0,DX,DY,GRF2D(ig)%NX,GRF2D(ig)%NY,RES2D)
      GRF2D(ig)%NX=MIN(RES2D%NX,MAXPIX)
      GRF2D(ig)%NY=MIN(RES2D%NY,MAXPIX)
      do J=1,GRF2D(ig)%NY
        do I=1,GRF2D(ig)%NX
          GRF2D(ig)%Z(i,j)=RES2D%Z(i,j)
        enddo
      enddo

  ! assign viewport
      call AssignGraph(ig,ip,grtyp_2D)
! VIEWPORT
  ! set viewport screen
      GRFPORT(ip)%DC=(/0.15,0.93,0.31,0.89/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
  ! set physical scale
      GRFPORT(ip)%WC=(/X0-DX/2,X0+DX/2,Y0-DY/2,Y0+DY/2/)
  ! viewport titles etc..
      call getCmdParamS('RES2D',1,'COM',GRFPORT(ip)%comment)
      GRFPORT(ip)%legend=' '
  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))
      END subroutine VIEWRES2D

!----------------------------------------------------------------------
      SUBROUTINE VIEWSCAN1D
! Prepare plot of 1D scan
! Plot parameters are taken from the SCAN1D object
!----------------------------------------------------------------------
      real :: dw
      integer :: ip,ig,i,IVAR,IVAL,is,il
      character(32) :: ESTR
      ip=1  ! port index
  ! set viewport screen
      GRFPORT(ip)%DC=(/0.1,0.9,0.3,0.8/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
! GRAPH DATA
      ig=getNextGraph(grtyp_1D)  ! graph index
  ! number of points
      GRF1D(ig)%NP=SCAN1D%NS
  ! variable index
      IVAR=min(SCAN1D%NVAR,NINT(getCmdParam('SCAN1D','IX')))
      IVAR=max(1,IVAR)
  ! value index
      IVAL=min(SCAN1D%EVAL%N,NINT(getCmdParam('SCAN1D','IV')))
      IVAL=max(1,IVAL)
  ! errors
      GRF1D(ig)%ER=1
  ! get scan data
      GRF1D(ig)%NP=MIN(SCAN1D%NS,MAXPNTS)
      do I=1,GRF1D(ig)%NP
        GRF1D(ig)%X(i)=SCAN1D%VARIABLES(IVAR,I)
        GRF1D(ig)%Y(i)=SCAN1D%VALUES(IVAL,I)
        GRF1D(ig)%DY(i)=SCAN1D%ERRORS(IVAL,I)
      enddo
! GRAPH styles
      GRF1D(ig)%COLOR=1
      GRF1D(ig)%POINT=4
      GRF1D(ig)%LINE=1
      GRF1D(ig)%PSIZE=MAX(0.3,1.5*dw)
! ASSIGN graph
      call AssignGraph(ig,ip,grtyp_1D)
! GRFPORT settings
  ! set physical x-scale
      GRFPORT(ip)%WC(1:2)=(/GRF1D(ig)%X(1),GRF1D(ig)%X(GRF1D(ig)%NP)/)
  ! set physical y-scale
      GRFPORT(ip)%LOGSCALE(2)=0.0
      call SetYscale(ip)
  ! captions
      call FINDSTRPAR(SCAN1D%HDRVAR,':',IVAR,IS,IL)
      GRFPORT(ip)%XTIT=SCAN1D%HDRVAR(IS:IS+IL-1)
      call FINDSTRPAR(SCAN1D%EVAL%HDR,'|',IVAL,IS,IL)
      GRFPORT(ip)%YTIT=SCAN1D%EVAL%HDR(IS:IS+IL-1)
      call GetEnumString('SCAN1D_TYPE',SCAN1D%EVAL%TYP,ESTR)
      GRFPORT(ip)%HEAD='1D scan, '//trim(ESTR)
      select case(SCAN1D%EVAL%TYP)
      case(rset_bpar,rset_rpar)
        GRFPORT(ip)%LEGEND='variable: '//trim(SCAN1D%EVAL%CAP)
      case default
        GRFPORT(ip)%LEGEND=' '
      end select
      GRFPORT(ip)%COMMENT=' '
  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))
! LEGENDS
  ! variable name
      ig=getNextGraph(grtyp_txt)
      GRFTXT(ig)%COLOR=1
      GRFTXT(ig)%DIST=10.0
      GRFTXT(ig)%SHIFT=-0.05
      GRFTXT(ig)%SIDE='B'
      GRFTXT(ig)%LWIDTH=2
      GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
      GRFTXT(ig)%TEXT=GRFPORT(ip)%LEGEND
      call AssignGraph(ig,ip,grtyp_txt)
      END subroutine VIEWSCAN1D

!----------------------------------------------------------------------
      SUBROUTINE VIEWSCAN2D
! Prepare plot of 2D scan
! Plot parameters are taken from the SCAN2D object
!----------------------------------------------------------------------
      integer :: ip,ig,i,j
      real :: dw
      character(32) :: ESTR
      ip=1  ! port index
! GRAPH DATA
      ig=getNextGraph(grtyp_2D)  ! graph index
  ! number of pixels
      GRF2D(ig)%NX=SCAN2D%NX
      GRF2D(ig)%NY=SCAN2D%NY
  ! logscale
      if (getCmdParam('BEAM2D','ZLOG').gt.0.D0) GRFPORT(ip)%LOGSCALE(3)=10.0
      do J=1,GRF2D(ig)%NY
        do I=1,GRF2D(ig)%NX
          GRF2D(ig)%Z(i,j)=SCAN2D%VALUES(i,j)
        enddo
      enddo
  ! assign to viewport
      call AssignGraph(ig,ip,grtyp_2D)
! PORT data
  ! set viewport screen
      GRFPORT(ip)%DC=(/0.15,0.93,0.31,0.89/)
      dw=GRFPORT(ip)%DC(2)-GRFPORT(ip)%DC(1) ! width of the viewport in screen units
  ! set physical scale
      GRFPORT(ip)%WC=(/SCAN2D%XMIN,SCAN2D%XMAX,SCAN2D%YMIN,SCAN2D%YMAX/)
  ! reference points
      GRFPORT(ip)%REF=0.0
  ! captions
      GRFPORT(ip)%XTIT=trim(SCAN2D%HXVAR)
      GRFPORT(ip)%YTIT=trim(SCAN2D%HYVAR)
      call GetEnumString('SCAN2D_TYPE',SCAN2D%EVAL%TYP,ESTR)
      GRFPORT(ip)%HEAD='2D scan, '//trim(ESTR)
      select case(SCAN2D%EVAL%TYP)
      case(rset_bpar,rset_rpar)
        GRFPORT(ip)%LEGEND='variable: '//trim(SCAN2D%EVAL%CAP)
      case default
        GRFPORT(ip)%LEGEND=' '
      end select
      GRFPORT(ip)%COMMENT=' '
  ! optimum character size
      GRFPORT(ip)%CHSIZE=MIN(MAX(0.7,2.3*dw),1.2)
  ! optimum line width
      GRFPORT(ip)%LWIDTH=MIN(1,NINT(3.0*dw))
! LEGENDS
  ! variable name
      ig=getNextGraph(grtyp_txt)
      GRFTXT(ig)%COLOR=1
      GRFTXT(ig)%DIST=10.0
      GRFTXT(ig)%SHIFT=-0.05
      GRFTXT(ig)%SIDE='B'
      GRFTXT(ig)%LWIDTH=2
      GRFTXT(ig)%CHSIZE=0.8*GRFPORT(ip)%CHSIZE
      GRFTXT(ig)%TEXT=GRFPORT(ip)%LEGEND
      call AssignGraph(ig,ip,grtyp_txt)

      END subroutine VIEWSCAN2D


C--------------------------------------------------------
      SUBROUTINE XML_GRFLIST
C Print a list of graphics devices
C format <GRFLIST selected="integer">...</GRFLIST>
C format of items:
C <GRFDEV name="string">/devicename</GRFDEV>
C <GRFFILE>filename</GRFFILE>
C--------------------------------------------------------
      INTEGER :: I,NS
      CHARACTER*10 :: CNUM
      CHARACTER*8 :: SELTYPE,SINT

10    format('<GRFLIST selected="',a,'">')
11    format('<GRFDEV name="',a,'" interactive="',a,'">',a,'</GRFDEV>')
20    FORMAT('    name      description ',/,79('-'))
21    FORMAT(a10,2x,2x,a)
22    FORMAT('selected graphics: ',a,2x,'file: ',a)

  !    write(*,*) 'XML_GRFLIST ', GRFSEL,GRFNUM
      NS=GRFSEL
      if (GRFSEL.GT.0) then
        NS=GRFSEL
        SELTYPE=GRFDEV(NS)%DTYPE
      else
        NS=-1
        SELTYPE=' '
      endif
C// header
      if (XMLOUT.GT.0) then
        call INT2STR(NS,CNUM)
        call XML_RSXDUMP(SMES,' ',1)
        write(SMES,10) trim(CNUM)
      else
        write(SMES,20)
      endif
C// list
      DO I=1,GRFNUM
        if (GRFDEV(I)%INTER.gt.0) then
          SINT='yes'
        else
          SINT='no'
        endif
        IF (XMLOUT.GT.0) then
          WRITE(SMES,11) trim(GRFDEV(I)%DESCR),trim(SINT),trim(GRFDEV(I)%DTYPE)
        else
          WRITE(SMES,21) trim(GRFDEV(I)%DTYPE),trim(GRFDEV(I)%DESCR)
        endif
      ENDDO
c// filename
      if (XMLOUT.GT.0) then
        call XML_TAG(SMES,'GRFFILE',trim(GRFFILE),2)
      else
        write(SMES,22) TRIM(SELTYPE),trim(GRFFILE)
      endif
C// end tags
      call XML_TAG(SMES,'GRFLIST',' ',0)
      call XML_RSXDUMP(SMES,' ',0)
      return

      END SUBROUTINE XML_GRFLIST

!-------------------------------------------------
      SUBROUTINE SELGRFDEV(SARG)
!  Select an existing graphics device
!  if SARG given, use it as requested device string,
!  otherwise start dialog
!  no info if IQUIET>0
!-------------------------------------------------
      CHARACTER*(*) SARG
      INTEGER*4 ISEL,I,ISS,LSS,IPOS,ILA,LLA
      CHARACTER*128 S
      integer IDCOMPARE

1     FORMAT(a)
3     FORMAT('selected graphics device: ',a,2x,'(',I3,')')
      ISEL=0 ! selection index

   !   write(*,*) 'SELGRFDEV ',trim(SARG)
!// with SARG empty, print the list and exit
      IF (SARG.EQ.' ') THEN
        call XML_GRFLIST
      else
        CALL BOUNDS(SARG,ISS,LSS)
        S=SARG(ISS:ISS+LSS-1)
        ISS=1
    ! find device type substring
        CALL LASTSUBSTR(S,'/',IPOS)
        ISEL=0
        IF (IPOS.GT.0) THEN
          ILA=MIN(LSS-IPOS+1,8) ! max. 8 characters for device ID string
          LLA=IPOS+ILA-1 ! last character of the device ID string
! search for matching device
          ISEL=0
          I=1
          DO WHILE ((ISEL.LE.0).AND.(I.LE.GRFNUM))
            if (IDCOMPARE(GRFDEV(I)%DTYPE,S(IPOS:LLA)).eq.0) ISEL=I
            I=I+1
          ENDDO
        ENDIF
    ! new device found
        GRFFILE=' '
        IF (ISEL.GT.0) THEN
      ! get filename if given
          GRFFILE=' '
          IF (IPOS.GT.ISS) GRFFILE=S(ISS:IPOS-1)
          GRFSEL=ISEL
        ELSE ! not found
          call MSG_WARN("device not available: "//TRIM(S),1)
          GRFSEL=0
        endif
        call getDeviceStr(GRFSEL,S)
        WRITE(*,*) "device is: ",TRIM(S)
! report the list with currently selected device
        if (XMLOUT.gt.0) call XML_GRFLIST
      ENDIF
      END SUBROUTINE SELGRFDEV

!----------------------------------------------------------------------
      SUBROUTINE SavePortData(iport,ig,fname)
! Saves graph data for given viewport
! iport ... viewport index
! ig ... graph part for file suffix
! fname  ... filename
!----------------------------------------------------------------------
      INTEGER, intent(in) :: iport,ig
      character*(*) :: fname
      integer :: IO
      integer :: ires
      character*256 :: S,MSG
      character*32 :: CNUM
1     format('CFG: "',a,'"')
!c// get filename interactively if needed
      ires=0
      if (len_trim(fname).LE.0) then
        call DLG_FILEOPEN(fname,' ','dat',0,1,ires,S)
        if (len_trim(S).gt.0) ires=1
      else
        S=trim(fname)
        ires=1
      endif
      if (ires.le.0) goto 200
!c// don't save empty viewports
      if(GRFPORT(iport)%NGRF.le.0) ires=0
      IF (ires.GT.0) THEN
        CNUM=' '
        if (ig>0) call INT2STR(ig,CNUM)
        !IO=OPENFILEUNIT(trim(S),.false.)
        IO=OPENFILEOUT_SFX(trim(S),CNUM)
        if (IO.gt.0) then
          write(IO,1) trim(XCONFIG_NAME)
          call WritePortData(IO,iport)
          if (IO.ne.6) close(IO)
        else
          goto 999
        endif
      ENDIF
      return

200   write(MSG,201) iport
201   format('No filename for port ',I2,', graph data not saved')
      call MSG_WARN(trim(MSG),1)
      return
999   write(MSG,*) 'Cannot open output file ',trim(S),'as unit ',IO
      call MSG_ERROR('SavePortData',trim(MSG),1,1)
      END SUBROUTINE SavePortData

!----------------------------------------------------------------------
      SUBROUTINE SaveAllPorts(fname)
! Saves all graph data from all viewports
! fname  ... base filename, index is appended for ports 2.
!----------------------------------------------------------------------
      character*(*) :: fname
      integer :: i,nip,iport(MAXPORTS)
!c/ no action without filename
      if (len_trim(fname).le.0) return
!c// collect indices of  ports with some data assigned
      nip=0
      do i=1,MAXPORTS
        if(GRFPORT(i)%NGRF.gt.0) then
          nip=nip+1
          iport(nip)=i
        endif
      enddo
      if (nip.le.0) goto 100
!C// scan all assigned ports and save each in a file
      do i=1,nip
    ! get output filename, append index for multiple ports
        if (nip.gt.1) then
          call SavePortData(iport(i),i,fname)
        else
          call SavePortData(iport(i),0,fname)
        endif
      enddo
      return
100   call MSG_WARN('No data to save ... ',1)
      END SUBROUTINE SaveAllPorts



!----------------------------------------------------------------------
      SUBROUTINE SaveAuxData(fname)
! Save auxilliary data from RES1D_BINS if available
! fname  ... filename
!----------------------------------------------------------------------
      character*(*) :: fname
      integer :: IO
      integer :: ires
      character*256 :: S,MSG
      character*32 :: CNUM
1     format('CFG: "',a,'"')
!c// get filename interactively if needed
      ires=0
      !c// don't save empty viewports
      if (RES1D_BINMAX<RES1D_BINMIN) return
      if (len_trim(fname).LE.0) then
        call DLG_FILEOPEN(fname,' ','dat',0,1,ires,S)
        if (len_trim(S).gt.0) ires=1
      else
        S=trim(fname)
        ires=1
      endif
      if (ires.le.0) goto 200
      IF (ires.GT.0) THEN
        CNUM='bins'
        !IO=OPENFILEUNIT(trim(S),.false.)
        IO=OPENFILEOUT_SFX(trim(S),CNUM)
        if (IO.gt.0) then
          write(IO,1) trim(XCONFIG_NAME)
          call WriteBinData(IO)
          if (IO.ne.6) close(IO)
        else
          goto 999
        endif
      ENDIF
      return

200   write(MSG,201)
201   format('No output filename for graph bin data, data not saved')
      call MSG_WARN(trim(MSG),1)
      return
999   write(MSG,*) 'Cannot open output file ',trim(S),'as unit ',IO
      call MSG_ERROR('SavePortData',trim(MSG),1,1)
      END SUBROUTINE SaveAuxData




      end MODULE GRFEXEC

