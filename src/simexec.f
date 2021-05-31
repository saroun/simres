!//////////////////////////////////////////////////////////////////////
!////  $Id: simexec.f,v 1.64 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2011, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.64 $
!////     $Date $
!//////////////////////////////////////////////////////////////////////
!////
!////  Top level execution commands for SIMRES
!////
!//////////////////////////////////////////////////////////////////////
      MODULE SIMEXEC
  !    use FIELDS
      USE FIELDDATA
      use RESULTS
      use REPORTS
      use DIALOGS
      use TRACING
      use COMPONENTS
      USE RESOL1D
      use BEAM1D
      use GRFEXEC

      implicit none

      contains

!----------------------------------------------------------------------
      SUBROUTINE DO_FMERIT(ID,valFM,errFM)
! adjust crystal with given ID (if exists)
!----------------------------------------------------------------------
      integer,intent(in) :: ID
      real(kind(1.D0)),intent(out) :: valFM(3),errFM(3)
      real(kind(1.D0)) :: norm,FM,EFM,VALUES(3),ERRORS(3),PAR(3),dPAR(3)
      real(kind(1.D0)) :: Z, XLIM(2)
      integer :: typ,fml,wexp,filt
      TYPE(TEVAL) :: EVAL
      character(32) :: CNUM1,CNUM2,CNUM3
      logical :: mpeak
      !REAL :: MPAR(3*MAXPEAKS),dMPAR(3*MAXPEAKS)

      typ=NINT(getCmdParam('FMERIT','TYPE'))
      EVAL%TYP=typ
      call DEFEVAL(EVAL)
      mpeak=((EVAL%TYP==rset_dpar).and.HAS_MULTIPEAK_DATA(EVAL%IV))
      if (mpeak) then
        call DO_FMERIT_MPK(ID,valFM,errFM)
        return
      endif
      VALUES=0.D0
      ERRORS=0.D0
      valFM=0.D0
      errFM=0.D0
      FM=0.D0
      EFM=0.D0
! run simulation
      call TRACING_PROC
! start of FM calculation
      if (TRACING_INT>0.D0) then
      fml=NINT(getCmdParam('FMERIT','FORMULA'))
      norm=getCmdParam('FMERIT','NORM')
      !write(*,*) 'FMERIT NORM=',NORM
      ! cost=getCmdParam('FMERIT','COST')
      wexp=NINT(getCmdParam('FMERIT','WEXP'))
      filt=NINT(getCmdParam('BEAM1D','FILT'))
      EVAL%TYP=typ
      call DEFEVAL(EVAL)
      call GET_EVAL_RANGE(EVAL, XLIM)
      select case(EVAL%TYP)
        case(rset_basic)
          VALUES(1:3)=(/TRACING_INT,TRACING_DE,TRACING_TOF/)
          ERRORS(1:3)=(/TRACING_dINT,TRACING_dDE,TRACING_dTOF/)
        case(rset_bpar)
          call GET_EVAL_RANGE(EVAL,XLIM)
          call BEAMPAR(EVAL%IV,show_monitor,filt,XLIM,PAR,dPAR)
          VALUES(1:3)=PAR(1:3)
          ERRORS(1:3)=dPAR(1:3)
        case(rset_dpar)
          call GET_EVAL_RANGE(EVAL,XLIM)
          call BEAMPAR(EVAL%IV,show_detector,filt,XLIM,PAR,dPAR)
          Z=1.0
          if ((wexp>0).and.(abs(PAR(3))>1.e-3)) Z=PAR(3)**wexp
          VALUES(1)=PAR(1)
          ERRORS(1)=dPAR(1)
          VALUES(2)=PAR(2)/Z
          ERRORS(2)=dPAR(2)/Z
          VALUES(3)=PAR(3)
          ERRORS(3)=dPAR(3)
        case(rset_rpar)
          call RESOLPAR(EVAL%IV,PAR,dPAR)
          VALUES(1:3)=PAR(1:3)
          ERRORS(1:3)=dPAR(1:3)
      end select
      if ((VALUES(1)>0.D0).and.(VALUES(2)>0.D0)) then
        select case (fml)
        case(fmerit_int)
          FM=VALUES(1)
          EFM=ERRORS(1)
        case(fmerit_w)
          FM=1.D0/VALUES(2)
          EFM=ERRORS(2)/VALUES(2)**2
        case(fmerit_intw)
          FM=VALUES(1)/VALUES(2)
          EFM=SQRT((ERRORS(1)/VALUES(1))**2 + (ERRORS(2)/VALUES(2))**2)*FM
        case(fmerit_intw2)
          FM=VALUES(1)/VALUES(2)**2
          EFM=SQRT((ERRORS(1)/VALUES(1))**2 + (2.D0*ERRORS(2)/VALUES(2))**2)*FM
        end select
      endif
      !write(*,*) 'FMERIT result=',FM, FM/NORM
      FM=FM/NORM
      EFM=EFM/NORM
      valFM(1)=FM
      errFM(1)=eFM
      valFM(2:3)=VALUES(1:2)
      errFM(2:3)=ERRORS(1:2)
      endif
! end of FM calculation
      if (id==-2) return
      call INT2STR(typ,CNUM1)
      call INT2STR(fml,CNUM2)
      call INT2STR(id,CNUM3)
      call XML_RSXDUMP(SOUT,' ',1)
1     format('<',a,' type="',a,'" formula="',a,'" id="',a,'">')
      write(SOUT,1) 'FMERIT',trim(CNUM1),trim(CNUM2),trim(CNUM3)
      call FLOAT2STR(FM,CNUM1)
      call FLOAT2STR(EFM,CNUM2)
      call XML_TAG(SOUT,'FM',trim(CNUM1),2)
      call XML_TAG(SOUT,'EFM',trim(CNUM2),2)
      call FLOAT2STR(VALUES(1),CNUM1)
      call FLOAT2STR(ERRORS(1),CNUM2)
      call XML_TAG(SOUT,'INT',trim(CNUM1),2)
      call XML_TAG(SOUT,'EINT',trim(CNUM2),2)
      call FLOAT2STR(VALUES(2),CNUM1)
      call FLOAT2STR(ERRORS(2),CNUM2)
      call XML_TAG(SOUT,'WIDTH',trim(CNUM1),2)
      call XML_TAG(SOUT,'EWIDTH',trim(CNUM2),2)
      call XML_RSXDUMP(SOUT,'FMERIT',0)

      end SUBROUTINE DO_FMERIT


!----------------------------------------------------------------------
      SUBROUTINE DO_FMERIT_MPK(ID,valFM,errFM)
! adjust crystal with given ID (if exists)
!----------------------------------------------------------------------
      integer,intent(in) :: ID
      real(kind(1.D0)),intent(out) :: valFM(3),errFM(3)
      real(kind(1.D0)) :: norm,FM,EFM,aFM,aEFM,VALUES(3),ERRORS(3),XSTEP,CHISQR,Z
      real(kind(1.D0)) :: SVAL(2),SERR(2),chisqlim
      integer :: typ,fml,NP,i,k,wexp
      TYPE(TEVAL) :: EVAL
      character(32) :: CNUM1,CNUM2,CNUM3
      REAL :: MPAR(3*MAXPEAKS),dMPAR(3*MAXPEAKS)
      VALUES=0.D0
      ERRORS=0.D0
      valFM=0.D0
      errFM=0.D0
      SVAL=0.D0
      SERR=0.D0
      aFM=0.D0
! run simulation
      call TRACING_PROC
! start of FM calculation
      if (TRACING_INT>0.D0) then
      FM=0.D0
      eFM=0.D0
      typ=NINT(getCmdParam('FMERIT','TYPE'))
      fml=NINT(getCmdParam('FMERIT','FORMULA'))
      norm=getCmdParam('FMERIT','NORM')
      ! cost=getCmdParam('FMERIT','COST')
      chisqlim=getCmdParam('FMERIT','CHI2LIM')
      wexp=NINT(getCmdParam('FMERIT','WEXP'))
      EVAL%TYP=typ
      call DEFEVAL(EVAL)
! fit multiple peaks
      CHISQR=BEAMPAR_MPK(EVAL%IV,MPAR,dMPAR,NP)
      if (CHISQR>=chisqlim) then
        call FLOAT2STR(CHISQR,CNUM1)
        call MSG_WARN('Poor fitting precision, CHI2='//trim(CNUM1),1)
        write(*,*) 'Poor fitting precision, CHI2='//trim(CNUM1)
        write(*,*) 'Try better statistics or smaller bin number.'
      endif
      ! sum of integral intensities for all peaks
      XSTEP=(RES1D%X(RES1D%N)-RES1D%X(1))/(RES1D%N-1)
! too high CHISQR signals an invalid result
      if (CHISQR<chisqlim) then
        do i=1,NP
        k=3*(i-1)
        if ((dMPAR(k+1)<MPAR(k+1)).and.(dMPAR(k+3)<MPAR(k+3))) then
      ! ignore peaks centered outside data range
        if ((RES1D%X(1)<MPAR(k+2)).and.(RES1D%X(RES1D%N)>MPAR(k+2))) then
            VALUES(1)=MPAR(k+1)*MPAR(k+3)*sqrt(TWOPI)/R8LN2/XSTEP
            ERRORS(1)=sqrt((dMPAR(k+1)*MPAR(k+3))**2 + (MPAR(k+1)*dMPAR(k+3))**2)*sqrt(TWOPI)/R8LN2/XSTEP
            Z=1.0
            if ((wexp>0).and.(abs(MPAR(k+3))>1.e-3)) Z=MPAR(k+3)**wexp
            VALUES(2)=MPAR(k+2)/Z
            ERRORS(2)=dMPAR(k+2)/Z
            select case (fml)
            case(fmerit_int)
                  aFM=VALUES(1)
                  aEFM=ERRORS(1)
            case(fmerit_w)
                  aFM=1.D0/VALUES(2)
                  aEFM=ERRORS(2)/VALUES(2)**2
            case(fmerit_intw)
                  aFM=VALUES(1)/VALUES(2)
                  aEFM=SQRT((ERRORS(1)/VALUES(1))**2 + (ERRORS(2)/VALUES(2))**2)*aFM
            case(fmerit_intw2)
                  aFM=VALUES(1)/VALUES(2)**2
                  aEFM=SQRT((ERRORS(1)/VALUES(1))**2 + (2.D0*ERRORS(2)/VALUES(2))**2)*aFM
            case default
                  aFM=VALUES(1)
                  aEFM=ERRORS(1)
            end select
2         format(a,2x,6(G13.5,2x))
        !write(*,2) 'DO_FMERIT_MPK ',i,MPAR(k+2),aFM,aEFM
            FM=FM+aFM
            eFM=eFM+aEFM**2
            SVAL(1)=SVAL(1)+VALUES(1)
            SVAL(2)=SVAL(2)+VALUES(2)
            SERR(1)=SERR(1)+ERRORS(1)**2
            SERR(2)=SERR(2)+ERRORS(2)**2
        endif
        endif
        enddo
        eFM=sqrt(eFM)
      !write(*,2) 'DO_FMERIT_MPK total',FM,eFM
        FM=FM/NORM
        EFM=EFM/NORM
        if (abs(FM)>1.D10) FM=0D0
        if (abs(eFM)>1.D10) eFM=0D0
        SERR(1)=sqrt(SERR(1))
        SERR(2)=sqrt(SERR(2))
        SVAL(2)=CHISQR
        SERR(2)=0
      endif
      valFM(1)=FM
      errFM(1)=eFM
      valFM(2:3)=SVAL(1:2)
      errFM(2:3)=SERR(1:2)
      endif
! end of FM calculation
      if (id==-2) return
      call INT2STR(typ,CNUM1)
      call INT2STR(fml,CNUM2)
      call INT2STR(id,CNUM3)
      call XML_RSXDUMP(SOUT,' ',1)
1     format('<',a,' type="',a,'" formula="',a,'" id="',a,'">')
      write(SOUT,1) 'FMERIT',trim(CNUM1),trim(CNUM2),trim(CNUM3)
      call FLOAT2STR(FM,CNUM1)
      call FLOAT2STR(EFM,CNUM2)
      call XML_TAG(SOUT,'FM',trim(CNUM1),2)
      call XML_TAG(SOUT,'EFM',trim(CNUM2),2)
      call FLOAT2STR(SVAL(1),CNUM1)
      call FLOAT2STR(SERR(1),CNUM2)
      call XML_TAG(SOUT,'INT',trim(CNUM1),2)
      call XML_TAG(SOUT,'EINT',trim(CNUM2),2)
      call FLOAT2STR(SVAL(2),CNUM1)
      call FLOAT2STR(SERR(2),CNUM2)
      call XML_TAG(SOUT,'WIDTH',trim(CNUM1),2)
      call XML_TAG(SOUT,'EWIDTH',trim(CNUM2),2)
      call XML_RSXDUMP(SOUT,'FMERIT',0)

      end SUBROUTINE DO_FMERIT_MPK

!----------------------------------------------------------------------
      SUBROUTINE DO_ADJCRYST(ID)
! adjust crystal with given ID (if exists)
!----------------------------------------------------------------------
      character(*) :: ID
      integer :: ierr,N,IS,IL
      REAL :: AUX(2)
      type(TCOMPID) :: COMPID
      TYPE(TCLASS) :: CLS
      character(LEN_NAME) :: CID
      character(32) :: CPAR
      if (.not.associated(INSOBJ%P_SPEC)) return
  ! split additiona arguments from the component ID
      call FINDSTRPAR(trim(ID),' ',1,IS,IL)
      CID=ID(IS:IS+IL-1)
      CPAR=' '
      if (IL.lt.len_trim(ID)) CPAR=trim(ID(IS+IL:))
  ! identify component by ID
      CLS=getClassByID(trim(CID))
      select case(CLS%TCLS)
      case(TCLS_COM)
        ierr=1
        if (CLS%P_COM%CCLS.eq.CCLS_CRYSTAL) then
          call GetComponentByID(CID,COMPID)
          select case(COMPID%IBEAM)
          case(2)
            call GENERATE_NOM(INSOBJ%P_SPEC%KF,NEUT)
          case default
            call GENERATE_NOM(INSOBJ%P_SPEC%KI,NEUT)
          end select
          call CRYSTAL_ADJUST(CLS%P_COM%INST,IERR)
        endif
      case(TCLS_SAM)
  ! additional parameter means scattering angle
        AUX=0.0
        call STR2ARRAY4(CPAR,' ',AUX,2,N)
        if (N.ge.1) then
          call INST_SET_THETAS(AUX(1)*PI/180.D0,IERR)
        else
          call INST_ADJUST_SAMPLE(IERR)
        endif
      end select
      if (ierr.ne.0) call MSG_WARN('Can''t adjust '//trim(ID),1)
      end SUBROUTINE DO_ADJCRYST

!----------------------------------------------------------------------
      SUBROUTINE REPCRYST(ID)
! report on crystal reflecting properties
!----------------------------------------------------------------------
      character(*) :: ID
      TYPE(TCOMOBJ) :: CLS
      real(KIND(1.D0)) :: AUX(4)
      TYPE(PCRYSTAL) :: PCR
      ! CLS=getClassByID(trim(ID))
      CLS=getComByID(ID)
      if (CLS%CCLS.eq.CCLS_CRYSTAL) then
          call CRYSTAL_INIT(CLS%INST)
          call CRYSTAL_GET(CLS%INST,PCR)
          AUX(1)=PCR%X%QHKL
          AUX(2)=PCR%X%REF
          AUX(3)=PCR%X%MI
          AUX(4)=PCR%X%DW
          call XML_ARRAY(SOUT,4,trim(CLS%ID)//' properties ','Qhkl|ref|mu|DW',AUX(1))
      endIF
      end SUBROUTINE REPCRYST

!----------------------------------------------------------------------
      SUBROUTINE REPXTAL(ID)
! report on XTAL reflecting properties
!----------------------------------------------------------------------
      character(*) :: ID
      TYPE(TCOMOBJ) :: CLS
      TYPE(PXTAL) :: PCR
      CLS=getComByID(ID)
      if (CLS%CCLS.eq.CCLS_XTAL) then
        call XTAL_REPORT_ABS(CLS%INST)
        call XTAL_GET(CLS%INST,PCR)
        call XREF_REPORT(PCR%X%ICR)
      endIF
      end SUBROUTINE REPXTAL

!----------------------------------------------------------------------
      SUBROUTINE DEFEVAL(EVAL)
! defines the quantity during evaluated scans
!----------------------------------------------------------------------
      TYPE(TEVAL),intent(out) :: EVAL
      integer :: auto
      EVAL%V0=0.0
      EVAL%DV=0.0
 ! select result type and headers
      select case(EVAL%TYP)
      case(rset_basic)
        EVAL%HDR='INT|DE|TOF'
        EVAL%N=3
        EVAL%CAP=' '
        EVAL%IV=0
      case(rset_bpar,rset_bpeak)
        ! variable index depends on the setting of BEAM1D command !
        EVAL%IV=MAX(1,NINT(getCmdParam('BEAM1D','X'))+1)
        EVAL%IV=MIN(EVAL%IV,RK_MAX)
        call VAR_GETNAME(EVAL%IV,.true.,EVAL%CAP)
        if (EVAL%TYP.eq.rset_bpeak) then
          EVAL%V0=getCmdParam('BEAM1D','X0')
          EVAL%DV=getCmdParam('BEAM1D','DX')
          EVAL%N=NINT(getCmdParam('BEAM1D','NP'))
          EVAL%HDR='INT'
        else
          EVAL%HDR='INT|SPREAD|CENTER'
          EVAL%N=3
        endif
      case(rset_dpar,rset_dpeak)
        ! variable index depends on the setting of DET1D command !
        EVAL%IV=MAX(1,NINT(getCmdParam('DET1D','X'))+1)
        EVAL%IV=MIN(EVAL%IV,DT_MAX)
        call DVAR_GETNAME(EVAL%IV,.true.,EVAL%CAP)
        auto=NINT(getCmdParam('DET1D','XAUTO'))
        if (auto==0) then
          EVAL%V0=getCmdParam('DET1D','X0')
          EVAL%DV=getCmdParam('DET1D','DX')
        endif
        if (EVAL%TYP.eq.rset_dpeak) then
          EVAL%N=NINT(getCmdParam('DET1D','NP'))
          EVAL%HDR='INT'
        else
          EVAL%N=3
          EVAL%HDR='INT|SPREAD|CENTER'
        endif
      case(rset_rpar,rset_rpeak)
        ! variable index depends on the setting of RES1D comand !
        EVAL%IV=MAX(1,NINT(getCmdParam('RES1D','X'))+1)
        EVAL%IV=MIN(EVAL%IV,QE_MAX)
        call QVAR_GETNAME(EVAL%IV,.true.,EVAL%CAP)
        if (EVAL%TYP.eq.rset_rpeak) then
          EVAL%V0=getCmdParam('RES1D','X0')
          EVAL%DV=getCmdParam('RES1D','DX')
          EVAL%N=NINT(getCmdParam('RES1D','NP'))
          EVAL%HDR='INT'
        else
          EVAL%HDR='INT|SPREAD|CENTER'
          EVAL%N=3
        endif
      case(rset_fmerit)
        EVAL%HDR='FM|INT|SPREAD'
        EVAL%N=3
        EVAL%CAP='Figure of merit'
        EVAL%IV=1
        EVAL%V0=0.0
        EVAL%DV=0.0
      end select
      end subroutine DEFEVAL

!----------------------------------------------------------------------
      SUBROUTINE DOSCAN1D
! Execute 1D scan based on the definition in SCAN1D comand
!----------------------------------------------------------------------
      integer :: i,is,il,ires, ierr
      logical :: log1
      character(32) :: CNUM
      character(MAX_FNAME_LENGTH) :: FN
      !write(*,*) 'DOSCAN1D'
      call DEFSCAN1D(ierr)
      if (ierr.ne.0) then
        call MSG_WARN('SCAN1D: cannot execute scan',1)
        return
      endif
      !write(*,*) 'DEFSCAN1D'
  ! check variables
      LOG1=.true.
      if (SCAN1D%NS.le.1) then
        call MSG_WARN('SCAN1D: number of steps must be > 1',1)
        LOG1=.false.
      endif
      do i=1,SCAN1D%NVAR
        ! call FINDSTRPAR(SCAN1D%HDRVAR,'|',i,IS,IL)
        call FINDSTRPAR(SCAN1D%HDRVAR,':',i,IS,IL)
        if (IL.gt.0) then
          ires=DLG_PARAM(SCAN1D%HDRVAR(IS:IS+IL-1),.false.)
          if (ires.le.0) then
            call MSG_WARN('SCAN1D: Unknown variable: '//SCAN1D%HDRVAR(IS:IS+IL-1),1)
            LOG1=.false.
            exit
          endif
        else
          LOG1=.false.
          call INT2STR(i,CNUM)
          call MSG_WARN('SCAN1D: Empty ID for variable '//trim(CNUM),1)
          exit
        endif
      enddo
      if (.not.log1) return
  ! run the scan
      !write(*,*) 'starting RUNSCAN1D'
      call RUNSCAN1D(.true.)
  ! and save it
      call getCmdParamS('SCAN1D',1,'FILE',FN)
      if (len_trim(FN).gt.0) then
        checkOverwrite=.false.
        call SAVESCAN1D(trim(FN))
      endif
      end SUBROUTINE DOSCAN1D


!----------------------------------------------------------------------
      SUBROUTINE DOSCAN2D
! Execute 2D scan based on the definition in SCAN2D comand
!----------------------------------------------------------------------
      integer :: i,ires,IL
      logical :: log1
      character(32) :: CNUM
      character(MAX_FNAME_LENGTH) :: FN
      call DEFSCAN2D
  ! check variables
      LOG1=.true.
      if ((SCAN2D%NX.le.1).or.(SCAN2D%NY.le.1)) then
        call MSG_WARN('SCAN2D: number of steps must be > 1',1)
        LOG1=.false.
      endif
      IL=len_trim(SCAN2D%HXVAR)
      if (IL.gt.0) then
        ires=DLG_PARAM(trim(SCAN2D%HXVAR),.false.)
        if (ires.le.0) then
          call MSG_WARN('SCAN2D: Unknown variable: '//trim(SCAN2D%HXVAR),1)
          LOG1=.false.
        endif
      else
        LOG1=.false.
        call INT2STR(i,CNUM)
        call MSG_WARN('SCAN2D: Empty ID for variable '//trim(CNUM),1)
      endif
      IL=len_trim(SCAN2D%HYVAR)
      if (IL.gt.0) then
        ires=DLG_PARAM(trim(SCAN2D%HYVAR),.false.)
        if (ires.le.0) then
          call MSG_WARN('SCAN2D: Unknown variable: '//trim(SCAN2D%HYVAR),1)
          LOG1=.false.
        endif
      else
        LOG1=.false.
        call INT2STR(i,CNUM)
        call MSG_WARN('SCAN2D: Empty ID for variable '//trim(CNUM),1)
      endif
      if (SCAN2D%EVAL%TYP.eq.rset_bpeak) then
        LOG1=.false.
        call MSG_WARN('SCAN2D: peak profiles can''t be recorded in 2D scan ',1)
      endif
      if (.not.log1) return
  ! run the scan
      call RUNSCAN2D
  ! and save it
      call getCmdParamS('SCAN2D',1,'FILE',FN)
      if (len_trim(FN).gt.0) then
        checkOverwrite=.false.
        call SAVESCAN2D(trim(FN))
      endif
      end SUBROUTINE DOSCAN2D

!----------------------------------------------------------------------
      SUBROUTINE DEFSCAN1D(ierr)
! Define 1D scan with current command settings
!----------------------------------------------------------------------
      integer, intent(out) :: ierr
      integer :: i,N1,N2,N0,ip
      character(FIELD_BUFFER_STR) :: CPAR
      character(32) :: CNUM
      real :: STA(MAXVARS),ENDS(MAXVARS),STP
      ierr=0
  ! variable ID strings, | delimitted
      call getCmdParamS('SCAN1D',1,'VARS',SCAN1D%HDRVAR)
  ! start values
      call getCmdParamS('SCAN1D',1,'START',CPAR)
      call COUNTPAR(CPAR,':',N0)
      call STR2ARRAY4(CPAR,':',STA,MAXVARS,N1)
      write(*,*) 'DEFSCAN1D: ',trim(SCAN1D%HDRVAR),trim(CPAR),N0,N1

      if (N1.ne.N0) then
        write(*,*) 'DEFSCAN1D: N0=',N0,' N1=',N1
        write(*,*) 'DEFSCAN1D: PAR=',trim(CPAR)
        call MSG_ERROR('SCAN1D','Syntax error on the start values list.',0,1)
        return
      endif
  ! end values
      call getCmdParamS('SCAN1D',1,'END',CPAR)
      call COUNTPAR(CPAR,':',N0)
      call STR2ARRAY4(CPAR,':',ENDS,MAXVARS,N2)
      if (N2.ne.N0) then
        call MSG_ERROR('SCAN1D','Syntax error on the end values list.',0,1)
        return
      endif
  ! check numer of variables
      if (N1.ne.N2) then
        call MSG_WARN('Unequal number of start and end values.',1)
        return
      endif
      SCAN1D%NVAR=N1
  ! number of steps
      SCAN1D%NS=NINT(getCmdParam('SCAN1D','NS'))
      SCAN1D%NS=MIN(SCAN1D%NS,MAXSCAN1D)
  ! value index (e.g. 1 for INT if TYP=BASIC
      SCAN1D%IV=NINT(getCmdParam('SCAN1D','IV'))
  ! coordinate index from the list of variables VARS)
      SCAN1D%IX=NINT(getCmdParam('SCAN1D','IX'))
      SCAN1D%IX=MIN(SCAN1D%IX,MAXVARS)
  ! result type
      SCAN1D%EVAL%TYP=NINT(getCmdParam('SCAN1D','TYPE'))
      call DEFEVAL(SCAN1D%EVAL)
  ! check IV value
      if ((SCAN1D%IV.le.0).or.(SCAN1D%IV>SCAN1D%EVAL%N)) then
        call INT2STR(SCAN1D%IV,CNUM)
        call MSG_WARN('Wrong value index, IV='//trim(CNUM)//'. IV is set to 1',1)
        SCAN1D%IV=1
      endif
  ! re-allocate variables
      call ALLOC_SCAN1D(SCAN1D%NVAR, SCAN1D%NS, SCAN1D%EVAL%N, ierr)
      if (ierr==0) then
  ! get X-values
        do i=1,SCAN1D%NVAR
          if (SCAN1D%NS.gt.1) then
            STP=(ENDS(i)-STA(i))/(SCAN1D%NS-1)
          else
            STP=0.
          endif
          do ip=1,SCAN1D%NS
            SCAN1D%VARIABLES(i,ip)=STA(i)+(ip-1)*STP
          enddo
        enddo
      endif
      end SUBROUTINE DEFSCAN1D

!----------------------------------------------------------------------
      SUBROUTINE DEFSCAN2D
! Define 2D scan with current command settings
!----------------------------------------------------------------------
      character(32) :: CNUM
      integer :: ierr
  ! variable ID strings
      call getCmdParamS('SCAN2D',1,'XVAR',SCAN2D%HXVAR)
      call getCmdParamS('SCAN2D',1,'YVAR',SCAN2D%HYVAR)
  ! variables range
      SCAN2D%XMIN=getCmdParam('SCAN2D','XMIN')
      SCAN2D%XMAX=getCmdParam('SCAN2D','XMAX')
      SCAN2D%YMIN=getCmdParam('SCAN2D','YMIN')
      SCAN2D%YMAX=getCmdParam('SCAN2D','YMAX')
      SCAN2D%NX=NINT(getCmdParam('SCAN2D','NX'))
      SCAN2D%NY=NINT(getCmdParam('SCAN2D','NY'))
      SCAN2D%IV=NINT(getCmdParam('SCAN2D','IV'))
      SCAN2D%NX=MIN(SCAN2D%NX,MAXSCAN2D)
      SCAN2D%NY=MIN(SCAN2D%NY,MAXSCAN2D)
  ! result type
      SCAN2D%EVAL%TYP=NINT(getCmdParam('SCAN2D','TYPE'))
      call DEFEVAL(SCAN2D%EVAL)
  ! check IV value
      if ((SCAN2D%IV.le.0).or.(SCAN2D%IV>SCAN2D%EVAL%N)) then
        call INT2STR(SCAN2D%IV,CNUM)
        call MSG_WARN('Wrong value index, IV='//trim(CNUM)//'. IV is set to 1',1)
        SCAN2D%IV=1
      endif
      select case (SCAN2D%EVAL%TYP)
      case(rset_basic,rset_bpar,rset_rpar,rset_dpar)
        if (SCAN2D%IV.GT.3) then
          call MSG_WARN('value index (IV) is too high, set to 1.',1)
          SCAN2D%IV=1
        endif
      end select
      call ALLOC_SCAN1D(1,SCAN2D%NX,SCAN2D%EVAL%N, ierr)
      end SUBROUTINE DEFSCAN2D

!----------------------------------------------------------------------
      SUBROUTINE RUNSCAN1D(showProgress)
! Run 1D scan with current settings
!----------------------------------------------------------------------
      logical,intent(in) :: showProgress
      type(TREPOPT) :: OPT
      type(TRES1D) :: RES
      real(KIND(1.D0)) :: PAR(3),dPAR(3), XLIM(2)
      integer :: iscan,ires,i,is,il,IS1,IL1,ierr,filt
      character*32 :: CNUM,CNUM0,CNUM1,CNUM2
1     format(a, ': ',6(G12.5,1x))
    ! switch reporting off
      OPT=REPOPT
      call REPOPT_MUTE(REPOPT)
      if (showProgress) call XML_PROGRESS(SMES,0,SCAN1D%NS,0,'1D scan')
  !    write(*,*) 'RUNSCAN1D' ,SCAN1D%IX,SCAN1D%IV,SCAN1D%EVAL%IV
      do iscan=1,SCAN1D%NS
      ! set parameers for iscan-th step
            ires=0
            do i=1,SCAN1D%NVAR
                  call FINDSTRPAR(SCAN1D%HDRVAR,':',i,IS,IL)
                  call FLOAT2STR(1.D0*SCAN1D%VARIABLES(i,iscan),CNUM)
                  ires=DLG_PARAM(SCAN1D%HDRVAR(IS:IS+IL-1)//' '//trim(CNUM),.true.)
                  if (ires.le.0) exit
            enddo
            if (ires.le.0) exit
      ! print header
            if (iscan.eq.1) then
                  call FINDSTRPAR(SCAN1D%EVAL%HDR,'|',SCAN1D%IV,IS,IL)
                  call FINDSTRPAR(SCAN1D%HDRVAR,':',SCAN1D%IX,IS1,IL1)
                  write(*,"(a10,a20,a20)") "step",SCAN1D%HDRVAR(IS1:IS1+IL1-1),SCAN1D%EVAL%HDR(IS:IS+IL-1)
            endif
      ! run MC
        call INST_ADJUST(IERR,'RUNSCAN1D')
        if (IERR.eq.0) then
          if (SCAN1D%EVAL%TYP==rset_fmerit) then
            call DO_FMERIT(-2,PAR,dPAR)
          else
            call TRACING_PROC
          endif
        endif
      ! store results
        select case(SCAN1D%EVAL%TYP)
        case(rset_basic)
          SCAN1D%VALUES(1:3,iscan)=(/TRACING_INT,TRACING_DE,TRACING_TOF/)
          SCAN1D%ERRORS(1:3,iscan)=(/TRACING_dINT,TRACING_dDE,TRACING_dTOF/)
        !  write(*,*) 'RUNSCAN1D ', iscan,TRACING_INT,TRACING_dINT
        !  read(*,*)
        case(rset_bpar)
          filt=NINT(getCmdParam('BEAM1D','FILT'))
          call GET_EVAL_RANGE(SCAN1D%EVAL, XLIM)
          call BEAMPAR(SCAN1D%EVAL%IV,show_monitor,filt,XLIM,PAR,dPAR)
          SCAN1D%VALUES(1:3,iscan)=PAR(1:3)
          SCAN1D%ERRORS(1:3,iscan)=dPAR(1:3)
        case(rset_dpar)
          filt=NINT(getCmdParam('DET1D','FILT'))
          call GET_EVAL_RANGE(SCAN1D%EVAL, XLIM)
          !write(*,1) 'RUNSCAN1D ', XLIM
          call BEAMPAR(SCAN1D%EVAL%IV,show_detector,filt,XLIM,PAR,dPAR)
          !write(*,1) 'RUNSCAN1D ', show_detector, PAR(3), dPAR(3)
          SCAN1D%VALUES(1:3,iscan)=PAR(1:3)
          SCAN1D%ERRORS(1:3,iscan)=dPAR(1:3)
        case(rset_rpar)
          call RESOLPAR(SCAN1D%EVAL%IV,PAR,dPAR)
          SCAN1D%VALUES(1:3,iscan)=PAR(1:3)
          SCAN1D%ERRORS(1:3,iscan)=dPAR(1:3)
        case(rset_bpeak)
          call GET_EVAL_RANGE(SCAN1D%EVAL,XLIM)
          call DOBEAM1D(SCAN1D%EVAL%IV,show_monitor,XLIM,SCAN1D%EVAL%N,RES)
          call setCmdParam('PLOT','TYPE',1.D0*pl_BEAM1D)
          call DOPLOT(.false.)
          do i=1,SCAN1D%EVAL%N
            SCAN1D%VALUES(i,iscan)=RES%Y(i)
            SCAN1D%ERRORS(i,iscan)=RES%DY(i)
          enddo
        case(rset_dpeak)
          call GET_EVAL_RANGE(SCAN1D%EVAL,XLIM)
          call DOBEAM1D(SCAN1D%EVAL%IV,show_detector,XLIM,SCAN1D%EVAL%N,RES)
          call setCmdParam('PLOT','TYPE',1.D0*pl_DET1D)
          call DOPLOT(.false.)
          do i=1,SCAN1D%EVAL%N
            SCAN1D%VALUES(i,iscan)=RES%Y(i)
            SCAN1D%ERRORS(i,iscan)=RES%DY(i)
          enddo
        case(rset_rpeak)
          call DORESOL1D(SCAN1D%EVAL%IV,1.D0*SCAN1D%EVAL%V0,1.D0*SCAN1D%EVAL%DV,SCAN1D%EVAL%N,RES)
          call setCmdParam('PLOT','TYPE',1.D0*pl_RES1D)
          call DOPLOT(.false.)
          do i=1,SCAN1D%EVAL%N
            SCAN1D%VALUES(i,iscan)=RES%Y(i)
            SCAN1D%ERRORS(i,iscan)=RES%DY(i)
          enddo
        case(rset_fmerit)
          SCAN1D%VALUES(1:3,iscan)=PAR(1:3)
          SCAN1D%ERRORS(1:3,iscan)=dPAR(1:3)
        end select
        call INT2STR(iscan,CNUM0)
        call FLOAT2STR(1.D0*SCAN1D%VARIABLES(SCAN1D%IX,iscan),CNUM1)
        call FLOAT2STR(1.D0*SCAN1D%VALUES(SCAN1D%IV,iscan),CNUM2)
        write(*,"(a10,a20,a20)") trim(CNUM0),trim(CNUM1),trim(CNUM2)
        if (showProgress) call XML_PROGRESS(SMES,1,SCAN1D%NS,iscan,'1D scan')
      enddo
  ! restore reporting options
      REPOPT=OPT
      write(*,*) 'SCAN1D finished'
      if (showProgress) call XML_PROGRESS(SMES,2,SCAN1D%NS,iscan,'1D scan')
      end SUBROUTINE RUNSCAN1D

!----------------------------------------------------------------------
      SUBROUTINE RUNSCAN2D
! Run 2D scan with current settings
!----------------------------------------------------------------------
      REAL :: DY
      integer :: ix,iy,ires
      character*32 :: CNUM
  ! cycle y-variable
      DY=(SCAN2D%YMAX-SCAN2D%YMIN)/(SCAN2D%NY-1)
      call XML_PROGRESS(SMES,0,SCAN2D%NY,0,'2D scan')
      do iy=1,SCAN2D%NY
        call FLOAT2STR(1.D0*(SCAN2D%YMIN+(iy-1)*DY),CNUM)
        ires=DLG_PARAM(trim(SCAN2D%HYVAR)//' '//trim(CNUM),.true.)
        if (ires.le.0) exit
      ! copy X-scan parameters to SCAN1D
        call Scan1DAssignX
      ! do 1D scan
        call RUNSCAN1D(.false.)
      ! copy results from SCAN1D
        do IX=1,SCAN2D%NX
      !    write(*,*) 'RUNSCAN2D ',IX,IY,SCAN1D%VALUES(SCAN1D%IV,IX)
      !    read(*,*)
          SCAN2D%VALUES(IX,IY)=SCAN1D%VALUES(SCAN1D%IV,IX)
        enddo
        call XML_PROGRESS(SMES,1,SCAN2D%NY,iy,'2D scan')
      enddo
      call XML_PROGRESS(SMES,2,SCAN2D%NY,iy,'2D scan')
      end SUBROUTINE RUNSCAN2D

!----------------------------------------------------------------------
      SUBROUTINE SAVESCAN1D(fname)
! Save 1D scan to a file
!----------------------------------------------------------------------
      character(*),intent(in) :: fname
      integer :: IU,IRES
      character(MAX_FNAME_LENGTH) :: FRES
2     format('CFG="',a,'"')

      call DLG_FILEOPEN(fname,OUTPATH,'dat',0,1,IRES,FRES)
      if (IRES<=0) return
      IU=OPENFILEUNIT(trim(FRES),.false.)
      if (IU.le.0) return
      write(IU,2) trim(XCONFIG_NAME)
      call WriteScan1D(IU)
      close(IU)
      end SUBROUTINE SAVESCAN1D

!----------------------------------------------------------------------
      SUBROUTINE SAVESCAN2D(fname)
! Save 1D scan to a file
!----------------------------------------------------------------------
      character(*),intent(in) :: fname
      integer :: IU,IRES
      character(MAX_FNAME_LENGTH) :: FRES
2     format('CFG="',a,'"')
      call DLG_FILEOPEN(fname,OUTPATH,'dat',0,1,IRES,FRES)
      if (IRES<=0) return
      IU=OPENFILEUNIT(trim(FRES),.false.)
      if (IU.le.0) return
      write(IU,2) trim(XCONFIG_NAME)
      call WriteScan2D(IU)
      close(IU)
      end SUBROUTINE SAVESCAN2D

      end module SIMEXEC