!//////////////////////////////////////////////////////////////////////
!////  $Id: beam2d.f,v 1.19 2019/08/15 15:02:06 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.19 $
!////     $Date: 2019/08/15 15:02:06 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for 2D beam profiles
!////
!//////////////////////////////////////////////////////////////////////
      MODULE BEAM2D
      use RESULTS
      use NSTORE
      USE TRACINGDATA
      use BEAM1D
      use COMPONENTS
      use EVENTMONITOR
      implicit none
      SAVE
      PRIVATE

      public GETRANGE2D,DOBEAM2D,DETRANGE2D,DOGAUGE2D

      contains

!---------------------------------------------------
      subroutine GETRANGE2D(IMON,IX,IY,FILT,X0,DX,Y0,DY)
! get range of a phase-space variable given by index
!---------------------------------------------------
      integer,intent(in) :: IMON,IX,IY,FILT
      real(kind(1.D0)),intent(out) :: X0,DX,Y0,DY
      real(kind(1.D0)) :: xmin,xmax,xval,ymin,ymax,yval,K0
      integer :: i,NCNT
      type(NEUTRON) :: NEU
      integer,parameter :: MC=10
      !NS=CSUBSET(ISTACK)
      !ND=DSET_ISEL
      integer :: ISTORE
      ISTORE=BMONITORS(IMON)%ISTORE
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      ymin=1.D20
      ymax=-1.D20
      X0=0.D0
      Y0=0.D0
      DX=10.D0
      DY=10.D0
      K0=NSTORE_GETK0(ISTORE)
      IF (NCNT.GT.MC) THEN
        DO I=1,NCNT
          call BMONITORS_GETEVENT(IMON,i,FILT,NEU)
          if (NEU%P>0.D0) then
            XVAL=GETVAR(IX,NEU,K0)
            xmin=min(xmin,XVAL)
            xmax=max(xmax,XVAL)
            YVAL=GETVAR(IY,NEU,K0)
            ymin=min(ymin,YVAL)
            ymax=max(ymax,YVAL)
          endif
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
        Y0=(ymin+ymax)/2.D0
        DY=(ymax-ymin)
      endif
      end subroutine GETRANGE2D

!---------------------------------------------------
      subroutine DETRANGE2D(IX,IY,DET,X0,DX,Y0,DY)
! get range of a phase-space variable given by index
!---------------------------------------------------
      integer,intent(in) :: IX,IY
      TYPE(DETECTOR),intent(in) :: DET
      real(kind(1.D0)),intent(out) :: X0,DX,Y0,DY
      real(kind(1.D0)) :: xmin,xmax,xval,ymin,ymax,yval,K0
      integer :: i,NCNT
      type(NEUTRON) :: NEU
      integer,parameter :: MC=10
      logical :: VALID
      integer :: ISTORE
      ISTORE=DET%FRAME%REGISTRY
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      ymin=1.D20
      ymax=-1.D20
      X0=0.D0
      Y0=0.D0
      DX=10.D0
      DY=10.D0
      K0=NSTORE_GETK0(ISTORE)
      IF (NCNT.GT.MC) THEN
        DO I=1,NCNT
          call NSTORE_GETEVENT(ISTORE,i,NEU)
          XVAL=DETVAR(IX,NEU,DET,VALID)
          if (VALID.and.(NEU%P>0.D0)) then
            xmin=min(xmin,XVAL)
            xmax=max(xmax,XVAL)
            YVAL=DETVAR(IY,NEU,DET,VALID)
            if (VALID) then
              ymin=min(ymin,YVAL)
              ymax=max(ymax,YVAL)
            endif
          endif
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
        Y0=(ymin+ymax)/2.D0
        DY=(ymax-ymin)
      endif
      end subroutine DETRANGE2D

!---------------------------------------------------------
      SUBROUTINE DOBEAM2D(IX,IY,what,X0,Y0,DX,DY,NX,NY,RES,IOUT)
! Fill in 2D histogram in the projection given by the IX,IY indices.
! what=0 ... beam profile
! what=1 ... detector data
! X0,Y0,DX,DY ... center and range of the plot
! NX,NY ... number of pixels
! RES ... structure to store the histogram in
! IOUT ... if > 0, dump event coordinates in this unit
!---------------------------------------------------------
      integer,intent(in) :: IX,IY,NX,NY,what,IOUT
      REAL(KIND(1.D0)),intent(in) :: X0,Y0,DX,DY
      TYPE(TRES2D) :: RES
      REAL(KIND(1.D0)) ::  xmax,xmin,ddx,xval,ymax,ymin,ddy,yval,CNORM,K0,PP
      integer :: I,JX,JY,NCNT,NNX,NNY !,NS,ND
      type(NEUTRON) :: NEU
      TYPE(PDETECTOR) :: DET
      logical :: VALID
      integer :: ISTORE,IMON
      !logical:: dbg,LGX,LGY
      !idbg=0
1     format(3(G12.5,2x))
2     format(a,6(2x,G12.5))
      select case(what)
      case(show_monitor)
        IMON=BMONITOR_REF
        ISTORE=BMONITORS(IMON)%ISTORE
      case(show_detector)
        if (GetDetector(DET)) then
          ISTORE=DET%X%FRAME%REGISTRY
        else
          return
        endif
      end select
      !ND=DSET_ISEL
      !NS=CSUBSET(IBEAM)
      CNORM=NSTORE_GETNORM(ISTORE)
      NNX=MIN(NX,MAX2DIM)
      NNY=MIN(NY,MAX2DIM)
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      K0=NSTORE_GETK0(ISTORE)
      IF (NCNT.GT.10) THEN
        xmin=X0-DX/2
        xmax=X0+DX/2
        ymin=Y0-DY/2
        ymax=Y0+DY/2
        ddx=(xmax-xmin)/NNX
        ddy=(ymax-ymin)/NNY
        if (.not.RES%APPEND) RES%Z=0.D0
        !idbg=0
        !dbg=.false.
        DO I=1,NCNT
          VALID=.true.
          if (what==show_monitor) then
            call BMONITORS_GETEVENT(IMON,i,RES%FILT,NEU)
            XVAL=GETVAR(IX,NEU,K0)
            YVAL=GETVAR(IY,NEU,K0)
          else if (what==show_detector) then
            call NSTORE_GETEVENT(ISTORE,i,NEU)
            XVAL=DETVAR(IX,NEU,DET%X,VALID)
            if (VALID) YVAL=DETVAR(IY,NEU,DET%X,VALID)
          endif
          if (VALID) then
            JX=INT((XVAL-xmin)/ddx)+1
            JY=INT((YVAL-ymin)/ddy)+1
            IF ((JX.GE.1).AND.(JX.LE.NNX)) THEN
              IF ((JY.GE.1).AND.(JY.LE.NNY)) THEN
                PP=CNORM*NEU%P
                RES%Z(JX,JY)=RES%Z(JX,JY)+PP
                if (IOUT>0) then
                  write (IOUT,1) XVAL,YVAL,PP
                endif
              ENDIF
            ENDIF
          endif
        enddo
        RES%NX=NNX
        RES%NY=NNY
      else
        RES%NX=0
        RES%NY=0
      ENDIF
      END SUBROUTINE DOBEAM2D

!---------------------------------------------------------
      SUBROUTINE DOGAUGE2D(IX,IY,ivar,X0,Y0,DX,DY,NX,NY,RES,IOUT)
! Prepare data for gauge volume plot.
! Fill in 2D histogram in the projection given by the IX,IY indices.
! ivar=0 ... intensity
! ivar>0 ... detector variable
! X0,Y0,DX,DY ... center and range of the plot
! NX,NY ... number of pixels
! RES ... structure to store the histogram in
! IOUT ... if > 0, dump event coordinates in this unit
!---------------------------------------------------------
      integer,intent(in) :: IX,IY,NX,NY,ivar,IOUT
      REAL(KIND(1.D0)),intent(in) :: X0,Y0,DX,DY
      TYPE(TRES2D) :: RES
      REAL :: P0(MAX2DIM,MAX2DIM)
      REAL(KIND(1.D0)) ::  xmax,xmin,ddx,xval,ymax,ymin,ddy,yval,CNORM,K0,PP
      integer :: I,JX,JY,NCNT,NNX,NNY
      type(PDATASET) :: DTS
      type(NEUTRON) :: NEU,NEU1,NEU2
      TYPE(PDETECTOR) :: DET
      logical :: VALID
      integer :: ISTORE
      character*32 :: SS

1     format(4(G12.5,2x))
2     format(a,7(2x,G12.5))
3     format(11(G12.5,1x), G13.6)

      !write(*,2) 'DOGAUGE2D: ',IX,IY,ivar,X0-DX/2,X0+DX/2,Y0-DY/2,Y0+DY/2

      ! get references for ki, kf and detector event registers
      if (GetDetector(DET)) then
        ISTORE=DET%X%FRAME%REGISTRY
        if (ISTORE*BMONITOR_PRI*BMONITOR_SEC==0) then
          if (ISTORE<=0) call MSG_WARN('No detector event storage registered.',1)
          if (BMONITOR_PRI<=0) call MSG_WARN('No sample ki-event registered.',1)
          if (BMONITOR_SEC<=0) call MSG_WARN('No sample kf-event registered.',1)
          return
        endif
      else
        ! should not happen, filtered out by calling procedure VIEWBEAM2D
        return
      endif
      CNORM=NSTORE_GETNORM(ISTORE)
      NNX=MIN(NX,MAX2DIM)
      NNY=MIN(NY,MAX2DIM)
      NCNT=NSTORE_MAXN(ISTORE) ! get number of events NCNT
      K0=NSTORE_GETK0(ISTORE)  ! reference K0 (nominal KF)
      IF (NCNT.GT.10) THEN
        xmin=X0-DX/2
        xmax=X0+DX/2
        ymin=Y0-DY/2
        ymax=Y0+DY/2
        ddx=(xmax-xmin)/NNX
        ddy=(ymax-ymin)/NNY
        if (.not.RES%APPEND) then
          RES%Z=0.0
          P0=0.0
        endif
        RES%NX=NNX
        RES%NY=NNY
        DTS=DSET_DATA(DSET_ISEL)
        !write(*,2) 'DOGAUGE2D: start cycle, ',NCNT, ISTORE, BMONITOR_PRI, BMONITOR_SEC, K0
        !write(*,2) 'registeres: ',DTS%P%REG1,DTS%P%REG2
        if (IOUT>0) then
          if (ivar>0) then
            call DVAR_GETNAME(ivar,.true.,SS)
            write(IOUT,2) '# columns: x, y, z, ki_x, ki_y, ki_z, , kf_x, kf_y, kf_z, p, '//trim(SS)
          else
            write(IOUT,2) '# columns: x, y, z, ki_x, ki_y, ki_z, , kf_x, kf_y, kf_z, p'
          endif
        endif
        DO I=1,NCNT
          VALID=.true.
          ! get event from the kf monitor
          call BMONITORS_GETEVENT(BMONITOR_SEC,i,RES%FILT,NEU2)
          !write(*,2) '   NEU2 ',I,NEU2%R(1),NEU2%R(3),NEU2%REF(1),NEU2%REF(2)
          ! get event from the ki monitor
          call NSTORE_GETEVENT(NEU2%REF(1),NEU2%REF(2),NEU1)
          !write(*,2) '   NEU1 ',NEU1%R(1),NEU1%R(3)
          ! get event from the detector
          call NSTORE_GETEVENT(ISTORE,i,NEU)
          !write(*,2) '   NEU ',NEU%R(1),NEU%R(3)
          ! calculatex x,y coordinates
          XVAL=GETVAR(IX,NEU2,K0)
          YVAL=GETVAR(IY,NEU2,K0)
          JX=INT((XVAL-xmin)/ddx)+1
          JY=INT((YVAL-ymin)/ddy)+1
          VALID=((JX.GE.1).AND.(JX.LE.NNX).and.(JY.GE.1).AND.(JY.LE.NNY))
          !write(*,2) '   event ',I,JX,JY,XVAL,YVAL, VALID
          if (VALID) then
            ! monitor detector coordinate (e.g. dhkl)
            if (ivar>0) then
              PP=DETVAR(ivar,NEU,DET%X,VALID)
              if (VALID) then
                RES%Z(JX,JY)=RES%Z(JX,JY)+CNORM*NEU%P*PP
                P0(JX,JY)=P0(JX,JY)+CNORM*NEU%P
                if (IOUT>0) then
                  write (IOUT,3) I, NEU2%R(1:3),NEU1%K(1:3),NEU2%K(1:3),CNORM*NEU%P, PP
                endif
                !write(*,2) '   ',PP,RES%Z(JX,JY),P0(JX,JY)
              endif
            ! monitor intensity
            else
              PP=CNORM*NEU%P
              RES%Z(JX,JY)=RES%Z(JX,JY)+PP
              if (IOUT>0) then
                write (IOUT,3) I, NEU2%R(1:3),NEU1%K(1:3),NEU2%K(1:3), PP
              endif
              !write(*,2) '   ',PP,RES%Z(JX,JY)
            endif

          endif
        enddo
        ! normalize to get the mean values
        if (ivar>0) then
          do JX=1,NNX
            DO JY=1,NNY
              if (P0(JX,JY).ne.0.0) RES%Z(JX,JY)=RES%Z(JX,JY)/P0(JX,JY)
            enddo
          enddo
        endif
      else
        RES%NX=0
        RES%NY=0
      ENDIF
      !write(*,2) 'DOGAUGE2D: done.'
      END SUBROUTINE DOGAUGE2D


      END MODULE BEAM2D
