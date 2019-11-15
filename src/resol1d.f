!//////////////////////////////////////////////////////////////////////
!////  $Id: resol1d.f,v 1.13 2019/08/16 17:16:25 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2009, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.13 $
!////     $Date: 2019/08/16 17:16:25 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Module for 1D resolution functions
!////
!//////////////////////////////////////////////////////////////////////
      MODULE RESOL1D
      use RESULTS
      use NSTORE
      use REPORTSOPT
    !  use COMMANDS
    !  USE COMPONENTS
      USE TRACINGDATA
      use DATASETS
      implicit none
      SAVE
      PRIVATE

      public GETQVAR,GETQRANGE,DORESOL1D,RESOLPAR

      contains

!---------------------------------------------------
      subroutine GETQVARSTR(IVAR,SARG)
! get resolution space variable by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      character(*),intent(out) :: SARG
      character(1),parameter :: CH(3)=(/'X','Y','Z'/)
      if (len(SARG).ge.16) then
        select case(IVAR)
    ! dQ/Q
        case(1)
          SARG='dQ/Q'
    ! theta
        case(2)
          SARG='theta [deg]'
    ! EN
        case(3)
          SARG='EN [meV]'
    ! dQ in C&N coordinates
        case(4,5,6)
          SARG='Q'//CH(IVAR-3)//' [A^-1]'
        end select
      endif
      end subroutine GETQVARSTR

!---------------------------------------------------
      real(KIND(1.D0)) function GETQVAR(IVAR,NEU1,NEU2)
! get resolution space variable by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      type(NEUTRON),intent(in) :: NEU1,NEU2
      real(KIND(1.D0)) :: XVAL,Q(3),s,c,Q0(3),Q00,QCN(3),AUX(3),K01,K02
1     format(a12,': ',6(G10.4,1x),I6)
      XVAL=0.D0
      Q0=SPEC%NOM%KF-(/0.D0,0.D0,SPEC%NOM%KI/)
      Q00=SQRT(Q0(1)**2+Q0(2)**2+Q0(3)**2)
      Q=NEU2%K-NEU1%K
      select case(IVAR)
    ! dQ/Q
        case(1)
          XVAL=(SQRT(Q(1)**2+Q(2)**2+Q(3)**2)-Q00)/Q00
          !write(*,1) 'GETQVAR I,Q0,Q00',IVAR,Q0,Q00
          !write(*,1) 'GETQVAR Q',Q
    ! theta
        case(2)
        !  c=(NEU1%K(1)*NEU2%K(1)+NEU1%K(2)*NEU2%K(2)+NEU1%K(3)*NEU2%K(3))/(NEU1%K0*NEU2%K0)
          K01=NEU1%K(1)**2+NEU1%K(3)**2
          K02=NEU2%K(1)**2+NEU2%K(3)**2
          c=(NEU1%K(1)*NEU2%K(1)+NEU1%K(3)*NEU2%K(3))/SQRT(K01*K02)
          s=(NEU1%K(3)*NEU2%K(1)-NEU1%K(1)*NEU2%K(3))/SQRT(K01*K02)
          XVAL=asin(s)/deg
          if (c.lt.0.D0) XVAL=SIGN(1.D0,s)*180.D0-XVAL

         ! XVAL=acos(c)/deg*SIGN(1.D0,NEU1%K(3)*NEU2%K(1)-NEU1%K(1)*NEU2%K(3))

        !  AUX=SPEC%NOM%KF
        !  K01=AUX(1)**2+AUX(3)**2
        !  K02=NEU2%K(1)**2+NEU2%K(3)**2
        !  c=(AUX(3)*NEU2%K(1)-AUX(1)*NEU2%K(3))/SQRT(K01*K02)
        !  XVAL=asin(c)/deg
    ! EN
        case(3)
          XVAL=HSQOV2M*(NEU1%K0**2-NEU2%K0**2)
    ! dQ in C&N coordinates
        case(4,5,6)
          AUX=Q-Q0
      !    write(*,1) 'GETQVAR PRE',IVAR,AUX(1:3)
          call M3xV3(1,SPEC%NOM%MCL,Q-Q0,QCN) ! convert to C&N
      !    write(*,1) 'GETQVAR POST',IVAR,QCN(1:3)
          XVAL=QCN(IVAR-3)
      end select
      GETQVAR=XVAL
      end function GETQVAR

!---------------------------------------------------
      subroutine GETQRANGE(IVAR,X0,DX)
! get range of a resolution variable given by index
!---------------------------------------------------
      integer,intent(in) :: IVAR
      real(kind(1.D0)),intent(out) :: X0,DX
      real(kind(1.D0)) :: xmin,xmax,xval
      integer :: i,NCNT2,IREF1,NEV1
      type(NEUTRON) :: NEU1,NEU2
      integer,parameter :: MC=10
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
      NCNT2=NSTORE_MAXN(DT%P%REG2) ! get number of events NCNT
      xmin=1.D20
      xmax=-1.D20
      X0=0.D0
      DX=10.D0
      IF (NCNT2.GT.MC) THEN
        DO I=1,NCNT2
          ! call NSTORE_REFERENCE(I,NS2,ND2,IREF)
          ! call NSTORE_GET_KF(I,NS2,ND2,NEU2)
          ! call NSTORE_GET_KI2(IREF,NEU1)
          call NSTORE_GETEVENT(DT%P%REG2,i,NEU2)
          IREF1=NEU2%REF(1)
          NEV1=NEU2%REF(2)
          if ((IREF1==DT%P%REG1).and.(NEV1>0)) then
            call NSTORE_GETEVENT(IREF1,NEV1,NEU1)
            XVAL=GETQVAR(IVAR,NEU1,NEU2)
            xmin=min(xmin,XVAL)
            xmax=max(xmax,XVAL)
          else
            write(*,*) 'Error in GETQRANGE ',IREF1,NEV1
          endif
        enddo
        X0=(xmin+xmax)/2.D0
        DX=(xmax-xmin)
      endif
      end subroutine GETQRANGE

C----------------------------------------------------
      SUBROUTINE DORESOL1D(IVAR,X0,XRANGE,N,RES)
C get resolution profile along 1 dimension
C----------------------------------------------------
      integer,intent(in) :: IVAR,N
      REAL(KIND(1.D0)),intent(in) :: X0,XRANGE
      TYPE(TRES1D) :: RES
      integer,parameter :: MC=10
      real :: EFY(SPCM,MC),Y,DY
      REAL(KIND(1.D0)) ::  xmax,xmin,dx,xval,CNORM
      integer :: I,J,IP,NCNT2,NCNT1,NN !,NS2,ND2
      type(NEUTRON) :: NEU1,NEU2
      REAL(KIND(1.D0)) :: PAR(3),dPAR(3)
      character(32) :: SARG
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
!1     format(a12,': ',6(G10.4,1x))
      !write(*,1) 'DORESOL1D NOM.KI,NOM.KF',SPEC%NOM%KI,SPEC%NOM%KF

      !NS2=CSUBSET(2)
      !ND2=DSET_ISEL
      CNORM=MC*NSTORE_GETNORM(DT%P%REG2)

      NN=MIN(N,SPCM)

      NCNT2=NSTORE_MAXN(DT%P%REG2) ! get number of events NCNT
      NCNT1=NSTORE_MAXN(DT%P%REG1) ! get number of events NCNT
      call CLEAR_RES1D(RES)

      if ((NCNT2<=0).or.(NCNT1<=0).or.(CNORM<=0.D0)) then
        call MSG_WARN('Can''t calculate resolution function, no neutrons stored in the scattering channel.',1)
        return
      endif

      IF (NCNT2.GT.MC) THEN
        EFY=0.0
        xmin=X0-XRANGE/2
        xmax=X0+XRANGE/2
        dx=(xmax-xmin)/NN
        DO I=1,NN
          RES%X(I)=xmin+(I-0.5D0)*dx
          RES%Y(I)=0.
        ENDDO
    ! scan events in STACK2 (final states)
        DO I=1,NCNT2
          ! call NSTORE_REFERENCE(I,NS2,ND2,IREF)
          ! call NSTORE_GET_KF(I,NS2,ND2,NEU2)
          ! call NSTORE_GET_KI2(IREF,NEU1)

          call NSTORE_GETEVENT(DT%P%REG2,i,NEU2)
          call NSTORE_GETEVENT(NEU2%REF(1),NEU2%REF(2),NEU1)
          XVAL=GETQVAR(IVAR,NEU1,NEU2)
      !    write(*,*) 'DORESOL1D IVAR=',IVAR,' CNT=',I
      !    write(*,"(a,3(1x,G10.3))") '   KI=', NEU1.K
      !    write(*,"(a,3(1x,G10.3))") '   KF=', NEU2.K
      !    write(*,"(a,3(1x,G10.3))") '  VAL=', XVAL
      !    read(*,*)
          J=INT((XVAL-xmin)/dx)+1
          IP=MIN(MC,INT(1.0*MC*I/NCNT2)+1)
          IF ((J.GE.1).AND.(J.LE.NN)) THEN
            EFY(J,IP)=EFY(J,IP)+NEU2%P
          ENDIF
        enddo
! renormalize
      !  if (NTRIALS.gt.0) then
          DO I=1,NN
            Y=0.D0
            DY=0.D0
            do IP=1,MC
              Y=Y+EFY(I,IP)
              DY=DY+EFY(I,IP)**2
            enddo
            Y=Y/MC
            DY=DY/MC
            RES%DY(I)=0.D0
            if(DY.gt.Y**2) RES%DY(I)=CNORM*SQRT(DY-Y**2)
            RES%Y(I)=CNORM*Y
          ENDDO
      !  endif
        RES%N=NN
        if (REPOPT%RES) then
          call RESOLPAR(IVAR,PAR,dPAR)
          call XML_RSXDUMP(SOUT,' ',1)
          call GETQVARSTR(IVAR,SARG)
          call XML_TITLE(SOUT,'Resolution function profile for '//trim(SARG))
          call XML_FVALUE(SOUT,'intensity',' ',PAR(1),dPAR(1))
          call XML_FVALUE(SOUT,'width',' ',PAR(2),dPAR(2))
          call XML_RSXDUMP(SOUT,' ',0)
        endif
      else
        RES%N=0
      ENDIF
      END SUBROUTINE DORESOL1D

!----------------------------------------------------
      SUBROUTINE RESOLPAR(IVAR,PAR,dPAR)
! get resolution profile parameters for the IVAR variable
! Similar to RESOL1D, but makes no histogram, just statistics
! the paramers in PAR are:
! Sum,dSum,Width,dWidth,Center,dCenter
!----------------------------------------------------
      integer,intent(in) :: IVAR
      REAL(KIND(1.D0)),intent(out) :: PAR(3),dPAR(3)
      integer,parameter :: MC=10
      REAL(KIND(1.D0)) ::  XVAL,CNORM,S0,S1,S2,C,W,Z
      integer :: I,IM,M(0:MC),NCNT2 !,NS2,ND2
      type(NEUTRON) :: NEU1,NEU2
      type(PDATASET) :: DT
      DT=DSET_DATA(DSET_ISEL)
      !NS2=CSUBSET(2)
      !ND2=DSET_ISEL
      CNORM=MC*NSTORE_GETNORM(DT%P%REG2)
      NCNT2=NSTORE_MAXN(DT%P%REG2) ! get number of events NCNT
      PAR=0.D0
      dPAR=0.D0
      if (NCNT2.le.MC) return
      M(0)=0
      do IM=1,MC-1
        M(IM)=M(IM-1)+NCNT2/MC
      enddo
      M(MC)=NCNT2
      IF (NCNT2.GT.MC) THEN
      !  write(*,*) 'RES1PAR NTRIALS,CNORM= ',NTRIALS,CNORM,NS2,ND2,IVAR
        do IM=1,MC
          S0=0.D0
          S1=0.D0
          S2=0.D0
          DO I=M(IM-1)+1,M(IM)
            ! call NSTORE_REFERENCE(I,NS2,ND2,IREF)
            ! call NSTORE_GET_KF(I,NS2,ND2,NEU2)
            ! call NSTORE_GET_KI2(IREF,NEU1)

            call NSTORE_GETEVENT(DT%P%REG2,i,NEU2)
            call NSTORE_GETEVENT(NEU2%REF(1),NEU2%REF(2),NEU1)
            XVAL=GETQVAR(IVAR,NEU1,NEU2)
            Z=CNORM*NEU2%P
            S0=S0+Z
            S1=S1+Z*XVAL
            S2=S2+Z*XVAL**2
          !  write(*,*) 'RES1PAR I,IREF = ',I,IREF
          !  write(*,*) 'RES1PAR KI = ',NEU1%K
          !  write(*,*) 'RES1PAR KF = ',NEU2%K
          !  write(*,*) 'RES1PAR X,P = ',XVAL,NEU2%P
          !  read(*,*)
          enddo
        !  write(*,*) 'RES1PAR S0,S1,S2= ',S0,S1,S2
          IF(S0.GT.0.D0) then
            C=S1/S0
            W=R8LN2*SQRT(ABS(S2/S0-C**2))
            PAR(1)=PAR(1)+S0
            dPAR(1)=dPAR(1)+S0**2
            PAR(2)=PAR(2)+W
            dPAR(2)=dPAR(2)+W**2
            PAR(3)=PAR(3)+C
            dPAR(3)=dPAR(3)+C**2
          endif
        !  write(*,*) 'RES1PAR PAR= ',PAR
        !  write(*,*) 'RES1PAR dPAR= ',dPAR
        !  read(*,*)
        enddo
        do I=1,3
          PAR(i)=PAR(i)/MC
          dPAR(i)=dPAR(i)/MC
          if (dPAR(i).gt.PAR(i)**2) dPAR(i)=SQRT(dPAR(i)-PAR(i)**2)
        enddo
      ENDIF
      END SUBROUTINE RESOLPAR


      END MODULE RESOL1D
