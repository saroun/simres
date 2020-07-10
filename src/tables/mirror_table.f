!C//////////////////////////////////////////////////////////////////////
!C////  $Id: mirror_table.f,v 1.35 2019/08/16 17:16:27 saroun Exp $
!C////
!C////  R E S T R A X - Simulation of neutron three-axis spectrometers
!C////
!C////     Copyright (C) 1995-2006, All rights reserved
!C////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!C////     Institut Laue Langevin, Grenoble, France
!C////
!C////     Written by:  Jan Saroun
!C////     $Revision: 1.35 $
!C////     $Date: 2019/08/16 17:16:27 $
!C//////////////////////////////////////////////////////////////////////
!C////
!C////  N E S S - Ray-tracing package for  RESTRAX
!C////
!C////  Table for mirror reflectivities
!C////
!C//////////////////////////////////////////////////////////////////////
      module MIRROR_TABLE
      use CONSTANTS
      use XMLINFO
      use RNDGEN
      use SIMATH
      use NITI
      IMPLICIT NONE
      private
      save

      ! defined parameter names
      character(31),parameter :: DEFMIRROR='THICK:RHO:SIG1:SIG2:MC:NITI:COH'
      integer, parameter :: MIRROR_DIM=25   ! max. number of tables
      INTEGER,parameter :: m_nrow=256 ! max. number of table rows

      integer :: m_num=0          ! number of actually used tables
      integer :: m_r(MIRROR_DIM)  ! m-value given as integer = nint(100*m)
      INTEGER :: m_n(MIRROR_DIM)  ! number of rows for each table
      character(512) :: m_files(MIRROR_DIM) ! filenames
      real(kind(1.D0)) :: m_alpha(m_nrow,MIRROR_DIM)  ! angle in nat. Ni
      real(kind(1.D0)) :: m_ref1(m_nrow,MIRROR_DIM)   ! reflectivity spin up
      real(kind(1.D0)) :: m_ref2(m_nrow,MIRROR_DIM)   ! reflectivity spin down

! data allowing for transmission calculations
! NOTE: only m_coating(5) and m_coating(7) is used. The rest was overridden by the m_niti model.
      real(kind(1.D0)) :: m_coating(7,MIRROR_DIM)

! data for niti layers, overrides the above settings, used by the NITI module
! includes (d_Ni, d_Ti, d_all, siga, sigi, sigs)
      real(kind(1.D0)) :: m_niti(6,MIRROR_DIM)


      ! storage for random points on a plane with normal distribution of radii
      ! used to add waviness to mirrors
      integer :: wav_np=0
      integer :: wav_idx=0
      integer,parameter :: wav_max_np=10000
      real(kind(1.D0)) :: wav_r(2,wav_max_np)

      public MIRROR_CLEARALL,MIRROR_GETTABLE,MIRROR_READ
      public MIRROR_REF,MIRROR_REF_EX
      public SET_MTABLE,MIRROR_WRITEXML,MIRROR_WAV_CLEAR
      public MIRROR_REFLECT,MIRROR_GET_PARAM
      public CalReflectivity, WriteReflectivity

      !public MIRROR_REF_EX_DBG

      logical :: mirror_dbg=.false.
      integer :: mirror_idbg

      public mirror_dbg

      contains



!---------------------------------------------------------
      subroutine CalReflectivity(mval, angle, lambda, q, prob, nrows)
! Calculate reflectivity curve.
! Input:
! mval = m-value of the supermirror
! angle = incidence angle [rad]
! lambda = wavelength
! q = array with scattering vectors in Ni_nat. units
! nrows = length of the q and prob arrays
! Output:
! prob(4) = array of probabilities:
!   (1) reflection
!   (2) capture
!   (3) scattering
!   (4) transmission into substrate
!
! if angle>0, calculate at fixed angle, else
! else if lambda>0, calculate at fixed wavelength
! else calculate at fixed wavelength=1.8 A
!--------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: mval, angle, lambda
      REAL(KIND(1.D0)),intent(in) :: q(nrows)
      integer,intent(in) :: nrows
      REAL(KIND(1.D0)),intent(out) :: prob(4,nrows)
      REAL(KIND(1.D0)) :: lam, P(4), QQ
      integer :: i, itab

      itab = MIRROR_GETTABLE(mval)
      if (angle>0) then
        do i=1,nrows
          QQ = q(i)*qTNi
          if (QQ>1.D-3) then
            lam = TWOPI/QQ*angle
          else
            lam = TWOPI/1.D-3*angle
          endif
          call MIRROR_GET_PROB(itab,.true.,QQ, 0.D0, mval, lam, P)
          prob(1:4,i) = P(1:4)
        enddo
      else
        if (lambda>0) then
          lam = lambda
        else
          lam = 1.8D0
        endif
        do i=1,nrows
          QQ = q(i)*qTNi
          call MIRROR_GET_PROB(itab,.true.,QQ, 0.D0, mval, lam, P)
          prob(1:4,i) = P(1:4)
        enddo
      endif
      end subroutine CalReflectivity

!---------------------------------------------------------
      subroutine WriteReflectivity(mval, angle, lambda, IU)
! write reflectivity data to given unit
! Input:
! mval = m-value of the supermirror
! angle = incidence angle [rad]
! lambda = wavelength
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: mval, angle, lambda
      integer,intent(in) :: IU
      integer :: i,nrows, itab
      real(kind(1.D0)) :: prob(4, 256), q(256), dq
      character(LEN_LINE) :: FN
8     FORMAT(a)
9     format('# m=',G12.5, ', angle=',G12.5,', lambda=', G12.5)
10    format(4(G12.5,1x),1x,G12.5)
      nrows = 256
      if (mval<=6.D0) then
        dq = 0.025
      else
        dq = 0.05
      endif
      do i=1,nrows
        q(i) = (i-1)*dq
      enddo
      itab = MIRROR_GETTABLE(mval)
      FN = ' '
      call MIRROR_GETFILE(itab, FN)
      call CalReflectivity(mval, angle, lambda, q, prob, nrows)
      write(IU,8) '# Event probabilities (reflection, capture, scattering, transmission)'
      write(IU,8) '# Table: '//trim(FN)
      write(IU,9) mval, angle, lambda
      write(IU,8) '# columns: q, PR, PC, PS, PT'
      do i=1,nrows
        write(IU,10) q(i), prob(1:4,i)
      enddo
      end subroutine WriteReflectivity

!---------------------------------------------------------
      REAL(KIND(1.D0)) function WAVEP0(theta,a)
! weighting function for importance sampling of waviness angles
! a = random angle generated with normal distribution
! theta = incidence angle
! All angles are in relative units, assuming waviness variance = 1
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: a,theta
      REAL(KIND(1.D0)), parameter :: RTPIHALF=0.5D0*SQRT2PI
      REAL(KIND(1.D0)), parameter :: RT2I=0.707106781186D0
      REAL(KIND(1.D0)) :: th1,dnom,g0,erfmin,erfp1,gn
1     format(a,': ',6(1x,G12.5))
      WAVEP0=0.D0
      th1=theta+a
      if (th1.lt.0) return
      if (theta<-7.0) then
        gn = abs(theta)*(1.D0 +2.D0/theta**2)/SQRT2PI
        if (mirror_dbg) write(*,1) 'WAVEP0   ',theta,gn
      else
        erfp1 = SERF(RT2I*theta)+1.D0
        dnom=exp(-0.5D0*theta*theta) + RTPIHALF*theta*erfp1
        g0=0.5D0*erfp1
        gn = g0/dnom
        if (mirror_dbg) write(*,1) 'WAVEP0   ',th1,erfp1,dnom,g0,gn
      endif
      WAVEP0=SQRT2PI*gn*th1

      end function WAVEP0

!---------------------------------------------------------
      REAL(KIND(1.D0)) function WAVEP01(theta,a)
! weighting function for importance sampling of waviness angles
! a = random angle generated with normal distribution
! theta = incidence angle
! All angles are in relative units, assuming waviness variance = 1
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: a,theta
      REAL(KIND(1.D0)), parameter :: RTPIHALF=0.5D0*SQRT2PI
      REAL(KIND(1.D0)), parameter :: RT2I=0.707106781186D0
      REAL(KIND(1.D0)) :: th1,dnom,g0,erfmin,x
1     format(a,': ',6(1x,G12.5))
      WAVEP01=0.D0
      th1=theta+2.D0*a
      if (th1.lt.0) return
      erfmin=SERF(0.5D0*RT2I*theta)
      x=0.5D0*(theta-1.D0)
      dnom=1+erfmin-(1.D0+SERF(RT2I*x))*exp(-x)
      g0=0.5D0*(1.D0 + erfmin)
      WAVEP01=SQRT2PI*g0*2*(1.D0-exp(-0.5D0*th1))/dnom
      if (mirror_dbg) write(*,1) 'WAVEP01   ',th1,erfmin,g0,SQRT2PI*th1/dnom
      end function WAVEP01


!---------------------------------------------------------
      subroutine WAVE_GEN_ANGLE(theta,a,w)
! Generates waviness angle for given incidence angle
! using importance sampling method
! Takes into account surface illumination and (approximately)
! multiple reflection effects
! theta = incidence angle
! a = random waviness angle
! w = corresponding weight
! All angles are in relative units, assuming waviness variance = 1
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: theta
      REAL(KIND(1.D0)),intent(out) :: a,w
      integer :: loop
      REAL(KIND(1.D0)) ::  th0,th1
      REAL(KIND(1.D0)), parameter ::  inf=5.D0
      REAL(KIND(1.D0)), parameter ::  eps=1.D-6
1     format(a,': ',6(1x,G12.5))
      th0=theta
      a=GASDEV2(-theta)
      w=WAVEP0(theta,a)
      if (mirror_dbg) write(*,1) 'WAVE_GEN_ANGLE',th0,a,w
      if (w<=eps) return
      loop=0
      ! exit angle
      th1 = theta + 2*a
      if (th1.lt.eps) then
        do while ((th1.le.0.D0).and.(loop<10).and.(w>eps))
          ! next incidence angle
          th0=-th1
          a=GASDEV2(-th0)
          w=w*WAVEP0(th0,a)
          ! new exit angle
          th1 = th0 + 2*a
          loop=loop+1
        enddo
        a = 0.5D0*(th1-theta)
        if (loop==10) w=0.D0

      endif
      if (mirror_dbg) write(*,1) '    end with',loop, a, w
      end subroutine WAVE_GEN_ANGLE


!---------------------------------------------------------
      subroutine WAVE_GEN_ANGLE1(theta,a,w)
! Generates waviness angle for given incidence angle
! using importance sampling method
! Takes into account surface illumination and (approximately)
! multiple reflection effects
! theta = incidence angle
! a = random waviness angle
! w = corresponding weight
! All angles are in relative units, assuming waviness variance = 1
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: theta
      REAL(KIND(1.D0)),intent(out) :: a,w
      REAL(KIND(1.D0)) ::  th0
1     format(a,': ',6(1x,G12.5))
      th0=theta
      a=GASDEV2(-th0/2.D0)
      w=WAVEP01(th0,a)
      if (mirror_dbg) write(*,1) 'WAVE_GEN_ANGLE1',th0,a,w
      if (w<=0.D0) return
      end subroutine WAVE_GEN_ANGLE1

!---------------------------------------------------------
      subroutine MIRROR_REFLECT(WAV,K,N,IC,KG,NW,QW,P)
! Calculate Q and K+Q for given surface normal N
! WAV = waviness (sigma) in rad
! K = incident k-vector
! N = surface normal without waviness.
! IC = side index, left,right,top,bottom = 1..4
! Include waviness. Avoid passage under surface.
! RETURN:
! Q = scattering vector
! Return Q<0 if reflection is not possible.
! DN = reflection surface normal (including waviness)
! KG = K + Q*DN
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: WAV,K(3),N(3)
      integer,intent(in) :: IC
      REAL(KIND(1.D0)),intent(out) :: KG(3),QW,NW(3),P
      integer :: i
      REAL(KIND(1.D0)) :: Q,K0,theta,W,KGN,KG0
      integer,parameter :: ID(4)=(/1,1,2,2/)
      integer,parameter :: IDINV(4)=(/2,2,1,1/)
      REAL(KIND(1.D0)),parameter :: THMIN=1.D-6
      REAL(KIND(1.D0)), parameter ::  eps=1.D-6
1     format(a,': ',6(1x,G14.7))
      !mirror_dbg=(mirror_idbg<30)
      !mirror_idbg=mirror_idbg+1

! get scattering vector for a smooth surface , Q=-2(k.n):
      P=0.D0
      Q=0.D0
      K0=0.D0
      DO I=1,3
        Q=Q-K(I)*N(I)
        K0=K0+K(I)**2
      ENDDO
      K0=SQRT(K0)
      theta=Q/K0 ! atan(q, k0) ~ q/k0
      Q=2.D0*Q

      !mirror_dbg=((theta>0.1).or.(theta<0))
      if (mirror_dbg) write(*,1) 'MIRROR_REFLECT ',IC
      if (mirror_dbg) write(*,1) 'N',N(1:3)
      if (mirror_dbg) write(*,1) 'K',K(1:3)
      if (mirror_dbg) write(*,1) 'theta, WAV',theta, WAV

      ! don't allow trajectories at grazing angles less than ~ 0.2 arcsec
      if (theta<THMIN) return

! get normal and Q with waviness
      if (WAV==0.D0) then
        NW=N
        QW=Q
        W=1.D0
      else
        call MIRROR_NWAVY(WAV,theta,K,N,NW,W)
        QW=0.D0
        DO I=1,3
          QW=QW-K(I)*NW(I)
        ENDDO
        theta=QW/K0
        QW=2.D0*QW
        if (mirror_dbg) write(*,1) '  NW,theta,W',NW(1:3),theta,W
      endif
      if (theta<THMIN) return
! get scattered wave vector
      KGN=0.D0
      KG0=0.D0
      do i=1,3
        KG(i)=K(i)+QW*NW(i)
        KGN=KGN+KG(i)*N(i)
        KG0=KG0+KG(i)**2
      enddo
      KG0=SQRT(KG0)
      if (mirror_dbg) write(*,1) '  KGN',KGN,W
      if (KGN<THMIN*K0) return
      ! normalize to K0
      DO i=1,3
        KG(i)=KG(i)*K0/KG0
      enddo
      if (W.lt.eps) W=0.D0
      P=W
      end subroutine MIRROR_REFLECT


!---------------------------------------------------------------
      subroutine MIRROR_NWAVY(sigma,theta,K,N,NW,P)
! Calculate surface normal vector corrected for random waviness
! sigma = waviness in radians
! N .. original surface normal
! K .. incident neutron direction vector
! return
! NW ...  normal with random waviness
! P ... weight factor
! Both N and k must be normalized !!
!---------------------------------------------------------------
      real(kind(1.D0)), intent(in) :: sigma,theta,K(3),N(3)
      real(kind(1.D0)), intent(out) :: NW(3),P
      real(kind(1.D0)) :: a,b,Y(3),Z(3),SG
      integer :: i
1     format(a,': ',6(1x,G14.7))
      ! generate waviness angles
      call WAVE_GEN_ANGLE(theta/sigma,a,P)
      if (mirror_dbg) write(*,1) 'MIRROR_NWAVY theta/sigma,a,P',theta/sigma,a,P
      a=sigma*a
      if (abs(theta+2*a)<1D-6) then
        P=0.D0
        NW=N
        return
      endif
      b=sigma*GASDEV1(0.D0,5.D0)
      ! unit vector along Y (transverse direction)
      ! unit vector along Z (beam propagation)
      ! both normal to N !
      if (mirror_dbg) write(*,1) 'MIRROR_NWAVY a,b, N',a, b, N
      SG=SIGN(1.D0,K(3))
      if (ABS(N(1)).eq.1.D0) then
        Y=(/0.D0, 1.D0, 0.D0/)
        Z=(/0.D0, 0.D0, SG/)
      else if (ABS(N(2)).eq.1.D0) then
        Y=(/1.D0, 0.D0, 0.D0/)
        Z=(/0.D0, 0.D0, SG/)
      else
        call V3cV3(K,N,Y)
        call NORM(Y)
        call V3cV3(N,Y,Z)
      endif
      ! NW = N - a*Z + b*Y
      do i=1,3
        NW(i)=N(i)-a*Z(i)+b*Y(i)
      enddo
      ! we need to normalize due to small-angle approximation
      call NORM(NW)
      end subroutine MIRROR_NWAVY


!---------------------------------------------------------------
      subroutine MIRROR_WAV(sigma,r)
! return coordinates for waviness correction (normalized)
!---------------------------------------------------------------
      real(kind(1.D0)), intent(in) :: sigma
      real(kind(1.D0)), intent(out) :: r(2)
      if (wav_np==0) call WAV_INIT
      if (wav_idx.ge.wav_max_np) wav_idx=0
      wav_idx=wav_idx+1
      r=sigma*wav_r(1:2,wav_idx)
      end subroutine MIRROR_WAV

!---------------------------------------------------------------
      subroutine MIRROR_WAV_N(sigma,N,IX,NW)
! calculate normal vector corrected for waviness
! sigma = waviness in radians
! N .. original surface normal
! IX ... indexes of N axes
! surface normal with random waviness
!---------------------------------------------------------------
      real(kind(1.D0)), intent(in) :: sigma,N(3)
      integer, intent(in) :: IX(3)
      real(kind(1.D0)), intent(out) :: NW(3)
      real(kind(1.D0)) :: a(2)
      call MIRROR_WAV(sigma,a)
      NW(IX(1))=N(IX(1))+N(IX(3))*a(1)
      NW(IX(2))=N(IX(2))+N(IX(3))*a(2)
      NW(IX(3))=N(IX(3))
      call NORM(NW)
      end subroutine MIRROR_WAV_N

!---------------------------------------------------------------
      subroutine MIRROR_WAV_CLEAR
! initialize the wav_r array
!---------------------------------------------------------------
      wav_np=0
      wav_idx=0
      end subroutine MIRROR_WAV_CLEAR

!---------------------------------------------------------------
      subroutine WAV_INIT
! initialize the wav_r array
!---------------------------------------------------------------
      integer :: i
      real(kind(1.D0)) :: x,b
      do i=1,wav_max_np
        x=abs(1.D0*gasdev())
        b=TWOPI*ran1()
        wav_r(1:2,i)=(/x*sin(b),x*cos(b)/)
      enddo
      wav_idx=0
      wav_np=wav_max_np
      end subroutine WAV_INIT


!---------------------------------------------------------------
      subroutine MIRROR_PARAM_GETDEFAULT(par)
! set default parameters of mirror coating
!---------------------------------------------------------------
      real(kind(1.D0)),intent(out) :: par(7)
      ! set defaut values for coating material (NiTi)
      ! thickness [mm]
      par(1)=6e-3
      ! average scattering density [Ang^-2]
      par(2)=4.13e-6
      ! att. cross-sections S1 and S2 (Stot=S0+lambda*S1) [mm-1,mm-1Ang-1]
      par(3)=0.1011
      par(4)=0.0211
      ! mc=0, there is no default
      par(5)=0
      ! using NiTi model
      par(6)=1
      ! accounts for coherent scattering  in NiTi
      par(7)=0

      end subroutine MIRROR_PARAM_GETDEFAULT

!---------------------------------------------------------------
      subroutine MIRROR_GET_PARAM(itab,par,MC)
! MC is the default m-value if there is no table (IC=0)
!---------------------------------------------------------------
      integer,intent(in) :: itab
      real(kind(1.D0)),intent(in) :: MC
      real(kind(1.D0)),intent(out) :: par(7)
      if ((itab>0).and.(itab<=m_num)) then
        par(1:7)=m_coating(1:7,itab)
      else
        call MIRROR_PARAM_GETDEFAULT(par)
        par(5)=MC
      endif
      end subroutine MIRROR_GET_PARAM

!---------------------------------------------------------------
      subroutine MIRROR_GET_PT(itab,MC,K0,QT,PT,PA,PS)
! Get transmission probability for the coating
! for given QT0 (normal k-component) outside
! itab ... table indexe
! MC ... default m-value
! K0 ... 2*PI/wavelength
! QT ... normal k-component of the incident wave
! return:
! PT = probability of transmission through the coating
! PA = probability of capture in the coating
! PS = probability of scattering in the coating
!---------------------------------------------------------------
      integer,intent(in) :: itab
      real(kind(1.D0)),intent(in) :: MC,K0,QT
      real(kind(1.D0)),intent(out) :: PT,PA,PS
      real(kind(1.D0)) :: par(7),MU2PI1,MU2PI2,A,sig1,sig2,sigtot
      call MIRROR_GET_PARAM(itab,par,MC)
      MU2PI1=par(3)*K0  ! scattering term
      MU2PI2=par(4)*TWOPI ! capture term
      ! MU2PI is already in 1/Ang/mm
      ! real and imaginary part of QT^2 - 4 PI rho
      A=QT**2-4.D0*PI*par(2)
      if (A>0.D0) then
        sig1=MU2PI1/SQRT(A)
        sig2=MU2PI2/SQRT(A)
        sigtot=sig1+sig2
      ! transmit
        PT=exp(-sigtot*par(1))
      ! capture
        PA=(1.D0-PT)*sig2/sigtot
      ! scatter
        PS=(1.D0-PT)*sig1/sigtot
      else
        PT=0.D0
        PA=1.D0
        PS=0.D0
      endif
      end subroutine MIRROR_GET_PT

!---------------------------------------------------------------
      subroutine MIRROR_PARAM_SETDEFAULT(ic)
! set default parameters of mirror coating
!---------------------------------------------------------------
      integer,intent(in) :: ic
      real(kind(1.D0)) :: par(7)
      integer :: i
      ! set defaut values for coating material (NiTi)
      call MIRROR_PARAM_GETDEFAULT(par)
      do i=1,7
        m_coating(i,ic)=par(i)
      enddo
      ! derive mc from the table m-value
      m_coating(5,ic)=m_r(ic)*1.D-2
      call MIRROR_WAV_CLEAR
      end subroutine MIRROR_PARAM_SETDEFAULT

!---------------------------------------------------------------
      subroutine MIRROR_CLEARALL
! clear all tables
!---------------------------------------------------------------
      integer :: i
      mirror_idbg=0
      m_num=0
      m_n=0
      do i=1,MIRROR_DIM
        m_files(i)=""
        call MIRROR_PARAM_SETDEFAULT(i)
      enddo
      call MIRROR_WAV_CLEAR
      end subroutine MIRROR_CLEARALL

!---------------------------------------------------------------
      subroutine MIRROR_WRITEXML(IU)
! write list of mirrors in XML
!---------------------------------------------------------------
      integer, intent(in) :: IU
      integer :: i
      character(32) :: CNUM1
1     FORMAT('<MIRRORLIST length="',a,'">')
2     FORMAT('<ITEM key="',a,'" >',a,'</ITEM>')
3     FORMAT('</MIRRORLIST>')

      call XML_RSXDUMP(IU,'',1)
      call INT2STR(m_num,CNUM1)
      write(IU,1) trim(CNUM1)
      do i=1,m_num
        call FLOAT2STR(m_r(i)/100.D0,CNUM1)
        write(IU,2) trim(CNUM1),trim(m_files(i))
      enddo
      write(IU,3)
      call XML_RSXDUMP(IU,'',0)
      end subroutine MIRROR_WRITEXML

!---------------------------------------------------------------
      INTEGER FUNCTION MIRROR_GETTABLE(mValue)
! get index to a table corresponding to given m-value
! return 0 if mValue<=0
! return -1 if no such table exists
!---------------------------------------------------------------
      real(kind(1.D0)),intent(in) :: mValue
      integer :: i,m,res
      if (mValue<=0.D0) then
        res=0
      else
        m=NINT(100*mValue)
        i=1
        do while(i<=m_num)
          if (m_r(i)==m) exit
          i=i+1
        enddo
        if (i<=m_num) then
          res=i
        else
          res=-1
        endif
      endif
      MIRROR_GETTABLE=res
      end FUNCTION MIRROR_GETTABLE


!---------------------------------------------------------------
      subroutine MIRROR_GETFILE(itab,FNAME)
! Return file name for given table
!---------------------------------------------------------------
      integer, intent(in) :: itab
      character*(*) :: FNAME
      if (itab>0) then
        FNAME = trim(m_files(itab))
      endif
      end subroutine MIRROR_GETFILE

!---------------------------------------------------------------
      subroutine MIRROR_READ_PARAM(LINE,IC,IERR)
! Read mirror parameter value from a text line, assumed format is ID=value
! ID      parameter name
! LINE    input text
! IC      table index
! return: ID and index of the parameter from DEFMIRROR
! IERR=-1  ... wrong syntax
! IERR=-2  ... wrong number format
! IERR=1   ... unknown ID
!---------------------------------------------------------------
      character(*),intent(in) :: LINE
      integer, intent(in) :: IC
      integer, intent(out) :: IERR
      character*(16) :: ID
      real(kind(1.D0)) :: Z(1)
      integer :: IE,NA,NPAR
      integer :: idx
1     format(a,' : ',3(1x,G12.6))
      ID=''
      IE=1
      idx=0

      call COUNTPAR(DEFMIRROR,':',NPAR)
      !write(*,1) trim(LINE), NPAR
      call READ_ID_R8(LINE,ID,Z,1,NA,IE)
      !write(*,1) trim(ID), Z, NA,IE
      if ((IE==0).and.(NA>0)) then
        call GETPARINDX(DEFMIRROR,':',ID,idx)
        if ((idx>0).and.(idx<=NPAR)) then
          m_coating(idx,IC)=Z(1)
          !write(*,1) 'MIRROR_READ_PARAM '//trim(ID),idx,Z(1)
        endif
      endif
      IERR=IE
      end subroutine MIRROR_READ_PARAM

!---------------------------------------------------------------
      INTEGER FUNCTION MIRROR_READ(FNAME,mValue)
! read reflectivity data for supemirror of given m-value from file
! use FILETOOLS wrappers - getting file data via GUI
! file structure: 1 line header + 3 columns: m, r(up), r(down)
! returns:
! if OK, index to lookup table
! -1, if table is full
! -2, if file format error occures (less than 3 valid lines)
! 0, if m outside valid range (0,10)
!---------------------------------------------------------------
      use FILETOOLS
      use IO
      IMPLICIT NONE
      real(kind(1.D0)),intent(in) :: mValue
      character*(*) :: FNAME
      character(128) :: MSG,LINE
      character(32) :: CNUM1,CNUM2,CNUM3,CNUM4
      INTEGER :: ierr,itab,res,IU,ird,nrow,iline
      real(kind(1.D0)) :: a(3), d(3), sig(3)
1     format(a,7(1x,G12.4))
2     format('MC=',1x,G12.4,'NITI=',1x,G12.4,'COH=',G12.4)
3     FORMAT('Mirror ref. table(',a,'/',a') read for m=',a,', ',a,' lines')
4     FORMAT('Error while reading mirror table ',a,', line ',a)
10    FORMAT(a)
! JS5/19 IU=OPENFILE_READ(FNAME,11)
      IU = OPENFILEUNIT(FNAME,.true.)
      nrow=0
      res=0
      iline=0
      MSG=''
      if ((mValue<=0.D0).or.(mValue>=1.D1)) then
        call FLOAT2STR(mValue,CNUM1)
        MSG='m-value outside range '//trim(CNUM1)
        goto 99
      endif
      itab=MIRROR_GETTABLE(mValue)
      if (itab==0) then
        goto 99
      else if (itab<0) then
        itab=m_num+1
      endif
      if (itab>=MIRROR_DIM) then
        res=-1
        MSG='Table space exceeded, all tables released'
        call MIRROR_CLEARALL
        goto 99
      endif
      m_r(itab)=NINT(100*mValue)
      call MIRROR_PARAM_SETDEFAULT(itab)
      !JS5/19 ird=FILE_READ(IU,LINE)
      Read(IU,10,iostat=IERR,end=88,err=99) LINE
      ird = len_trim(LINE)
      iline=iline+1
      ierr=0
      do while((ird>0).and.(nrow<m_nrow).and.(ierr==0))
        !JS5/19 ird=FILE_READ(IU,LINE)
        Read(IU,10,iostat=IERR,end=88,err=99) LINE
        iline=iline+1
        ! skip comments
        if (LINE(1:1).ne.'#') then
        ! read parameter
          if (INDEX(LINE,'=')>1) then
            call MIRROR_READ_PARAM(LINE,itab,IERR)
        ! read table values
          else
            call READ_ARRAY8(LINE,3,a,ird)
            if (ird>2) then
              nrow=nrow+1
              m_alpha(nrow,itab)=a(1)
              m_ref1(nrow,itab)=a(2)
              m_ref2(nrow,itab)=a(3)
            else if (ird.ne.0) then
              ierr=1
              ! write(*,*) 'error: ',ierr,ird
              ! write(*,*) 'line:  ',trim(LINE)
            endif
          endif
        endif
      enddo
88    m_n(itab)=nrow
! JS5/19       call CLOSEFILE_READ(IU)
      CLOSE(IU)
      if ((ierr==0).and.(nrow>2)) then
    ! write message
        m_num=max(itab,m_num)
        call INT2STR(itab,CNUM1)
        call FLOAT2STR(mValue,CNUM2)
        call INT2STR(nrow,CNUM3)
        call INT2STR(m_num,CNUM4)
        write(*,3) trim(CNUM1),trim(CNUM4),trim(CNUM2),trim(CNUM3)
        write(*,1) 'file '//trim(FNAME)
        res=itab
      else
        res=-2
        m_n(itab)=0
        goto 99
      endif
      MIRROR_READ=res
      if (res>0) then
        call NITI_GETD(mValue, d)
        call NITI_GETSIG(d, sig)
        m_niti(1:3,itab)=d(1:3)
        m_niti(4:6,itab)=sig(1:3)
        !write(*,1) 'MIRROR_READ m= ',mValue
        !write(*,2) m_coating(5:7,res)
        !write(*,1) 'niti_par= ',m_niti(1:6,res)
        !write(*,1) '---------------------'
        !call NITI_SET_COH((m_coating(7,res)>0))
        !coh = NITI_COH(d, 2.D0)
        !write(*,1) 'sigc(2A)=',coh
      endif
      return
! error:
! JS5/19
!99    call CLOSEFILE_READ(IU)
99    CLOSE(IU)
      call INT2STR(iline,CNUM1)
      write(*,4) trim(FNAME),trim(CNUM1)
      if (len_trim(MSG)>0) write(*,1) trim(MSG)
      MIRROR_READ=res

      end FUNCTION MIRROR_READ


!---------------------------------------------------------------
      INTEGER FUNCTION MIRROR_READ_1(FNAME,mValue)
! read reflectivity data for supemirror of given m-value from file
! use FILETOOLS wrappers - getting file data via GUI
! file structure: 1 line header + 3 columns: m, r(up), r(down)
! returns:
! if OK, index to lookup table
! -1, if table is full
! -2, if file format error occures (less than 3 valid lines)
! 0, if m outside valid range (0,10)
!---------------------------------------------------------------
      use FILETOOLS
      use IO
      IMPLICIT NONE
      real(kind(1.D0)),intent(in) :: mValue
      character*(*) :: FNAME
      character(128) :: MSG,LINE
      character(32) :: CNUM1,CNUM2,CNUM3,CNUM4
      INTEGER :: ierr,itab,res,IU,ird,nrow,iline
      real(kind(1.D0)) :: a(3), d(3), sig(3)
1     format(a,7(1x,G12.4))
2     format('MC=',1x,G12.4,'NITI=',1x,G12.4,'COH=',G12.4)
3     FORMAT('Mirror ref. table(',a,'/',a') read for m=',a,', ',a,' lines')
4     FORMAT('Error while reading mirror table ',a,', line ',a)
      IU=OPENFILE_READ(FNAME,11)
      nrow=0
      res=0
      iline=0
      MSG=''
      if ((mValue<=0.D0).or.(mValue>=1.D1)) then
        call FLOAT2STR(mValue,CNUM1)
        MSG='m-value outside range '//trim(CNUM1)
        goto 99
      endif
      itab=MIRROR_GETTABLE(mValue)
      if (itab==0) then
        goto 99
      else if (itab<0) then
        itab=m_num+1
      endif
      if (itab>=MIRROR_DIM) then
        res=-1
        MSG='Table space exceeded, all tables released'
        call MIRROR_CLEARALL
        goto 99
      endif
      m_r(itab)=NINT(100*mValue)
      call MIRROR_PARAM_SETDEFAULT(itab)
      ird=FILE_READ(IU,LINE)
      iline=iline+1
      ierr=0
      do while((ird>0).and.(nrow<m_nrow).and.(ierr==0))
        ird=FILE_READ(IU,LINE)
        iline=iline+1
        ! skip comments
        if (LINE(1:1).ne.'#') then
        ! read parameter
          if (INDEX(LINE,'=')>1) then
            call MIRROR_READ_PARAM(LINE,itab,IERR)
        ! read table values
          else
            call READ_ARRAY8(LINE,3,a,ird)
            if (ird>2) then
              nrow=nrow+1
              m_alpha(nrow,itab)=a(1)
              m_ref1(nrow,itab)=a(2)
              m_ref2(nrow,itab)=a(3)
            else if (ird.ne.0) then
              ierr=1
              ! write(*,*) 'error: ',ierr,ird
              ! write(*,*) 'line:  ',trim(LINE)
            endif
          endif
        endif
      enddo
      m_n(itab)=nrow
      call CLOSEFILE_READ(IU)
      if ((ierr==0).and.(nrow>2)) then
    ! write message
        m_num=max(itab,m_num)
        call INT2STR(itab,CNUM1)
        call FLOAT2STR(mValue,CNUM2)
        call INT2STR(nrow,CNUM3)
        call INT2STR(m_num,CNUM4)
        write(*,3) trim(CNUM1),trim(CNUM4),trim(CNUM2),trim(CNUM3)
        write(*,1) 'file '//trim(FNAME)
        res=itab
      else
        res=-2
        m_n(itab)=0
        goto 99
      endif
      MIRROR_READ_1=res
      if (res>0) then
        call NITI_GETD(mValue, d)
        call NITI_GETSIG(d, sig)
        m_niti(1:3,itab)=d(1:3)
        m_niti(4:6,itab)=sig(1:3)
        write(*,1) 'MIRROR_READ m= ',mValue
        write(*,2) m_coating(5:7,res)
        write(*,1) 'niti_par= ',m_niti(1:6,res)
        !write(*,1) '---------------------'
        !call NITI_SET_COH((m_coating(7,res)>0))
        !coh = NITI_COH(d, 2.D0)
        !write(*,1) 'sigc(2A)=',coh
      endif
      return
! error:
99    call CLOSEFILE_READ(IU)
      call INT2STR(iline,CNUM1)
      write(*,4) trim(FNAME),trim(CNUM1)
      if (len_trim(MSG)>0) write(*,1) trim(MSG)
      MIRROR_READ_1=res

      end FUNCTION MIRROR_READ_1

!-----------------------------------------------------------
      SUBROUTINE MIRROR_REF_EX_OLD(itab,Q,S,MC,K0,PR,PT,PA)
! extended procedure for calculation of mirror reflectivity
! itab .. mirror table index
! Q    .. one half of the normal component of the scattering vector, in [1/A]
! S    .. projection of neutron polarization vector on mirror magnetizaation (-1..1)
! MC   .. critical m-value (used if itab=0)
! K0   .. |k-vector|=2PI/lambda
! RETURNS:
! PR  .. reflectivity
! PT  .. probability of transmission under the coating
! PA  .. probability of capture in the coating
! Then probability of scattering in the coating is 1 - PR - PT - PA
!-----------------------------------------------------------
      integer,intent(in) :: itab
      REAL(kind(1.D0)),intent(in) :: Q,S,MC,K0
      REAL(kind(1.D0)),intent(out) :: PR,PT,PA
      REAL(kind(1.D0)) :: QC,PT0,p0,PS

      PR=MIRROR_REF(itab,MC,Q,S)
      QC=MC*QTNi
      if (MC.le.0.D0) then
        PR=0.D0
        PA=0.D0
        PT=1.D0
        ! PS=0.D0
      else if (Q<QC) then
        call MIRROR_GET_PT(itab,MC,K0,Q,PT0,PA,PS)
      ! below reflectivity edge, derive PT from the reflectivity curve after QC
        PT=MIRROR_REF(itab,MC,2.D0*QC-Q,S)
        p0=(1.D0-PR-PT)/(1.D0-PT0)
        PA=p0*PA
        ! PS=p0*PS
        ! if (PA<0.D0) write(*,*) 'MIRROR_REF_EX Error, ',PA
      else
      ! otherwise calculate PT from the coating material data
        call MIRROR_GET_PT(itab,MC,K0,Q,PT,PA,PS)
        p0=(1.D0-PR)
        PA=p0*PA
        PT=p0*PT
        ! PS=p0*PS
      endif
      end SUBROUTINE MIRROR_REF_EX_OLD


!-----------------------------------------------------------
      SUBROUTINE MIRROR_GET_PROB(itab,isMonitor,Q,S,MC, lambda, prob)
! extended procedure for calculation of mirror reflectivity
! itab .. mirror table index
! isMonitor .. flag to set monitoring on/off (pass in the monitor value for GUIDE and SGUIDE)
! Q    .. one half of the normal component of the scattering vector, in [1/A]
! S    .. projection of neutron polarization vector on mirror magnetizaation (-1..1)
! MC   .. critical m-value (used if itab=0)
! lambda   .. wavelength [A]
! RETURNS:
! prob(4) = array of probabilities:
!   (1) reflection
!   (2) capture
!   (3) scattering
!   (4) transmission into substrate
!-----------------------------------------------------------
      integer,intent(in) :: itab
      logical,intent(in) :: isMonitor
      REAL(kind(1.D0)),intent(in) :: Q,S,MC,lambda
      REAL(kind(1.D0)),intent(out) :: prob(4)
      REAL(kind(1.D0)) :: qT, r, r0, rc, PR, PA, PS, PT
      REAL(kind(1.D0)) :: P(4),d(3),sig(3), mValue, mpar(2)
      logical isc, validTable
1     format(a,7(1x,G12.5))

      validTable = ((itab>0).and.(itab<=m_num))
      if (validTable) then
        isc = (m_coating(7,itab) > 0.1)
        mValue = m_coating(5,itab)
      else
        isc = .false.
        mValue = MC
      endif
      r=MIRROR_REF(itab,mValue,Q,S)
      if (.not. (isMonitor .or. isc)) then
      ! no need to evaluate scattering and capture
        PR =r ! reflection
        PA = 0.D0
        PS = 0.D0
        PT = 1.D0 - r ! transmission
      else
        if (validTable) then
          d(1:3)=m_niti(1:3,itab)
          sig(1:3)=m_niti(4:6,itab)
          r0 = m_ref1(1,itab)
        else
          call NITI_GETD(mValue, d)
          call NITI_GETSIG(d, sig)
          r0 = 0.99
        endif
        mpar = (/mValue, r0/)
        qT = Q/QTNi
        call NITI_SET_MODE(isc, isMonitor)
        call NITI_PROB(qT, lambda, d, r, mpar, rc, P)
        PR = rc ! reflection
        PA = P(1) ! capture
        PS = P(2)+P(3) ! scattering
        PT = P(4) ! transmission into substrate
      endif
      prob = (/PR, PA, PS, PT/)
      end SUBROUTINE MIRROR_GET_PROB



!-----------------------------------------------------------
      SUBROUTINE MIRROR_REF_EX(itab,isMonitor,Q,S,MC,K0,PROC)
! extended procedure for calculation of mirror reflectivity
! itab .. mirror table index
! isMonitor .. flag to set monitoring on/off (pass in the monitor value for GUIDE and SGUIDE)
! Q    .. one half of the normal component of the scattering vector, in [1/A]
! S    .. projection of neutron polarization vector on mirror magnetizaation (-1..1)
! MC   .. critical m-value (used if itab=0)
! K0   .. |k-vector|=2PI/lambda
! RETURNS:
! PROC: what happens with the neutron:
! 0  .. reflection
! 1  .. transmission under the coating
! 2  .. capture in the coating
! 4  .. scattering in the coating
!-----------------------------------------------------------
      integer,intent(in) :: itab
      logical,intent(in) :: isMonitor
      REAL(kind(1.D0)),intent(in) :: Q,S,MC,K0
      integer,intent(out) :: PROC
      REAL(kind(1.D0)) :: SUMA, PR, PA, PS, PT, RN, prob(4)
1     format(a,7(1x,G14.7))

      call MIRROR_GET_PROB(itab,isMonitor,Q,S,MC,TWOPI/K0,prob)
      PR = prob(1)
      PA = prob(2)
      PS = prob(3)
      PT = prob(4)
      SUMA=PR+PA+PS+PT
      if (abs(SUMA-1.D0)>1.D-6) then
        write(*,1) 'MIRROR_REF_EX error: P, SUMA = ',PR, PA, PS, PT, SUMA
      endif
! play roulette to decide about reflection, transmission or absorption by coating
      RN=1.D0*RAN1()

      if (mirror_dbg) then
         write(*,1) 'MIRROR_REF_EX RN, prob = ', RN, PR, PR+PT, PR+PT+PA, SUMA
         write(*,1) '     test: ',(RN<=PR), (RN<PR+PT), (RN<PR+PT+PA)
      endif

      IF (RN<=PR) THEN
      ! reflected
        PROC=0
      else if (RN<PR+PT) then
      ! transmitted
        PROC=1
      else if (RN<PR+PT+PA) then
      ! absorbed in coating
        PROC=2
      else
      ! scattered in coating
        PROC=4
      endif

      end SUBROUTINE MIRROR_REF_EX

!-----------------------------------------------------------
      REAL(kind(1.D0)) FUNCTION MIRROR_REF(itab,MC,Q,S)
! returns reflectivity  for given momentum transfer Q
! itab .. table index
! Q    .. TWOPI/lambda*incident_angle (assume small-angle approx.)
! S    .. projection of neutron polarization vector on mirror magnetizaation (-1..1)
!-----------------------------------------------------------
      integer,intent(in) :: itab
      REAL(kind(1.D0)),intent(in) :: MC,Q,S
      REAL(kind(1.D0)) :: Q1,res,dQ,z,rnd
      integer :: iz,NR

      Q1=Q/qTNi
      ! no table - use analytical form
      if (itab<0) then
        res=MIRROR_REF_ANAL(Q1,MC)
      ! udefined table or MC=0
      else if ((itab==0).or.(itab>m_num)) then
        res=0.D0
      ! Q1 out of table range
      else if ((Q1<m_alpha(1,itab)).or.(Q1>m_alpha(m_n(itab),itab))) then
        res=0.D0
      ! valid input, interpolate in the table
      else
        NR=itab
        dQ=(m_alpha(m_n(itab),NR)-m_alpha(1,NR))/(m_n(NR)-1)
        z=(Q1-m_alpha(1,NR))/dQ
        iz=INT(z)+1
        IF ((z<0).OR.(z>=m_n(NR)).OR.(iz>=m_n(NR))) THEN
          res=0.D0
        else
          rnd=2.D0*ran1()-1.D0
          IF (S>=rnd) THEN
            res=m_ref1(iz,NR)+(z-iz+1)*(m_ref1(iz+1,NR)-m_ref1(iz,NR))
          ELSE
            res=m_ref2(iz,NR)+(z-iz+1)*(m_ref2(iz+1,NR)-m_ref2(iz,NR))
          ENDIF
        endif
      ENDIF
      MIRROR_REF=res
      END FUNCTION MIRROR_REF


!-----------------------------------------------------------
      REAL(kind(1.D0)) FUNCTION MIRROR_REF_ANAL_OLD(Q,mc,w,r1,r2)
! returns reflectivity for analytical curve
! Q  .. incident angle, m-value
! mc .. critical mirror angle [m-value]
! w  .. cut-off width (m-value)
! r1     .. top reflectivity for m<1
! r2     .. edge reflectivity for m=gamma
!-----------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: Q,mc,w,r1,r2
      REAL(kind(1.D0)) :: res,Z
      if ((Q<0.D0).or.(Q>mc+10.D0*w)) then
        res=0.D0
      else if ((Q<=1.D0).or.(mc<=1.D0)) then
        res=r1
      else
        if (w<=1.D-5) then
          Z=-1.D0
        else
          Z=tanh((Q-mc)/w)
        endif
        res=0.5D0*(1.D0-Z)*(r1-(Q-1.D0)/(mc-1.D0)*(r1-r2))
      endif
      MIRROR_REF_ANAL_OLD=res
      END FUNCTION MIRROR_REF_ANAL_OLD

!-----------------------------------------------------------
      REAL(kind(1.D0)) FUNCTION MIRROR_REF_ANAL(Q,mc)
! returns reflectivity for analytical curve
! sets r(q) parameters automatically for given mc
! using fits of average r(mc) reported by SwissNeutronics
! http://www.swissneutronics.ch, status 2016
! Q  .. incident angle, m-value
! mc .. critical mirror angle [m-value]
!-----------------------------------------------------------
      REAL(kind(1.D0)),intent(in) :: Q,mc
      REAL(kind(1.D0)) :: res,Z,r2,mval
      REAL(kind(1.D0)),parameter :: R0=0.99
      REAL(kind(1.D0)),parameter :: W=0.05
      REAL(kind(1.D0)),parameter :: PM(3)=(/1.02D0,-2.86D-2,-6.25D-3/)
      mval=mc+1.5*W
      if ((Q<0.D0).or.(Q>mval+10.D0*W)) then
        res=0.D0
      else if ((Q<=1.D0).or.(mc<=1.D0)) then
        res=R0
      else
        Z=tanh((Q-mval)/W)
        r2=PM(1)+PM(2)*mc+PM(3)*mc**2
        res=0.5D0*(1.D0-Z)*(R0-(Q-1.D0)/(mc-1.D0)*(R0-r2))
      endif
      MIRROR_REF_ANAL=res
      END FUNCTION MIRROR_REF_ANAL

!--------------------------------------------------------------------
      subroutine SET_MTABLE(CMD)
! wrapper for MIRROR_READ function, to be called from command interpreter
! interpret 1st argument as m-value and the 2nd as filename
!--------------------------------------------------------------------
      character(*),intent(in) :: CMD
      real(kind(1.D0)) :: mValue=0.D0
      character*(256) :: FNAME
      integer :: IS,IL,IRES,ierr
      ! write(*,*) 'SET_MTABLE: ',trim(CMD)
      IRES=-3
      IS=1
      CALL FINDPAR(CMD,1,IS,IL)
      if (CMD(IS:IS+IL-1).eq.'CLEAR') then
        call MIRROR_CLEARALL
        return
      else if (CMD(IS:IS+IL-1).eq.'DUMP') then
        ! write(*,*) 'SET_MTABLE ',CMD(IS:IS+IL-1)
        call MIRROR_WRITEXML(smes)
        return
      endif
      read(CMD(IS:IS+IL-1),*,ERR=99,END=10,IOSTAT=ierr) mValue
10    CALL FINDPAR(CMD,2,IS,IL)
      FNAME=CMD(IS:IS+IL-1)

      IRES=MIRROR_READ(FNAME,mValue)
      ! write(*,*) 'SET_MTABLE ires=',ires
      if (IRES>0) then
        m_files(IRES)=trim(FNAME)
        return
      endif

99    select case(IRES)
      case(-3)
        call MSG_ERROR('MIRROR_READ','Wrong format of mirror table key, must be a number',0,1)
      case(-2)
        call MSG_ERROR('MIRROR_READ','Wrong mirror table format: '//trim(FNAME),0,1)
      case(-1)
        call MSG_ERROR('MIRROR_READ','No more space for mirror tables.',0,1)
      case(0)
        call MSG_ERROR('MIRROR_READ','m-value for the table is out of range (0..10)',0,1)
      !case default

      end select
      end subroutine SET_MTABLE

      end  module MIRROR_TABLE

