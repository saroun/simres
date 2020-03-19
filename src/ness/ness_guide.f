!       module GUIDE
!//////////////////////////////////////////////////////////////////////
!////  $Id: ness_guide.f,v 1.52 2019/08/16 17:16:26 saroun Exp $
!////
!////  S I M R E S - Ray-tracing simulation of neutron spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.52 $
!////     $Date: 2019/08/16 17:16:26 $
!//////////////////////////////////////////////////////////////////////
!////
!////  N E S S - Ray-tracing package for RESTRAX
!////
!////  Neutron guide and bender
!////  Can be multi-channeled and curved to one side in both horizontal
!////  and vertical directions. Lamellae are equidistant.
!////  Handles transmitting (and polarizing) lamellae if OBJ%MU.le.1.D5
!//////////////////////////////////////////////////////////////////////

!       contains
!---------------------------------------------------------------------
      SUBROUTINE GDE_INIT(OBJ)
! initializes some calculated parameters and transformation matrices
!----------------------------------------------------------------------
      USE GUIDES
      USE FRAMES_TRACE
      IMPLICIT NONE

      REAL*8 ZL
      TYPE(GUIDE) :: OBJ
      INTEGER*4 I

1     FORMAT(a,$)
      CALL FRAME_INIT(OBJ%FRAME)
! CH,CV are curvatures
      OBJ%CH=OBJ%ROH
      OBJ%CV=OBJ%ROV

! limit for number of slits = 255
      IF (OBJ%NLH.GT.255) OBJ%NLH=255
      IF (OBJ%NLV.GT.255) OBJ%NLV=255
! set A(I) = angles for lamellae
! set L(I) = positions of lamellae at the entry
      IF(OBJ%FRAME%SIZE(3).GT.0) THEN
        DO I=0,OBJ%NLH
          ZL=I*1.D0/OBJ%NLH - 0.5D0
          OBJ%AH(I)=ZL*(OBJ%W2-OBJ%FRAME%SIZE(1))/OBJ%FRAME%SIZE(3)
          OBJ%LH(I)=ZL*OBJ%FRAME%SIZE(1)
        ENDDO
        DO I=0,OBJ%NLV
          ZL=I*1.D0/OBJ%NLV - 0.5D0
          OBJ%AV(I)=(OBJ%H2-OBJ%FRAME%SIZE(2))/OBJ%FRAME%SIZE(3)*ZL
          OBJ%LV(I)=ZL*OBJ%FRAME%SIZE(2)
        ENDDO
      ENDIF

      END SUBROUTINE GDE_INIT

!---------------------------------------------------------------------
      REAL*8 FUNCTION GUIDE_LAM(OBJ,ID,IL,Z,IORI)
! Function describing lamellas profile along guide axis (z-coordinate)
! Returns lateral distance of the IL-th lamella (middle)
! or its derivative along z
! ID   ... derivative
! IL   ... lamella index
! Z    ... distance from guide entry
! IORI  ... horizontal (1) or vertical (2)
!----------------------------------------------------------------------
      USE GUIDES
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ

      REAL*8 a,f,x0
      REAL*8 Z
      INTEGER*4 ID,IL,IORI

      IF (IORI.EQ.2) THEN   ! vertical slit
        f=OBJ%CV
        a=OBJ%AV(IL)
        x0=OBJ%LV(IL)
      ELSE                 ! horizontal slit
        f=OBJ%CH
        a=OBJ%AH(IL)
        x0=OBJ%LH(IL)
      ENDIF

! zero deriv.
      IF (ID.LE.0) THEN
        GUIDE_LAM=x0+a*z+0.5*f*z**2
! 1st deriv.
      ELSE IF (ID.EQ.1) THEN
        GUIDE_LAM=a+f*z
! 2nd deriv.
      ELSE IF (ID.EQ.2) THEN
        GUIDE_LAM=f
      ELSE
        GUIDE_LAM=0.D0
      ENDIF

      END FUNCTION GUIDE_LAM


!----------------------------------------------------------------------------
      REAL*8 FUNCTION GUIDE_CROSS1(kx,kz,x,z,gx,alpha,ro)
! return time of flight to a cross-point with a parabolic surface
! surface is assumed in the form x=0.5*ro*z**2 + alpha*z
! if there are 2 solutions, take the minimum positive one
! return 10^30 if there is no such solution
! alpha  .. tan(angle) between lamella and guide axis)
! ro     .. curvature
! x,z    .. starting point (z // guide axis)
! kx,kz  .. ray direction
!----------------------------------------------------------------------------
      use TRACELIB
      IMPLICIT NONE
      real(kind(1.D0)) :: kx,kz,x,z,alpha,ro,gx
      real(kind(1.D0)) :: t1,t2
      real(kind(1.D0)),parameter :: EPS=1.D-8
      call TLIB_CROSS(kx,kz,x,z,(/0.5D0*ro,alpha,0.D0/),gx,t1,t2)
      if (t1>EPS) then
        GUIDE_CROSS1=t1
      else if (t2>EPS) then
        GUIDE_CROSS1=t2
      else
        GUIDE_CROSS1=TLIB_INF
      endif
      END FUNCTION GUIDE_CROSS1

!---------------------------------------------------------------------
      SUBROUTINE GUIDE_CON(OBJ,R,K,IH,IV,IC,DT,Q,NW,KG,P)
! get transport data inside CHANNEL
! - find next contact with a channel wall
! - moves to the next contact point, or to the exit
! *** INPUT:
! IH,IV   ... channel indices
! R(3)    ... neutron coordinate
! K(3)    ... neutron k-vector
! *** RETURN:
! IC    .. which wall:  left(1), right(2), top(3), bottom (4)
!          if no contact, pass through, IC=0
! DT    .. time to reach the next contact point
! Q     .. reflection vector (magnitude)
! N     .. unit vector normal to the surface
!----------------------------------------------------------------------
      USE GUIDES
      USE TRACINGDATA
      USE CONSTANTS
      use MIRROR_TABLE
      use SIMATH
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      real(kind(1.D0)), intent(in) :: K(3)
      real(kind(1.D0)), intent(out) :: DT,Q,NW(3),KG(3),P
      real(kind(1.D0)), intent(inout) :: R(3)
      INTEGER,intent(in) :: IH,IV
      INTEGER,intent(out) :: IC
      INTEGER :: ISG
      real(kind(1.D0)) :: GET_CROSSTIME,GUIDE_LAM
      real(kind(1.D0)) :: LZ(4),ANG,T0,LPOS,KGR(3),N(3)
      INTEGER :: I,IX,ILAM(4)
      integer,parameter :: IL(4)=(/1,0,1,0/)
      integer,parameter :: ID(4)=(/1,1,2,2/)
      integer,parameter :: IDINV(4)=(/2,2,1,1/)
! shifts for an oscillating collimator
      real(kind(1.D0)) :: OSCD,OSCA
      COMMON /OSCBEND/ OSCD,OSCA

      ! dbg=(trim(OBJ%FRAME%ID).eq.'GDB')

2     format(a,': ',8(G12.5,1x))

! TOF to the guide end
      if (TRACING_UP) then
        T0=-R(3)/K(3)
      else
        T0=(OBJ%FRAME%SIZE(3)-R(3))/K(3)
      endif

! TOF to the next boundary in units h_bar/m=1, i.e. path=k*LZ
      ILAM=(/IH,IH,IV,IV/)
      DO I=1,4
        LZ(I)=GET_CROSSTIME(OBJ,R,K,ILAM(I),I,.FALSE.)
      ENDDO

! select the shortest time to cross a surface
      DT=MIN(LZ(1),LZ(2),LZ(3),LZ(4),T0)
! DT<0 should never happen, but to be sure, ...
      IF (DT.LE.0D0) goto 99 ! signal from GET_CROSSTIME to stop

 ! and find which of the 4 surfaces it was
      IC=0
      DO I=1,4
        IF (DT.EQ.LZ(I)) IC=I
      ENDDO
    ! index of the left, right, top, bottom lamella
      ILAM=(/IH+1,IH,IV+1,IV/)
! if no contact => pass through
      IF (IC.EQ.0) THEN
        DT=T0
        return
      ENDIF
      IX=ID(IC)

! get the surface normal vector at the contact point
      N(1:3)=0.D0
      ISG=1-2*IL(IC) ! /1,0,1,0/ ==> /-1,1,-1,1/
      ANG=GUIDE_LAM(OBJ,1,ILAM(IC),R(3)+DT*K(3),IX)
      LPOS=GUIDE_LAM(OBJ,0,ILAM(IC),R(3)+DT*K(3),IX)
      if (IX==1) then
         LPOS=LPOS + OSCD + (R(3)+DT*K(3))*OSCA
         ANG=ANG+OSCA
      endif
      N(3)=-ISG*ANG/SQRT(1.D0+ANG**2)
      N(IX)=ISG*SQRT(1.D0-N(3)**2)

! stop events reflected from convex side
      if ((OBJ%ONESIDE.GT.0).and.(DT.LT.T0)) then
        select case (IX)
        case(1)
          IF (N(1)*OBJ%ROH.LT.0) goto 99
        case(2)
          IF (N(2)*OBJ%ROV.LT.0) goto 99
        end select
      endif

      ! get K changed by gravity
      do i=1,3
        KGR(i)=K(i)+NEUT_GR(i)*DT
      enddo
      ! write(*,2) 'GUIDE_CON '//trim(OBJ%FRAME%ID),IC,K
      call MIRROR_REFLECT(OBJ%WAV,KGR,N,IC,KG,NW,Q,P)
      if (P<=0.D0) goto 99

      Return
! signal to stop the event
99    NEUT%STATE=1
      IC=0
      END SUBROUTINE GUIDE_CON

!---------------------------------------------------------------------
      SUBROUTINE GUIDE_TRANS(OBJ,R,K,IH,IV,IC,DT,Q,NW,KG,P)
! get transport data inside LAMELLA
! - find next contact with a channel wall
! - moves to the next contact point, or to the exit
! *** INPUT:
! IH,IV   ... channel indices
! R(3)    ... neutron coordinate
! K(3)    ... neutron k-vector
! *** RETURN:
! IC    .. which wall:  left(1), right(2), top(3), bottom (4)
!          if no contact, pass through, IC=0
! DT    .. time to reach the next contact point
! Q     .. reflection vector (magnitude)
! N     .. unit vector normal to the surface
!----------------------------------------------------------------------
      USE GUIDES
      USE TRACINGDATA
      USE CONSTANTS
      use MIRROR_TABLE
      use SIMATH
      IMPLICIT NONE
      TYPE(GUIDE) :: OBJ
      real(kind(1.D0)), intent(in) :: K(3)
      real(kind(1.D0)), intent(out) :: DT,Q,NW(3),KG(3),P
      real(kind(1.D0)), intent(inout) :: R(3)
      INTEGER,intent(in) :: IH,IV
      INTEGER,intent(inout) :: IC
      real(kind(1.D0)) :: GET_CROSSTIME,GUIDE_LAM
      real(kind(1.D0)) :: LZ(4),ANG,T0,LPOS,N(3)
      INTEGER :: I,ISG,IC0
      INTEGER :: IX,ILAM(4)
      logical :: TRANS
      integer,parameter :: IL(4)=(/1,0,1,0/)
      integer,parameter :: ID(4)=(/1,1,2,2/)
      integer,parameter :: IDINV(4)=(/2,2,1,1/)

! shifts for an oscillating collimator
      real(kind(1.D0)) :: OSCD,OSCA
      COMMON /OSCBEND/ OSCD,OSCA
1     format('GUIDE_TRANS ',a,2x,6(G10.4,2x))

!     remember input IC
      IC0=IC
! TOF to the guide end
      if (TRACING_UP) then
        T0=-R(3)/K(3)
      else
        T0=(OBJ%FRAME%SIZE(3)-R(3))/K(3)
      endif
! TOF to the next boundary in units h_bar/m=1
      ILAM=(/IH,IH,IV,IV/)
      DO I=1,4
        TRANS=(ID(IC).EQ.ID(I))
        LZ(I)=GET_CROSSTIME(OBJ,R,K,ILAM(I),I,TRANS)
      ENDDO

! select the shortest time to cross a surface
      DT=MIN(LZ(1),LZ(2),LZ(3),LZ(4),T0)
! DT<0 should never happen, but to be sure, ...
      if (GUIDES_dbg)  write(*,1) '   DT,T0',DT,T0
      IF (DT.LE.0D0) goto 99 ! signal from GET_CROSSTIME to stop

 ! and find which of the 4 surfaces it was
      IC=0
      DO I=1,4
        IF (DT.EQ.LZ(I)) IC=I
      ENDDO
      if (GUIDES_dbg)  write(*,1) '   IC0,IC',IC0,IC

! if no contact => pass through
      IF (IC.EQ.0) THEN
        DT=T0
        return
      ENDIF
      IX=ID(IC)

! stop events trying to migrate between horiz. and vertical directions
      if (ID(IC).NE.ID(IC0)) then
        NEUT%P=0.D0
        goto 99
      endif

! get the surface normal vector at the contact point
      N(1:3)=0.D0
      ISG=1-2*IL(IC) ! /1,0,1,0/ ==> /-1,1,-1,1/
      ANG=GUIDE_LAM(OBJ,1,ILAM(IC),R(3)+DT*K(3),IX)
      LPOS=GUIDE_LAM(OBJ,0,ILAM(IC),R(3)+DT*K(3),IX)
      if (IX==1) then
         LPOS=LPOS + OSCD + (R(3)+DT*K(3))*OSCA
         ANG=ANG+OSCA
      endif

      N(3)=-ISG*ANG/SQRT(1.D0+ANG**2)
      N(IX)=ISG*SQRT(1.D0-N(3)**2)
      ! add waviness

! stop events reflected from concave side
      if ((OBJ%ONESIDE.GT.0).and.(DT.LT.T0)) then
        select case (IX)
        case(1)
          IF (N(1)*OBJ%ROH.GT.0) goto 99
        case(2)
          IF (N(2)*OBJ%ROV.GT.0) goto 99
        end select
      endif

      ! write(*,1) trim(OBJ%FRAME%ID),IC,K
      call MIRROR_REFLECT(OBJ%WAV,K,N,IC,KG,NW,Q,P)
      if (GUIDES_dbg)  write(*,1) '    Q,P,STATE',Q,P,NEUT%STATE

      if (P<=0.D0) goto 99
      return


! signal to stop the event
99    NEUT%STATE=1
      IC=0
      END SUBROUTINE GUIDE_TRANS

!---------------------------------------------------------------------
      REAL*8 FUNCTION GET_CROSSTIME(OBJ,R,K,ICHAN,ISIDE,TRANS)
! Get time-of-flight to the cross-section with the closest wall of a channel
! If OBJ%MU<1D5, handles also passage inside associated lamella on the right/bottom
! The time must be positive and is expressed in units h_bar/m=1
! 1) no solution .. return 10^30
! 2) 1 solution  .. return this
! 3) 2 solutions .. return the smaller one
! *** INPUT:
! R(3)    ... neutron coordinate
! K(3)    ... neutron k-vector
! ICHAN   ... index of the channel: (0..OBJ%NL-1)
! ISIDE   ... side of the wall in the current channel:
!             LEFT (1) RIGHT (2) TOP (3) BOTTOM (4)
! TRANS   ... true for tracing inside a lamella
! *** RETURN:
! Z-distance (not time !) to reach the next contact point
! time = result/K(3)
!----------------------------------------------------------------------
      USE GUIDES
      USE CONSTANTS
      use TRACINGDATA
      use TRACELIB
      IMPLICIT NONE
      REAL*8,parameter :: EPS=1.D-8
      TYPE(GUIDE) :: OBJ
      REAL*8,intent(in) :: R(3),K(3)
      INTEGER,intent(in) :: ICHAN,ISIDE
      logical, intent(in) :: TRANS
      integer :: ISG,INEXT,IX,NLAM(2)
      logical :: IN_LAM,SIDE_WALL
      REAL*8 :: RHO(2),DLAM(2),DL,LPOS,LANG,WPOS,WANG,GUIDE_CROSS1,GUIDE_LAM
      REAL*8 :: XTIME,GX
! depending on index ISIDE, the following flags determine:
  ! the side of a lamella, IL: left and top (1) , right and bottom (0)
      integer,parameter :: IL(4)=(/1,0,1,0/)
  ! direction, ID: horizontal(1), vertical(2)
      integer,parameter :: ID(4)=(/1,1,2,2/)

! shifts for an oscillating collimator
      REAL*8 :: OSCD,OSCA
      COMMON /OSCBEND/ OSCD,OSCA

      WPOS=0.D0
      LPOS=0.D0
      LANG=0.D0


2     format(a12,' ',8(G12.6,1x))

! to simplify code, store horizontal and vertical guide parameters in arays
      NLAM=(/OBJ%NLH,OBJ%NLV/)
      DLAM=(/OBJ%DLH,OBJ%DLV/)
      RHO =(/OBJ%CH,OBJ%CV/)
      IX=ID(ISIDE) ! normal coordinate (1 or 2)

! for outer wall, do not include lamella thickness !

! INEXT in a channel: = 1 for LEFT/TOP wall, else = 0
! ISG in a channel: -1 for LEFT,TOP wall, else +1
! ISG in a lamella: oposite to channel
      INEXT=IL(ISIDE)
      ISG=1-2*INEXT ! /1,0,1,0/ ==> /-1,1,-1,1/
      IF (TRANS) THEN
        INEXT=0
        ISG=-ISG
      ENDIF

!       if (dbg) write(*,2) 'IX,IW,RX:',IX,ICHAN+INEXT,R(IX)

! check whether the wall is an outer wall of the collimator
      SIDE_WALL=(ICHAN+INEXT.EQ.0).OR.(ICHAN+INEXT.EQ.NLAM(IX))
      ! there are no side walls if OBJ%LAMONLY=true
      if (OBJ%LAMONLY.and.SIDE_WALL) then
        GET_CROSSTIME=TLIB_INF
        return
      endif

! Get lamella parameters:
! RHO is curvatures
! DL   .. the half-width of the lamella
! LANG .. lamella angle at z=R(3)
! LPOS .. lamella middle position at z=R(3)
! WPOS .. wall position at z=0
! WANG .. wall angle at z=0

      DL=DLAM(IX)/2.D0
      IF (SIDE_WALL) DL=0.D0
      select case (IX)
      case(1)
        WPOS=OBJ%LH(ICHAN+INEXT) + ISG*DL + OSCD + R(3)*OSCA
        WANG=OBJ%AH(ICHAN+INEXT) + OSCA
        LPOS=GUIDE_LAM(OBJ,0,ICHAN+INEXT,R(3),IX) + OSCD + R(3)*OSCA
        LANG=GUIDE_LAM(OBJ,1,ICHAN+INEXT,R(3),IX) + OSCA
      case(2)
        WPOS=OBJ%LV(ICHAN+INEXT) + ISG*DL
        WANG=OBJ%AV(ICHAN+INEXT)
        LPOS=GUIDE_LAM(OBJ,0,ICHAN+INEXT,R(3),IX)
        LANG=GUIDE_LAM(OBJ,1,ICHAN+INEXT,R(3),IX)
      end select

! do not allow tangential neutrons to pass
      if ((ABS(RHO(IX)).LT.EPS).and.(ABS(LANG-K(IX)/K(3)).LT.EPS)) then
!         write(*,*) 'rejected tangential event in ',trim(OBJ%FRAME%NAME),OBJ%FRAME%COUNT
        GET_CROSSTIME=TLIB_INF
        Return
      endif

! calculate TOF to reach the wall
      !XTIME=GUIDE_CROSS(K(IX),K(3),R(IX)-WPOS,R(3),WANG,RHO(IX))

      gx=NEUT_GR(IX)
      XTIME=GUIDE_CROSS1(K(IX),K(3),R(IX)-WPOS,R(3),gx,WANG,RHO(IX))


      GET_CROSSTIME=XTIME
      ! if (dbg) write(*,2) 'IX,IW,RX,WX,T:',IX,ICHAN+INEXT,R(IX),WPOS,XTIME

! --------------- DEBUG option -----------------------
! for debuging, check consistency between TRANS and IN_LAM:
! set IN_LAM=true if we are inside an inner lamella
      IN_LAM=(ABS(R(IX)-LPOS).LT.DL)
      if (IN_LAM.neqv.TRANS) then
        GET_CROSSTIME=0.D0
!         write(*,*) 'rejected edge event in ',trim(OBJ%FRAME%NAME),OBJ%FRAME%COUNT

!         MSG='channel'
!         if (TRANS) MSG='lamella'
!         write(*,*) 'Fatal error in '//trim(OBJ%FRAME%NAME)//': '
!         write(*,*) '   Unexpected position, claimed to be inside ',trim(MSG),' ',ICHAN
!         write(*,*) '   X(neutron)=',R(IX),' X(lamella)=',LPOS,' X(wall)=',WPOS
!         write(*,*) '   I(wall)=',ICHAN+INEXT,' IX=',IX, 'x-time=',XTIME
!         write(*,*) '   distance=',R(3) ,'  thickness=+-',DL,' count=',OBJ%FRAME%count
!         write(*,*) '   program exits ...'
!         read(*,*)
!         STOP
      endif
! --------------- DEBUG option -----------------------

      END FUNCTION GET_CROSSTIME

!       end MODULE
