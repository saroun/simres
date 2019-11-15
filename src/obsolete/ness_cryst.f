	
C	---------------------------------------
	SUBROUTINE CRYST_INIT2(CR)
C	---------------------------------------
      implicit none
      INCLUDE 'structures.inc'	                
      RECORD /CRYSTAL/ CR
      REAL*8 Z,sec
      PARAMETER (sec=4.85D-6)
      REAL*8 GETQKIN,GETREFDYN,GETMI
      INTEGER*4 I,J

C///  THB,DHKL,etc.. must be specified before !

        CALL SLIT_INIT(CR.FRAME)
        
        CR.GTOT=2*PI/ABS(CR.DHKL)
        CR.G(1)=CR.GTOT*SIN(CR.CHI)
        CR.G(2)=0.
        CR.G(3)=CR.GTOT*COS(CR.CHI)
        CR.STMCH=SIN(CR.THB-CR.CHI)
        CR.CTMCH=COS(CR.THB-CR.CHI)
	CR.LAMBDA=2.*CR.DHKL*SIN(CR.THB)
        DO 10 I=1,3
          CR.MAPG(I)=.FALSE.
          DO 10 J=1,3 
            CR.DG_DR(I,J)=0
10      CONTINUE               

c G-gradient for elastically bent crystal      
        IF(CR.HMOS.LT.sec.AND.CR.NH.LE.1) THEN 
          CR.DG_DR(1,1)=-COS(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(1,3)=SIN(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(3,1)=SIN(CR.CHI)*CR.GTOT*CR.RH
          CR.DG_DR(2,2)=0    ! no vertical bending
          CR.DG_DR(3,3)=-CR.POI*COS(CR.CHI)*CR.GTOT*CR.RH
	  CR.MAPG(1)=.TRUE.        
	  CR.MAPG(3)=.TRUE.        
	ENDIF
c d-gradient
        IF(CR.DGR.NE.0.) THEN 	
	  Z=1.D-4*CR.GTOT*CR.DGR
          CR.DG_DR(1,1)=CR.DG_DR(1,1)+Z*cos(CR.DGA+CR.CHI)
          CR.DG_DR(1,3)=CR.DG_DR(1,3)-Z*sin(CR.DGA+CR.CHI)
          CR.DG_DR(3,1)=CR.DG_DR(3,1)-Z*sin(CR.DGA+CR.CHI)
          CR.DG_DR(3,3)=CR.DG_DR(3,3)-Z*cos(CR.DGA+CR.CHI)
	  CR.MAPG(1)=.TRUE.        
	  CR.MAPG(3)=.TRUE.        
	ENDIF  		
c unit vector |- to G
	CR.gama(1)=COS(CR.CHI)
        CR.gama(3)=-SIN(CR.CHI)
        CR.gama(2)=0 

        CR.QHKL=GETQKIN(CR,CR.lambda)  ! kin. reflectivity
        CR.dext=CR.DHKL/CR.LAMBDA*SQRT(4*PI/CR.QML) ! extinction length
        Z=CR.dlam/CR.Dext
        IF(Z.GT.1D-5) THEN 
          CR.Ext1=TANH(Z)/Z         ! primary extinction
        ELSE
          CR.Ext1=1.
        ENDIF	 
        CR.MI=GETMI(CR,CR.lambda,300.D0)  ! absorption coefficient
        CR.REF=GETREFDYN(CR,CR.LAMBDA)    ! peak reflectivity (no mosaic)         
        CR.DELTA=CR.QHKL*CR.Dext*1D-4/PI   ! Darwin box width
	IF(CR.HMOS.GT.SEC) CR.HMOS=MAX(CR.HMOS,CR.DELTA)

        RETURN        
        END

C	-------------------------------------------
      LOGICAL FUNCTION CRYST_GO2(CRYST,NEUI,NEUF)
C	-------------------------------------------
      implicit none
      REAL*8 sec
      PARAMETER(sec=4.85D-6)

      INCLUDE 'ness_common.inc'
	                
      RECORD /CRYSTAL/ CRYST
      RECORD /NEUTRON/ NEUI,NEUF
      REAL*8 V(3),K(3),R(3),KaG(3)
      REAL*8 PP,DT,H
      real*8 CRYST_MOS,CRYST_BENT
      INTEGER*4 i
      LOGICAL*4 REP

101   format('CR: ',7(G11.5,2x),a1)
      REP=(CRYST.FRAME.NAME(1:2).EQ.'Si')
      write(*,*) 'CR: ',CRYST.FRAME.NAME

C/// if nh=0 then accept neutrons without transformations
      IF (CRYST.NH.EQ.0) then           
	neuf=neui
	dt=(cryst.frame.dist-neui.r(3))/neui.k(3)
	do i=1,3
	      neuf.r(i)=neui.r(i)+dt*neui.k(i)
	end do
        NEUF.T=NEUI.T+DT/HOVM
	CRYST_GO2=.TRUE.
	return
      endif   

      if (CRYST.MAG*NEUI.S.LT.0) GOTO 99      

        
        IF (REP) WRITE(*,101) (NEUI.R(I),I=1,3),(NEUI.K(I),I=1,3)
      CALL SLIT_PRE(CRYST.FRAME,NEUI.R,NEUI.K,V,K)
        IF(REP) WRITE(*,101) (V(I),I=1,3),(K(I),I=1,3)
      V(1)=V(1)+CMX              
C test right barrier
      IF (CBAR.GT.0) THEN
         H=V(3)-(CRYST.FRAME.SIZE(1)/2.+V(1)+40.)/K(1)*K(3)
         IF (ABS(H).LT.CBAR.AND.K(1).GT.0) THEN
           CRYST_GO2=.FALSE.
           NEUF.P=0.
           RETURN
         ENDIF          
      ENDIF
c choose appropriate routine for reflection simulation
      IF (CRYST.HMOS.GE.SEC.OR.
     *   (CRYST.NH.GT.1.AND.CRYST.DGR.NE.0.)) THEN
         PP=CRYST_MOS(CRYST,V,K,R,KaG,DT)
      ELSE
         PP=CRYST_BENT(CRYST,V,K,R,KaG,DT)
      ENDIF
      IF (PP.GT.1.D-4) THEN
C test right barrier
          IF (CBAR.GT.0) THEN
            H=R(3)-(CRYST.FRAME.SIZE(1)/2.+R(1)+40.)/KaG(1)*KaG(3)
            IF (ABS(H).LT.CBAR.AND.KaG(1).LT.0) THEN
              CRYST_GO2=.FALSE.
              NEUF.P=0.
              RETURN
            ENDIF          
          ENDIF
C transform to local axis coordinate and return
        IF(REP) WRITE(*,101) (R(I),I=1,3),(KaG(I),I=1,3)
          CALL SLIT_POST(CRYST.FRAME,R,KaG,NEUF.R,NEUF.K)
        IF (REP) WRITE(*,101) (NEUF.R(I),I=1,3),(NEUF.K(I),I=1,3)
          NEUF.R(1)=NEUF.R(1)-CMX              	           
	  NEUF.S=NEUI.S
          NEUF.T=NEUF.T+DT/HOVM
          NEUF.P=NEUI.P*PP
          CRYST.FRAME.COUNT=CRYST.FRAME.COUNT+1
          CRYST_GO2=.TRUE.
      ELSE
          CRYST_GO2=.FALSE.
      ENDIF
        
      RETURN
99    NEUF.P=0
      CRYST_GO2=.FALSE.
      RETURN        
      END        

C----------------------------------------------------------------
      REAL*8 FUNCTION CRYST_BENT(CR,R1,K1,R2,K2,DT)
C Transmission function for elastically bent crystal 
C (eta<")      
C----------------------------------------------------------------
      implicit none
      INCLUDE 'structures.inc'
      RECORD /CRYSTAL/ CR
      REAL*8 Z,C,B1,KK,T0,TIN,TOUT,gabs 
      REAL*8 KaG(3),G(3),DGK(3),G0(3),alfv ,dv
      REAL*8 R1(3),R2(3),K1(3),K2(3),DT,A(3),K0,REF,sigma,costh
      INTEGER*4 i,j,iv
        
      call CR_BORDER(CR,R1,K1,TIN,TOUT)
      TIN=TIN+0.00001*(TOUT-TIN)   ! get inside the crystal, avoid rounding error
      DT=TIN
      IF (TIN.GE.TOUT) GOTO 99     ! No intersection with the crystal
      DO I=1,3 
           R2(I)=R1(I)+TIN*K1(I)
      END DO

      iv=nint((0.5+R2(2)/CR.FRAME.SIZE(2))*CR.nv+0.5)
                
      IF((iv.LT.1).OR.(iv.GT.CR.nv)) then 
         write(*,*) 'Unexpected error in CR_BENT: ',iv,R2(2)
	 goto 99
      ENDIF 
      
      dv=(2*iv-CR.nv-1)*CR.FRAME.SIZE(2)/2./CR.nv
C Calculate local G-vector 
      DO I=1,3
        G0(I)=CR.G(I) 
        IF (CR.MAPG(I)) THEN
          A(I)=0.
	  DO J=1,3
	    A(I)=A(I)+CR.DG_DR(I,J)*R2(J)
	  ENDDO
	  G0(I)=G0(I)+A(I)
	ENDIF  
      END DO
      GABS=SQRT(G0(1)**2+G0(2)**2+G0(3)**2)
C Add segment tilt angle and vertical mosaic spread
      alfv=dv*CR.RV
      G(1)=G0(1)   
      G(3)=G0(3) 
      IF (ABS(G0(3)).LT.0.01*GABS) THEN  ! symmetric Laue case          
         G(2)=G0(2)-G0(1)*alfv
      ELSE
         G(2)=G0(2)-G0(3)*alfv  ! other case       
      ENDIF
      Z=GABS/SQRT(G(1)**2+G(2)**2+G(3)**2)
      DO i=1,3
        G(i)=G(i)*Z
      END DO
      DO I=1,3
	KaG(I)=K1(I)+G(I)
      END DO

      KK=K1(1)**2+K1(2)**2+K1(3)**2  
      K0=SQRT(KK)
      C=KaG(1)**2+KaG(2)**2+KaG(3)**2-KK       
      B1=0.
      DO I=1,3
        DGK(I)=0.
	IF(CR.MAPG(I)) THEN
	  DO J=1,3
	    DGK(I)=DGK(I)+CR.DG_DR(I,J)*K1(J)
	  END DO 
	  B1=B1+KaG(I)*DGK(I)
	ENDIF
      END DO
      B1=2.*B1  

c1     format(a10,6(1x,G12.6))	
c      write(*,1) 'Cryst: ',C,B1,-C/B1*K0,(TOUT-TIN)*K0  
c      write(*,1) 'R1,K1: ',(R2(I),I=1,3),(K1(I),I=1,3)       
        
      IF(ABS(B1).LT.ABS(C)*1.0D-10) goto 99   ! avoid floating overflow
      T0=-C/B1	   
      IF (T0.LE.0.D+0.OR.T0.GE.(TOUT-TIN)) GOTO 99
c      write(*,1) 'G,DG: ',(G(I),I=1,3),(T0*DGK(I),I=1,3)
      DO I=1,3
          R2(I)=R2(I)+K1(I)*T0
          G(I)=G(I)+DGK(I)*T0
      END DO 
       
      DO I=1,3
          K2(I)=K1(I)+G(I)
      END DO    
      Z=sqrt(KK/(K2(1)**2+K2(2)**2+K2(3)**2))
      do i=1,3
         K2(i)=K2(i)*Z
      end do  

c      write(*,1) 'R2,K2: ',(R2(I),I=1,3),(K2(I),I=1,3)

      call CR_BORDER(CR,R2,K2,TIN,TOUT) 
      DO I=1,3
           R2(I)=R2(I)+K2(I)*TOUT
      end do   
      T0=T0+TOUT
      DT=DT+T0
      COSTH=K2(1)*CR.GAMA(1)+K2(2)*CR.GAMA(2)+K2(3)*CR.GAMA(3)
      B1=ABS(B1/2./COSTH)
      sigma=K0*CR.QHKL*CR.DW*CR.EXT1/10./B1*GABS      
      REF=1-EXP(-sigma)
      CRYST_BENT=REF*EXP(-K0*T0/10.*CR.MI)
      
c      write(*,1) 'R3,DR: ',(R2(I),I=1,3),(TOUT*K2(I),I=1,3)
c      write(*,1) 'Ref: ',REF*EXP(-K0*T0/10.*CR.MI)
c      pause


c      if (CR.FRAME.COUNT.EQ.1) then
c         write(*,*) REF,EXP(-K0*T0/10.*CR.MI),K0*T0/10.
c      endif
      RETURN
99    CRYST_BENT=0.D+0  
      END


C-----------------------------------------------------------
      REAL*8 FUNCTION CRYST_MOS(CR,R1,K1,R2,K2,DT)
C Transmission function for mosaic crystal 
C (eta>")  
C simplified version (Darwin width=0)    
C-----------------------------------------------------------
      implicit none
      INCLUDE 'structures.inc'
c      INCLUDE 'ness_const.inc'
c      INCLUDE 'crystal.inc'
      REAL*8 EPS
      PARAMETER (EPS=1D-5)
      RECORD /CRYSTAL/ CR
      REAL*8 R1(3),K1(3),R2(3),K2(3),DT
      REAL*8 G(3)
      INTEGER*4 I
      REAL*8 PP,TIN,TOUT,T0,tau
      REAL*8 KK,K0,Qhkl,Z,PATH,MI

      PP=1.D0
C/// calculate values, which are constant in the course of the iteration      
      KK=K1(1)**2+K1(2)**2+K1(3)**2
      K0=SQRT(KK)      
      Qhkl=CR.Qhkl/10.*CR.DW*CR.Ext1
      MI=CR.MI/10.
      
      call CR_BORDER(CR,R1,K1,TIN,TOUT)
      IF (TIN.GE.TOUT) GOTO 99     ! No intersection with the crystal
      T0=TIN
      DT=T0
      TOUT=TOUT-TIN
C/// move to the crystal entry point         
      DO I=1,3
         R2(I)=R1(I)+(TIN+0.00000001*TOUT)*K1(I)
      END DO

C ****  Begin of multiple scattering iteration cycle  ****

C// generate random walk step in K direction
C---------------------------------------------------
50    CALL WALKSTEP(CR,1,R2,K1,K0,Qhkl,TAU,TOUT,G,PP)      
      if (PP.LT.EPS) goto 99   ! don't follow trajectories with low probability         
C// move to the new reflection point and set K=K+G
      DT=DT+tau
      DO I=1,3
	 K2(I)=K1(I)+G(I)
      END DO
      Z=SQRT(KK/(K2(1)**2+K2(2)**2+K2(3)**2))  ! set |K2| = |K1|
      DO I=1,3
	   K2(I)=K2(I)*Z
      END DO
c      IF(CR.FRAME.COUNT.LT.5) write(*,1) 'K0: ',K0*TAU,K0*TOUT,K0*DT

C// generate 2nd random walk step in K+G direction
C-------------------------------------------------
      CALL WALKSTEP(CR,-1,R2,K2,K0,Qhkl,TAU,TOUT,G,PP)      
c      IF(CR.FRAME.COUNT.LT.5) 
c     *  write(*,1) 'KG: ',K0*TAU,K0*TOUT,K0*(DT+TAU)
      IF (PP.LT.EPS) GOTO 99
      if (tau.GE.TOUT) THEN    
	DT=DT+TOUT
	PATH=(DT-T0)*K0
	PP=PP*EXP(-MI*PATH)
      else       
C// move to the new reflection point and set K=K+G
        DT=DT+tau
        DO I=1,3
           K1(I)=K2(I)+G(I)
        END DO
        Z=SQRT(KK/(K1(1)**2+K1(2)**2+K1(3)**2))  ! set |K2| = |K1|
        DO I=1,3
	   K1(I)=K1(I)*Z
        END DO	

c      write(*,*) Z
c        pause 'Error'
 	GOTO 50  ! repeat cycle
      endif
	
C ****  End of multiple scattering cycle  ****
      
      CRYST_MOS=PP
c      IF(CR.FRAME.COUNT.LT.5) THEN
c1       format(a,5(G15.8,1x))
c        write(*,1) 'R2: ',(R2(I),I=1,3)
c        write(*,1) 'K2: ',(K2(I),I=1,3)
c        write(*,1) 'DT,P,P0: ',PATH,PP/EXP(-MI*PATH),PP
c        IF(CR.FRAME.COUNT.EQ.4) THEN
c          pause
c	endif
c      ENDIF
      RETURN
        
99    CRYST_MOS=0.D0
      END	 




C-----------------------------------------------------------------
      REAL*8 FUNCTION CRYST_MOS_SIMPLE(CR,R1,K1,R2,K2,DT)
C Transmission function for mosaic crystal 
C (eta>") Simplified version according to Bacon-Lowde model     
C-----------------------------------------------------------------
      implicit none
      INCLUDE 'structures.inc'
c      INCLUDE 'crystal.inc'
c      INCLUDE 'ness_MC.inc'
      REAL*8 R1(3),K1(3),R2(3),K2(3),DT
      REAL*8 PP,GABS,Z,C,B2,TIN,TOUT,KK
      REAL*8 ETA1,G(3),KaG(3)
      REAL*8 CRYSTAL_RC
      INTEGER*4 I
      RECORD /CRYSTAL/ CR
      REAL*4 RAN1
      REAL*8 V3XV3,CRYST_TILT
      
      PP=1
C/// move to the crystal interior         
      call CR_BORDER(CR,R1,K1,TIN,TOUT)
      IF (TIN.GE.TOUT) then     ! No intersection with the crystal
	 PP=0
	 GOTO 99
      ENDIF	
      DT=(TIN+TOUT)/2+(RAN1()-0.5)*(TOUT-TIN)
      DO I=1,3
         R2(I)=R1(I)+DT*K1(I)
      END DO
      

      GABS=CRYST_TILT(CR,R2,G)

      CALL V3AV3(1,K1,G,KaG)            ! K+G
      KK=V3XV3(K1,K1)
      C=V3XV3(KaG,KaG)-KK     ! (K+G)^2-K^2        
      B2=2.*V3xV3(KaG,CR.gama)*GABS
      
      ETA1=-C/B2
      IF(ABS(ETA1).GT.CR.DETA) THEN
        PP=0
	goto 99
      ENDIF
      DO  I=1,3
         G(I)=G(I)+GABS*CR.GAMA(I)*ETA1
      END DO
      Z=SQRT(V3xV3(G,G))
      DO  I=1,3
           G(I)=G(I)*GABS/Z
           K2(I)=K1(I)+G(I)
      END DO
      Z=sqrt(KK/v3xv3(K2,K2))
      DO i=1,3
         K2(i)=K2(i)*Z
      ENDDO   
      PP=CRYSTAL_RC(CR,ETA1,0.D0)
99    CRYST_MOS_SIMPLE=PP
      END
      

CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCx
C
C  ****  Auxilliary  **** 
C
CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCx

C---------------------------------------------------------------
      REAL*8 FUNCTION HMOS_DIST(X)
C Distribution of horizontal angular deviation of mosaic blocks
C --------------------------------------------------------------
      IMPLICIT NONE
      include 'ness_common.inc'
      REAL*8 X,SQRT2PI,SQRT2LN2,DD0
      PARAMETER (SQRT2PI=2.506628274631, SQRT2LN2=1.177410023, 
     *           DD0=0.270347526)      

C Gauss
      IF (MDIST.EQ.0) THEN	
	  HMOS_DIST=EXP(-0.5*X**2)/SQRT2PI
      ELSE  
C pseudo-Voigt
        HMOS_DIST=0.5*DD0/(1. + (X/SQRT2LN2)**2)+
     *   0.5*EXP(-0.5*X**2)/SQRT2PI
      ENDIF
      END 

C ---------------------------------------------------------------------C
      SUBROUTINE CR_BORDERS(CR,DIR,R,K,K0,Q,PHI)
C Traces through all segments along the neutron path
C PHI= vertical mosaic tilt. angle 
C Q = Qhkl*DW*Ext1
C K0=|K|
C R,K position and K-vector
C DIR=1,-1 for directions K,K+G
C IN /CRBORDERS/ returns for each segment crossed:
C TSEG    time of intersection with segment rear border
C ALPHA   dThetaB, deviation from Bragg angle 
C GRAD    grad(ThetaB), gradient of Bragg angle along K 
C PSEG    scattering probability
C NSEG    number of crossed segments
C ---------------------------------------------------------------------C
      implicit none
      INCLUDE 'structures.inc'
      RECORD /CRYSTAL/ CR
      REAL*8 EPS,SEC
      PARAMETER (EPS=1D-8,SEC=4.85D-6)
      INTEGER*4 MSEG
      PARAMETER (MSEG=50)
      REAL*8 TSEG(0:MSEG),ALPHA(0:MSEG),GRAD(0:MSEG),PSEG(0:MSEG),
     *       PSEG0(0:MSEG)
      INTEGER*4 NSEG
      INTEGER*4 I,IH,IV,DIR,J,I1
      REAL*8 R(3),K(3),Q,PHI,DX,X0,DUM,T1(3),T2(3),V(3),A,B,sigma,DT
      REAL*8 KK,K0,G(3),KaG(3),alfh,alfv,GABS,Z,G0(3),V0(3),W(3)
      REAL*8 HMOS_DIST,ERF,Z1,Z2
      COMMON /CRBORDERS/ TSEG,PSEG,PSEG0,ALPHA,GRAD,NSEG 
c      REAL*4 RAN1
      DO I=1,3
         V(I)=R(I)
      END DO
      J=0
      TSEG(0)=0
      PSEG(0)=0.
      KK=K0**2
      DO WHILE ((ABS(V(1)/CR.FRAME.SIZE(1)).LT.0.5).AND.
     *          (ABS(V(2)/CR.FRAME.SIZE(2)).LT.0.5).AND.
     *          (ABS(V(3)/CR.FRAME.SIZE(3)).LT.0.5).AND.J.LT.MSEG)
      
C Find time of intersection with exit surface 
	J=J+1
        ih=nint((0.5+V(1)/CR.FRAME.SIZE(1))*CR.nh+0.5)
        iv=nint((0.5+V(2)/CR.FRAME.SIZE(2))*CR.nv+0.5)
        DO I=1,3
         IF (ABS(K(I)).GT.EPS) THEN
           if (I.EQ.1) then
             x0=ih-(CR.nh+1)/2.0D0
             dx=CR.FRAME.SIZE(I)/CR.nh
           else if (I.EQ.2) then
             x0=iv-(CR.nv+1)/2.0D0
             dx=CR.FRAME.SIZE(I)/CR.nv
           else if (I.EQ.3) then
             x0=0.D0
             dx=CR.FRAME.SIZE(I)
           endif
       	   T2(I)=((x0+0.5D0)*dx - V(I))/K(I)
       	   T1(I)=((x0-0.5D0)*dx - V(I))/K(I)
           IF (T1(I).GT.T2(I)) THEN
	     DUM=T1(I)
             T1(I)=T2(I)
	     T2(I)=DUM
           ENDIF 
         ELSE
	   T2(I)=1.0D30
	   T1(I)=-1.0D30
         ENDIF
        END DO
        TSEG(J)=TSEG(J-1)+MIN(T2(1),T2(2),T2(3))
	DT=TSEG(J)-TSEG(J-1)
C// Coordinates of the segment centre
        V0(1)=(2*ih-CR.nh-1)*CR.FRAME.SIZE(1)/2./CR.nh
        V0(2)=(2*iv-CR.nv-1)*CR.FRAME.SIZE(2)/2./CR.nv
	V0(3)=0.
	
C Calculate local G-vector 
        DO I=1,3
          G0(I)=CR.G(I) 
          IF (CR.MAPG(I)) THEN
            W(I)=0.
	    DO I1=1,3
	      W(I)=W(I)+CR.DG_DR(I,I1)*(V(I1)-V0(I1))
	    ENDDO
	    G0(I)=G0(I)+W(I)
	  ENDIF  
        END DO
        GABS=SQRT(G0(1)**2+G0(2)**2+G0(3)**2)
C Add segment tilt angle and vertical mosaic spread
        alfh=V0(1)*CR.RH
        alfv=V0(2)*CR.RV
        G(1)=G0(1)-G0(3)*alfh   
        G(3)=G0(3)+G0(1)*alfh  
        IF (ABS(G0(3)).LT.0.01*GABS) THEN              
           G(2)=G0(2)-G0(1)*alfv+GABS*PHI ! only in Laue case
        ELSE
           G(2)=G0(2)-G0(3)*alfv+GABS*PHI
        ENDIF
        Z=GABS/SQRT(G(1)**2+G(2)**2+G(3)**2)
        DO i=1,3
          G(i)=DIR*G(i)*Z
        ENDDO
C get delta(ThetaB) (=alpha)
        DO I=1,3
	   KaG(I)=K(I)+G(I)
	END DO
        a=KaG(1)**2+ KaG(2)**2+KaG(3)**2-KK    ! (K+G)^2-K^2 
        b=DIR*GABS*(KaG(1)*CR.GAMA(1)+KaG(2)*CR.GAMA(2)+
     *         KaG(3)*CR.GAMA(3))   
        alpha(J)=-a/(2.*b)
C// add random angle within Darwin box width
c	alpha(J)=alpha(J)+(RAN1(1)-0.5)*Q*CR.DEXT*1D-3/PI
C get grad(ThetaB) (=grad)
        grad(J)=0.
        DO I=1,3
          IF (CR.MAPG(I)) THEN
            Z=0.D0
	    DO I1=1,3
	      Z=Z+DIR*CR.DG_DR(I,I1)*K(I1)
	    ENDDO
	    grad(J)=grad(J)+KaG(I)*Z
	  ENDIF  
        END DO
	grad(J)=-grad(J)/b
C get scattering probability
        IF(abs(grad(J)*DT).LE.1.D-6*ABS(alpha(J))) THEN   ! mosaic only
	  sigma=K0*Q*HMOS_DIST(alpha(J)/CR.HMOS)/CR.HMOS*DT
	ELSE  ! mosaic with gradient
	  IF(CR.HMOS.GT.SEC) THEN
	    Z1=ERF((alpha(J)+grad(J)*DT)/CR.HMOS,0)
	    Z2=ERF(alpha(J)/CR.HMOS,0)
	    sigma=K0*Q*(Z1-Z2)/grad(J)
	    IF(DIR.LT.0) THEN
	      IF(J.EQ.1.AND.ABS(ALPHA(J)).LT.SEC) sigma=0.D0 ! no back refl. in the same segment
	    ENDIF
	  ELSE   ! no mosaic
	    Z1=-alpha(J)/grad(J)
	    IF(Z1.GT.0.AND.Z1.LT.DT) THEN  
	       sigma=K0*Q/abs(grad(J))
	       IF (DIR.LT.0.AND.J.EQ.1) sigma=0.D0 ! no back refl. in the same segment
	    ELSE
	       sigma=0.D0
	    ENDIF
	  ENDIF      	    	    
        ENDIF	
	IF(sigma.gt.14) then
	  PSEG0(J)=1.D0
	  PSEG(J)=1.D0
	else  
  	  PSEG0(J)=1.D0-exp(-sigma)
	  PSEG(J)=1.D0-(1.D0-PSEG(J-1))*(1.D0-PSEG0(J))
	endif  
	
C// Move to the next segment        
	DO I=1,3
          V(I)=V(I)+DT*1.00000001D0*K(I)
        ENDDO       
c      write(*,*) 'V: ',(V(I),I=1,3)
      ENDDO
      NSEG=J
      
c1     format(a10,5(G10.4,1x))
c      write(*,*) 'R: ',(R(I),I=1,3)
c      write(*,*) 'K: ',(K(I),I=1,3)
c      write(*,*) 'DIR: ',DIR,CR.FRAME.COUNT
c      write(*,1) 'PSEG : ',(PSEG(I),I=1,NSEG)       
c      write(*,1) 'PSEG0: ',(PSEG0(I),I=1,NSEG)       
c      write(*,1) 'ALPHA: ',(ALPHA(I),I=1,NSEG)       
c      write(*,1) 'TSEG: ',(TSEG(I),I=1,NSEG)       
c      pause
      END

C--------------------------------------------------------------
      SUBROUTINE WALKSTEP(CR,DIR,R,K,K0,Q,TAU,TOUT,G,PP)
C Generate random walk step in the crystal
C//   Q = Qhkl*DW*Ext1
C//   K0=|K|
C//   R,K position and K-vector
C//   DIR=1,-1 for directions K,K+G
C//   TAU.. Time to reach reflection point (in h/m units)
C//   TOUT  .Time to reach crystal exit
C//   G(3) ...Local diff. vector in the point of reflection
C//   P ... weight to be attributed to the event
C--------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'structures.inc'
c      INCLUDE 'crystal.inc'
      RECORD /CRYSTAL/ CR
      INTEGER*4 I,J,M,DIR,IH,IV
      REAL*8 SEC,EPS,ksi,Z,Z1,sigma,P0,PHI,alfh,alfv,GABS,G0(3)
      REAL*8 dTheta,P,DT,AA,V0(3)
      PARAMETER(SEC=4.85D-6,EPS=1D-5)
      REAL*8 R(3),K(3),K0,Q,TAU,TOUT,G(3),PP,A(3)
      REAL*4 RAN1,GASDEV1
      REAL*8 ERF,HMOS_DIST
      INTEGER*4 MSEG,NSEG
      PARAMETER (MSEG=50)
      REAL*8 TSEG(0:MSEG),ALPHA(0:MSEG),GRAD(0:MSEG),PSEG(0:MSEG),
     *       PSEG0(0:MSEG)
      COMMON /CRBORDERS/ TSEG,PSEG,PSEG0,ALPHA,GRAD,NSEG 
      
C/// add random vertical mosaic angle 
      PHI=CR.VMOS*GASDEV1(1,0.,3.)       
C/// trace through segments and find times, scatt. probabilities, angular dev. etc..
      call CR_BORDERS(CR,DIR,R,K,K0,Q,PHI) 
            
      dTheta=0
      TOUT=TSEG(NSEG)
      P=PSEG(NSEG)   ! scattering probability for the whole rest of the crystal
      IF(DIR.LT.0) P=1.  
      IF (P*PP.LT.EPS) GOTO 99 ! don't follow trajectories with low probability
      ksi=RAN1()
      Z=P*RAN1()
C Select segment for next reflection according to the scattering probabilities
      M=1
      DO WHILE (M.LE.NSEG.AND.Z.GE.PSEG(M))
	  M=M+1
      END DO

      
C No reflection in any segment
      IF(M.GT.NSEG) THEN
        IF(DIR.LT.0) THEN
	    GOTO 90
	ELSE    
	    GOTO 99
	ENDIF   
      ENDIF	  

C/// Find point of reflection inside the I-th segment      
C Mosaic crystal
      IF(ABS(GRAD(M)).LT.1.D-6*K0*Q) THEN
          P0=PSEG0(M)
          AA=ALPHA(M)/CR.HMOS
	  sigma=K0*Q*HMOS_DIST(AA)/CR.HMOS
	  DT=-Log(1-ksi*P0)/sigma

c        write(*,1) 'MOS : ',DT,grad(M),AA,sigma
C Mosaic crystal with deformation gradient
      ELSE
        IF(CR.HMOS.GT.SEC) THEN
          P0=PSEG0(M)
	  Z=GRAD(M)*Log(1-ksi*P0)/(K0*Q)
C ERF shoud be used only on (-inf.;0) because of num. precision
          AA=ALPHA(M)/CR.HMOS
	  IF(AA.GT.0) THEN
	    Z1=ERF(-AA,0)
	    IF(1.-Z1-Z.GT.0.5) THEN
	      DT=-ERF(Z1+Z,1)-AA
            ELSE
	      DT=ERF(1.-Z1-Z,1)-AA
	    ENDIF
          ELSE
	    Z1=ERF(AA,0)
	    IF(Z1-Z.GT.0.5) THEN
	      DT=-ERF(1.-Z1+Z,1)-AA
            ELSE
	      DT=ERF(Z1-Z,1)-AA
	    ENDIF	     
	  ENDIF
 	  DT=DT*CR.HMOS/GRAD(M)		     
	ELSE  ! no mosaic
	  DT=-ALPHA(M)/GRAD(M)
	  IF (DT.LE.0) GOTO 99
	ENDIF
      ENDIF
      TAU=DT+TSEG(M-1) 
            
C Move to the point of reflection
      DO I=1,3
         R(I)=R(I)+TAU*K(I)
      ENDDO

c1     format(a,5(G12.4,1x))
2     format(a,I3,2x,5(G12.4,1x))

C// Coordinates of the segment centre
      ih=nint((0.5+R(1)/CR.FRAME.SIZE(1))*CR.nh+0.5)
      iv=nint((0.5+R(2)/CR.FRAME.SIZE(2))*CR.nv+0.5)
      V0(1)=(2*ih-CR.nh-1)*CR.FRAME.SIZE(1)/2./CR.nh
      V0(2)=(2*iv-CR.nv-1)*CR.FRAME.SIZE(2)/2./CR.nv
      V0(3)=0.
c      write(*,2) 'WALK: ',ih,PSEG0(M),(R(I),I=1,3),K(3)
	
C Calculate local G-vector 
      DO I=1,3
        G0(I)=CR.G(I) 
        IF (CR.MAPG(I)) THEN
          A(I)=0.
	  DO J=1,3
	    A(I)=A(I)+CR.DG_DR(I,J)*(R(J)-V0(J))
	  ENDDO
	  G0(I)=G0(I)+A(I)
	ENDIF  
      END DO
      GABS=SQRT(G0(1)**2+G0(2)**2+G0(3)**2)  	    
C Add segment tilt angle and vertical mosaic spread
      alfh=V0(1)*CR.RH
      alfv=V0(2)*CR.RV
      G(1)=G0(1)-G0(3)*alfh   
      G(3)=G0(3)+G0(1)*alfh  
      IF (ABS(G0(3)).LT.0.01*GABS) THEN  ! symmetric Laue case          
         G(2)=G0(2)-G0(1)*alfv+GABS*PHI 
      ELSE
         G(2)=G0(2)-G0(3)*alfv+GABS*PHI  ! other case       
      ENDIF
C Add the angle of the mosaic block
      DO I=1,3
         G(I)=G(I)+GABS*CR.GAMA(I)*(ALPHA(M)+GRAD(M)*DT)
      END DO
      Z=GABS/SQRT(G(1)**2+G(2)**2+G(3)**2)
      DO i=1,3
        G(i)=G(i)*Z*DIR
      END DO

      PP=PP*P
      RETURN
      
      
90    TAU=TOUT
      PP=PP*P
      DO I=1,3
         R(I)=R(I)+TAU*K(I)
      ENDDO
      RETURN
      
99    PP=0
      RETURN
      
      
      END



C -------------------------------------------------------------
      SUBROUTINE CR_BORDER(CR,R,K,TIN,TOUT)
C     Returns times of intersection with the crystal borders, 
C     starting at current position R and measured along K.
C     All in crystal local coordinate.
C  !! Time is in units [sec*h/m] everywhere in NESS99 !!!   i.e. length=time*K
C--------------------------------------------------------------
      implicit none
      INCLUDE 'structures.inc'
c      INCLUDE 'crystal.inc'
      RECORD /CRYSTAL/ CR
        
      REAL*8 R(3),K(3),TIN,TOUT,T1(3),T2(3),DUM
      INTEGER*4 I
	
      DO I=1,3
      IF (ABS(K(I)).GT.1.0D-8) THEN
       	T2(I)=(CR.FRAME.SIZE(I)/2 - R(I))/K(I)
	T1(I)=(-CR.FRAME.SIZE(I)/2 - R(I))/K(I)
        IF (T1(I).GT.T2(I)) THEN
	  DUM=T1(I)
          T1(I)=T2(I)
	  T2(I)=DUM
        ENDIF   	 
      ELSE
	  T2(I)=1.0D30
	  T1(I)=-1.0D30
      ENDIF
      END DO
      TIN=MAX(T1(1),T1(2),T1(3))
      TOUT=MIN(T2(1),T2(2),T2(3))
      IF (TIN.GT.TOUT) THEN
	 TIN=1D30
	 TOUT=1D30
      ENDIF
      END            

C -----------------------------------------------------------------------
      SUBROUTINE CR_BORDER1(CR,R,K,TIN,TOUT,TSIN,TSOUT)
C     Returns times of intersection with the crystal borders, 
C     starting at current position R and measured along K.
C     All in crystal local coordinate.
C New! Extended for the case of migration from one segment to another one.  
C New! TSIN,TSOUT ... times of intersection with the nearest segment borders.
C Caution! Time is in units [sec*h/m] everywhere in NESS99 !!!   i.e. length=time*K
C------------------------------------------------------------------------
      implicit none
      INCLUDE 'structures.inc'
c      INCLUDE 'crystal.inc'
      RECORD /CRYSTAL/ CR
        
      REAL*8 R(3),K(3),TIN,TOUT,TSIN,TSOUT
      REAL*8 DX,X0,DUM,T1(3),T2(3),EPS
      PARAMETER (EPS=1D-8)
      INTEGER*4 I,IH,IV
C/// check for crystal intersectionss 
      DO I=1,3
      IF (ABS(K(I)).GT.EPS) THEN
        DX=CR.FRAME.SIZE(I)/2
       	T2(I)=(DX - R(I))/K(I)
	T1(I)=(-DX - R(I))/K(I)
        IF (T1(I).GT.T2(I)) THEN
	  DUM=T1(I)
          T1(I)=T2(I)
	  T2(I)=DUM
        ENDIF   	 
      ELSE
	  T2(I)=1.0D30
	  T1(I)=-1.0D30
      ENDIF
      END DO
      TIN=MAX(T1(1),T1(2),T1(3))
      TOUT=MIN(T2(1),T2(2),T2(3))
      IF (TIN.GT.TOUT) THEN
	 TIN=1D30
	 TOUT=1D30
      ENDIF
C/// determine segment coordinates 
      ih=nint((0.5+R(1)/CR.FRAME.SIZE(1))*CR.nh+0.5)
      iv=nint((0.5+R(2)/CR.FRAME.SIZE(2))*CR.nv+0.5)
C/// check for segment intersectionss 
      DO I=1,2
      IF (ABS(K(I)).GT.EPS) THEN
        if (I.EQ.1) then
           x0=ih-(CR.nh+1)/2.0
           dx=CR.FRAME.SIZE(I)/CR.nh
        else 
           x0=iv-(CR.nv+1)/2.0
           dx=CR.FRAME.SIZE(I)/CR.nv
        endif
       	T2(I)=((x0+0.5)*dx - R(I))/K(I)
       	T1(I)=((x0-0.5)*dx - R(I))/K(I)
        IF (T1(I).GT.T2(I)) THEN
	  DUM=T1(I)
          T1(I)=T2(I)
	  T2(I)=DUM
        ENDIF   	 
      ELSE
	  T2(I)=1.0D30
	  T1(I)=-1.0D30
      ENDIF
      END DO
      TSIN=MAX(T1(1),T1(2),T1(3))
      TSOUT=MIN(T2(1),T2(2),T2(3))
      IF (TSIN.GT.TSOUT) THEN
	 TSIN=1D30
	 TSOUT=1D30
      ENDIF

      if (ABS(TSOUT-TOUT).LE.EPS) TSOUT=TOUT


      END            

C     ----------------------------------------- 
      REAL*8 FUNCTION CRYSTAL_RC(CR,ETA1,ETA2)
C     ----------------------------------------- 
      IMPLICIT NONE
      INCLUDE 'structures.inc'
c      INCLUDE 'ness_const.inc'
c      INCLUDE 'crystal.inc'
      
      
      REAL*4 Z,LINTERP4,DUM
      REAL*8 ETA1,ETA2,SEC
      PARAMETER(SEC=4.8481368D-6)
      INTEGER*4 MMOS
      PARAMETER(MMOS=128)
      REAL*4 X_CR(MMOS,2),Y_CR(MMOS,2),T_CR(MMOS,2)
      INTEGER*4 N_CR
      COMMON /CR_REF/ X_CR,Y_CR,T_CR,N_CR
	
      RECORD /CRYSTAL/ CR 
        
c/// horizontal
      IF (CR.HMOS.LT.sec) THEN
	 Z=CR.REF
      ELSE IF (CR.NREF.LT.1.OR.CR.NREF.GT.2) THEN 
         DUM=(ETA1/CR.HMOS)**2
	 IF (DUM.GT.16) THEN
            Z=0
         ELSE
            Z=EXP(-0.5*DUM)
         ENDIF
      ELSE 
         IF (ABS(ETA1).GT.CR.DETA) THEN
            Z=0
         ELSE
	    DUM=ETA1
	    Z=LINTERP4(X_CR(1,CR.NREF),Y_CR(1,CR.NREF),N_CR,DUM)
	 ENDIF
      ENDIF      

c// vertical
      IF ((abs(ETA2).GT.sec).AND.(CR.VMOS.GT.sec)) THEN 
        DUM=(ETA2/CR.VMOS)**2
        IF (DUM.GT.16) THEN
          Z=0
        ELSE
          Z=Z*EXP(-0.5*DUM)
        ENDIF
      ENDIF
      CRYSTAL_RC=Z
        
      RETURN
      END           

C--------------------------------------------------------------
      REAL*8 FUNCTION CRYST_TILT(CR,R2,G)
C Add tilt angle + vertical mosaicity + d-gradient to G, returns |G|
C small-angle approx.   
C--------------------------------------------------------------
      IMPLICIT NONE
      INCLUDE 'structures.inc'
c      INCLUDE 'crystal.inc'
      RECORD /CRYSTAL/ CR
      REAL*8 R2(3),G(3),alfh,alfv,dah,dav,GABS,GABS0,Z,ETA2
      INTEGER*4 ih,iv,i
      REAL*8 V3XV3
      REAL*4 GASDEV1
      
      
C/// determine segment coordinates and tilt angles
      if (CR.nh.GT.1) THEN
         ih=int((0.5+R2(1)/CR.FRAME.SIZE(1))*CR.nh)+1
         dah=CR.FRAME.SIZE(1)*CR.RH/CR.nh        
         alfh=(-(CR.nh-1)/2+ih-1)*dah
      else
         ih=1
	 alfh=0
      endif	 
      if (CR.nv.GT.1) THEN	 	 
        iv=int((0.5+R2(2)/CR.FRAME.SIZE(2))*CR.nv)+1
        dav=CR.FRAME.SIZE(2)*CR.RV/CR.nv
        alfv=(-(CR.nv-1)/2+iv-1)*dav
      else
         iv=1
	 alfv=0
      endif
      
C/// calculate diffraction vector including the tilt angle       
      GABS0=SQRT(V3xV3(CR.G,CR.G))      
      G(1)=CR.G(1)-CR.G(3)*alfh   ! +CR.DG_DR(1,1)*R2(1)+CR.DG_DR(1,3)*R2(3)
      G(3)=CR.G(3)+CR.G(1)*alfh   ! +CR.DG_DR(3,1)*R2(1)+CR.DG_DR(3,3)*R2(3)
      IF (ABS(CR.G(3)).LT.0.01*GABS0) THEN              
           G(2)=CR.G(2)-CR.G(1)*alfv ! only in Laue case
      ELSE
           G(2)=CR.G(2)-CR.G(3)*alfv
      ENDIF
C/// add random vertical mosaic angle 
      ETA2=CR.VMOS*GASDEV1(1,0.,3.)       
      G(2)=G(2)+GABS0*ETA2            
C/// normalize: |G|=const.
      GABS=SQRT(V3xV3(G,G))
      Z=GABS0/GABS
      DO i=1,3
          G(i)=G(i)*Z
      ENDDO
      CRYST_TILT=GABS0
      
      END


C-----------------------------------------------------
      SUBROUTINE ERF_INIT(F,AMIN,AMAX)
C Calculate lookup tables for ERF function
C-----------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 I,J
      REAL*8 SUM,Z1,Z2,Z3,A,B,DET,X1I,XJ,AMIN,AMAX
      INTEGER*4 DIM
      PARAMETER(DIM=1025)
      REAL*8 XMIN,DX,Y(DIM),XMIN1,DX1,Y1(DIM)
      COMMON /ERFCOM/ XMIN,DX,Y,XMIN1,DX1,Y1      
      REAL*8 F
      EXTERNAL F

C Generate cumulative function from F(X)
      SUM=0.
      DX=(AMAX-AMIN)/(DIM-1)
      XMIN=AMIN
      Y(1)=0.
      DO I=1,DIM-1
        Z1=XMIN+(I-1)*DX
        Z2=Z1+DX/2
        Z3=Z1+DX  
        SUM=SUM+F(Z1)+4*F(Z2)+F(Z3)
        Y(I+1)=SUM
      END DO  
      DO I=1,DIM
        Y(I)=Y(I)/Y(DIM)
      ENDDO
C Generate inverse cumulative function
      DX1=1.D+0/(DIM-1)
      XMIN1=Y(1)
      Y1(1)=XMIN
      Y1(DIM)=XMIN+(DIM-1)*DX
      J=1
      DO I=2,DIM-1
        X1I=XMIN1+(I-1)*DX1
        DO WHILE (Y(J).LT.X1I.AND.J.LT.DIM-1)
          J=J+1
        END DO
10      XJ=XMIN+(J-1)*DX
        A=(Y(J+1)+Y(J-1)-2*Y(J))/2
        B=(Y(J+1)-Y(J-1))/2
        IF (ABS(A).LT.1D-30) THEN
           J=J-1
           GOTO 10
        ELSE
          DET=B**2-4*A*(Y(J)-X1I)
          IF (DET.LE.0) then
            write(*,*) 'Error in ERF_INIT: ',DET,A,B
            pause
          ENDIF
          Z1=XJ+DX*(-B+SQRT(DET))/2/A
          Z2=XJ+DX*(-B-SQRT(DET))/2/A
          IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
          Y1(I)=Z1
        ENDIF
      END DO
  
      END

C----------------------------------------------------------------
      REAL*8 FUNCTION ERF(ARG,INV)
c Return cumulative function (or inverse, if INV=0)
c Uses lookup table generated by CUM_INIT
C----------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 DIM,INV,J1,J2,J,I
      REAL*8 ARG,A,B,Z,XJ,DET,Z1,Z2,ARG1
      PARAMETER(DIM=1025)
      REAL*8 XMIN,DX,Y(DIM),XMIN1,DX1,Y1(DIM)
      REAL*8 ERF_INTERP
      COMMON /ERFCOM/ XMIN,DX,Y,XMIN1,DX1,Y1      

C Cumul. function
      IF (INV.NE.1) THEN  
        ERF=ERF_INTERP(XMIN,DX,Y,DIM,ARG)
C else Iverse cumul. function
      ELSE
        IF(ARG.LE.XMIN1) THEN
           ERF=Y1(1) ! left limit
        ELSE IF(ARG.GE.XMIN1+(DIM-1)*DX1) THEN
           ERF=Y1(DIM)  ! right limit
        ELSE
C Find J1,J2 so that Y(J1) > A >= Y(J2) and J2=J1+1
          ARG1=ARG
          IF(ARG.GT.9.D-1) ARG1=1.D0-ARG
          Z=(ARG1-XMIN1)/DX1+1
          I=INT(Z)
          Z1=(Y1(I)-XMIN)/DX+1
          Z2=(Y1(I+1)-XMIN)/DX+1   
          J1=INT(Z1)
          J2=INT(Z2)+1
          DO WHILE(J2.GT.J1+1)
            J=(J2+J1)/2
            IF(Y(J).LT.ARG1) THEN
               J1=J
            ELSE
               J2=J
            ENDIF
          ENDDO
C Set J so that Y(J) is close to ARG 
          J=J1
          IF(ARG1-Y(J1).GT.Y(J2)-ARG1) J=J2
C but avoid J=1 or J=DIM
          IF(J.LT.2) J=J+1
          IF(J.GT.DIM-1) J=J-1
C interpolate quadratically between Y(J-1) and Y(J+1) 
C return inverse value for Y=ARG
10        XJ=XMIN+(J-1)*DX
          A=(Y(J+1)+Y(J-1)-2*Y(J))/2
          B=(Y(J+1)-Y(J-1))/2
          IF(ABS(B).LT.1.D-30) THEN
             ERF=XJ
          ELSE IF (ABS(A).LT.1D-30) THEN
            ERF=XJ+Y(J+1)
            J=J-1
            GOTO 10
          ELSE
            DET=B**2-4*A*(Y(J)-ARG1)
            DET=SQRT(DET)
            Z1=XJ+DX*(-B+DET)/2/A
            Z2=XJ+DX*(-B-DET)/2/A
            IF (ABS(Z2-XJ).LT.ABS(Z1-XJ)) Z1=Z2
            IF(ARG.GT.9.D-1) Z1=-Z1
            ERF=Z1
          ENDIF
        ENDIF            
      ENDIF
      END 

C-------------------------------------------------------------------
      REAL*8 FUNCTION ERF_INTERP(XMIN,DX,Y,DIM,A)
C Quadratic interpolation in Y array with equidistant X
C given by XMIN,DX. MIN,MAX are values to be returned outside limits
C-------------------------------------------------------------------
      IMPLICIT NONE
      INTEGER*4 DIM,I
      REAL*8 XMIN,DX,A,Z,XMAX
      REAL*8 Y(DIM)   

      XMAX=XMIN+(DIM-1)*DX
      if (A.LE.XMIN) THEN 
         ERF_INTERP=Y(1)
      elseif (A.GE.XMAX) THEN         
         ERF_INTERP=Y(DIM)
      elseif (A.LE.XMIN+DX) THEN         
         ERF_INTERP=Y(1)+(Y(2)-Y(2))*(A-XMIN)/DX
      elseif (A.GE.XMAX-DX) THEN         
         ERF_INTERP=Y(DIM-1)+(Y(DIM)-Y(DIM-1))*(A-XMAX+DX)/DX
      else
        Z=(A-XMIN)/DX+1
        I=NINT(Z)
        if (I.EQ.DIM-1) I=I-1
        IF (I.EQ.2) I=3
        ERF_INTERP=Y(I)+(Y(I+1)-Y(I-1))/2*(Z-I)+
     *             (Y(I+1)+Y(I-1)-2*Y(I))/2*(Z-I)**2
      endif
      end

