C//////////////////////////////////////////////////////////////////////
C////                                                              //// 
C////  NEutron Scattering Simulation - v.2.0, (c) J.Saroun, 2000   ////
C////                                                              //// 
C//////////////////////////////////////////////////////////////////////
C////
C////  Subroutines describing objects - BENDER
C////  
C////                          
C//////////////////////////////////////////////////////////////////////

C-------------------------------------------------------------------
        SUBROUTINE QUADREQ(A,B,C,X)
C Solve quadratic equation A*X^2 + B*X + C = 0
C Try to find a solution > EPS=1E-10
C 1) no solution .. return 10^30
C 2) 1 solution  .. return this
C 3) 2 solutions .. return the smaller one
C-------------------------------------------------------------------
      IMPLICIT NONE
      REAL*8 EPS
      PARAMETER (EPS=1.D-10)
      REAL*8 A,B,C,X,X1,X2,DET

      IF (A.EQ.0) THEN          
        IF (ABS(B).LT.EPS) THEN 
          GOTO 20
        ELSE
          X=-C/B
          GOTO 30
        ENDIF
      ELSE          
        DET=B**2-4*A*C
        IF (DET.EQ.0.) THEN
           X=-B/2./A
           GOTO 30
        ELSE IF (DET.LT.0.) THEN
           GOTO 20
        ELSE
           DET=SQRT(DET) 
           X1=(-B+DET)/2./A
           X2=(-B-DET)/2./A
           IF (X1.GT.EPS.AND.X2.GT.EPS) THEN
              X=MIN(X1,X2)
           ELSE IF (X1.GT.EPS) THEN
              X=X1
           ELSE IF (X2.GT.EPS) THEN
              X=X2
           ELSE 
              GOTO 20
           ENDIF
        ENDIF
      ENDIF 
                  
30    IF (X.LT.EPS) GOTO 20
      RETURN  

20    X=1.D30

      END
	
C--------------------------------------------------------
      LOGICAL*4 FUNCTION BENDER_PASS(BENDER,R,IH,IV,ex)
c Checks, whether neutron fits inside any slit and 
c returns slit indices
C--------------------------------------------------------	
      IMPLICIT NONE
      INCLUDE 'structures.inc'
      INTEGER*4 IH,IV,ex
      LOGICAL*4 LOG1
      REAL*8 R(3),jh,jv,W,H
      RECORD /BENDER/ BENDER
      REAL*8 DJ
      COMMON /OSCBEND/ DJ

      if (ex.EQ.1) then 
         W=BENDER.W2
         H=BENDER.H2
      else
         W=BENDER.FRAME.SIZE(1)
         H=BENDER.FRAME.SIZE(2)
      endif
      jh=(R(1)/W+0.5-DJ)*BENDER.NLH
      jv=(R(2)/H+0.5)*BENDER.NLV
      IH=NINT(jh-0.5)
      IV=NINT(jv-0.5)
      LOG1=((ABS(jh-NINT(jh))*W/BENDER.NLH.GE.BENDER.DLH/2.).AND.
     1      (ABS(jv-NINT(jv))*H/BENDER.NLV.GE.BENDER.DLV/2.).AND.
     1      (jh.GT.0.).AND.
     1      (jv.GT.0.).AND.
     1      (jh.LT.BENDER.NLH).AND.
     1      (jv.LT.BENDER.NLV))
              
      BENDER_PASS=LOG1
      END  

        
C	----------------------------------------------------------------------
	LOGICAL*4 FUNCTION CONTACT(X0,Z0,alpha,RO,kx,kz,T)
c	finds next contact with a lamella and returns time T needed 
C       to reach it
C       X0,Z0 ... initial coordinates of the neutron with respect
C                 to the lamella front end
C       alpha ... lamella inclination angle (without curvature)
C       R0    ... lamella curvature
C       kx,kz ... transversal and longitudinal components of neutron 
C                 k vector
C	-----------------------------------------------------------------------	
        IMPLICIT NONE
        LOGICAL*4 LOG1, QUADREQ
        REAL*8 X0,Z0,alpha,RO,kx,kz,T,A,B,C

        A=0.5*ro*kz**2
        B=alpha*kz-kx+ro*kz*z0
        C=alpha*z0+0.5*ro*z0**2-x0
        LOG1=QUADREQ(A,B,C,T)
        CONTACT=LOG1
        END



C	---------------------------------------------------
	REAL*8 FUNCTION BENDER_REF(ID,BENDER,Q,S)
c	returns reflectivity for given momentum transfer Q
C       ID identifies which surface is touched
C       ID=0  left 
C       ID=1  right 
C       ID=2  top 
C       ID=3  bottom 
C	---------------------------------------------------	
        IMPLICIT NONE
        INCLUDE 'structures.inc'
        INTEGER*4 m_n(5),iz,ID,NR
        CHARACTER*3 m_name(5)
        REAL*8 Q,S,z,dQ,Q1,gamma,R	
        REAL*8 m_alpha(128,5), m_ref1(128,5),m_ref2(128,5)
      
        RECORD /BENDER/ BENDER
        COMMON /MIRROR/ m_alpha,m_ref1,m_ref2,m_n,m_name    

        if (ID.EQ.0) then 
          if (S.ge.0) then
            gamma=BENDER.GHLU
            R=BENDER.RHLU
            NR=BENDER.NHLU
          else  
            gamma=BENDER.GHLD
            R=BENDER.RHLD
            NR=BENDER.NHLD
          endif 
        else if (ID.EQ.1) then 
          if (S.ge.0) then
            gamma=BENDER.GHRU
            R=BENDER.RHRU
            NR=BENDER.NHRU
          else  
            gamma=BENDER.GHRD
            R=BENDER.RHRD
            NR=BENDER.NHRD
          endif 
        else if (ID.EQ.2) then 
            gamma=BENDER.GVT
            R=BENDER.RVT
            NR=BENDER.NVT
        else if (ID.EQ.3) then 
            gamma=BENDER.GVB
            R=BENDER.RVB
            NR=BENDER.NVB
        else
            gamma=0.   
            R=0. 
        endif    
        
        if (gamma.le.0.OR.Q.LE.0) then
           BENDER_REF=0
           RETURN
        endif   
        
	if(NR.LE.0.OR.NR.GT.5) then
              if (Q.ge.0.AND.Q.lt.2*PI*gamma) then	                         
                 BENDER_REF=R
              else
                 BENDER_REF=0                 
              endif 
        else
              Q1=Q/2/PI/GammaNi
              dQ=m_alpha(2,NR)-m_alpha(1,NR)
              z=(Q1-m_alpha(1,NR))/dQ
              
              iz=INT(z)+1
              if(z.LT.0.OR.z.GE.m_n(NR).OR.iz.GE.m_n(NR)) then 
                 BENDER_REF=0                 
              else if (S.GE.0) then
               BENDER_REF=m_ref1(iz,NR)+(z-iz+1)*
     1                (m_ref1(iz+1,NR)-m_ref1(iz,NR))  
              else if (S.LT.0) then
               BENDER_REF=m_ref2(iz,NR)+(z-iz+1)*
     1                (m_ref2(iz+1,NR)-m_ref2(iz,NR))
              else     
                 BENDER_REF=0 
              endif 
        endif                         
	END
	
C     -------------------------------
      SUBROUTINE BENDER_INIT(BENDER)
C     -------------------------------	
      IMPLICIT NONE
      
      INCLUDE 'structures.inc'
      RECORD /BENDER/ BENDER
      
      CALL SLIT_INIT(BENDER.FRAME)
      BENDER.TYP=0
      IF (BENDER.CURV.NE.0) THEN 
         BENDER.TYP=BENDER.TYP+1
      ENDIF   
      IF (BENDER.GHLU.NE.0.OR.
     *     BENDER.GHLD.NE.0.OR.
     *     BENDER.GHRU.NE.0.OR.
     *     BENDER.GHRD.NE.0.OR.
     *     BENDER.GVT.NE.0.OR.
     *     BENDER.GVB.NE.0) THEN 
        BENDER.TYP=BENDER.TYP+2
      ENDIF
      END  
      


CxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCxxxxxxxxxCx

C----------------------------------------------------
      LOGICAL*4 FUNCTION BENDER_GO(BENDER,NEUI,NEUF)
C----------------------------------------------------	
      IMPLICIT NONE

      INCLUDE 'structures.inc'

      
      RECORD /BENDER/ BENDER
      RECORD /NEUTRON/ NEUI,NEUF
      LOGICAL*4 LOG1, BENDER_PASS,CONTACT
      REAL*8 BENDER_REF
      INTEGER*4 IH,IV,IH1,IV1,I
      REAL*8 V(3),K(3),R(3),R2(3)
      REAL*8 AL,AR,AT,AB,TL,TR,TT,TB,ZL,ZR,ZT,ZB,XL,XR,XT,XB
      REAL*8 KK,DUM,BETA0,DELTA0,P,PP,DT,T,Q
      REAL*8 DJ
      REAL*4 RAN1
      COMMON /OSCBEND/ DJ

      dj=0
D      logical*4 rep

D      CHARACTER*4 dname
D      CHARACTER*1 CH    
D 100    format('tra',I3,a4,'.dat')
D 101    format(a,2x,7(F11.5,2x),a1)
D 102    format(1x,5(E10.3,2x),2x,a10)
D 103    format(1x,3(E10.3,2x),2x,I4,2x,a10)

c      rep=(bender.frame.size(3).eq.240.)      
c      rep=(rep.and.bender.frame.count.gt.1000)
      
c      dname=bender.frame.name(1:4)
c      rep=(dname.eq.'col2'.and.bender.typ.gt.0)
      
D      IF(REP) WRITE(*,101) (NEUI.R(I),I=1,3),(NEUI.K(I),I=1,3),NEUI.P
      
C Convert to local coordinate system and move to the entry
      call SLIT_PRE(BENDER.FRAME,NEUI.R,NEUI.K,V,K)
      DO I=1,2
         R(I)=V(I)-V(3)/K(3)*K(I)
      ENDDO
      R(3)=0.
      P=1.D0
      T=-V(3)/HOVM/K(3)
      
      IF (BENDER.TYP.LT.0.OR.BENDER.FRAME.SIZE(3).LE.0) GOTO 700  ! collimator ignored
      
      dj=0
      if (BENDER.OSCILATE.GT.0) dj=(RAN1()-0.5)/BENDER.NLH 
      dj=dj*MAX(BENDER.FRAME.SIZE(1),BENDER.W2)
      
      IF (BENDER.TYP.EQ.0) THEN
         CALL SOLLER_GO(BENDER,R,K,P,T)
      ELSE IF (BENDER.TYP.EQ.1) THEN
         CALL GUIDE_GO(BENDER,R,K,P,T)
      ELSE IF (BENDER.TYP.EQ.2) THEN
         CALL PGUIDE_GO(BENDER,R,K,P,T)
      ELSE 
         GOTO 300
      ENDIF      
      IF (P.LE.0.D0) GOTO 300
      
      
            
      IF (BENDER.FRAME.SIZE(3).LE.0) GOTO 200  ! collimator ignored
      
      PP=1.
      BETA0=0.

C///  check the pass through the entry

c      IF(REP) then
c        write(*,101) (NEUI.R(i),i=1,3),(NEUI.K(i),i=1,3)
c        write(*,101) (BENDER.FRAME.ROT(1,i),i=1,3),BENDER.FRAME.AXI 
c        write(*,101) (R(i),i=1,3),(K(i),i=1,3),NEUI.P 
c        write(*,*) BENDER_PASS(BENDER,R,IH,IV,0), IH,IV
c        pause  
c      ENDIF 

      dj=0
      if (BENDER.OSCILATE.GT.0) dj=(RAN1()-0.5)/BENDER.NLH 
      IF (.NOT.BENDER_PASS(BENDER,R,IH,IV,0)) GOTO 300      
            
      BETA0=BENDER.CURV*BENDER.FRAME.SIZE(3)


D      CH='I'
D      write(*,103) (R(I),I=1,3),BENDER.FRAME.COUNT,
D    1 BENDER.FRAME.NAME 
D        write(dname,100) BENDER.FRAME.COUNT,BENDER.FRAME.NAME
D        if (dname(4:4).eq.' ') dname(4:4)='0'
D        if (dname(5:5).eq.' ') dname(5:5)='0'
D 	 Open(Unit=15,File=dname,Status='Unknown')
D        write(15,101) (R(i),i=1,3),(K(i),i=1,3),NEUI.P*PP,CH 

D      if (rep) then
D        write(*,*) IH,IV
D        write(*,101) dname,(R(i),i=1,3),(K(i),i=1,3),NEUI.P*PP,CH 	
D      endif
	


C ****  neutron guide  ****
      
      KK=SQRT(K(1)**2+K(2)**2+K(3)**2)

C at the guide entry:
C Z ... coordinate of the lamella center rel. to the guide center [units of guide width/height]
C X ... coordinate of the lamella wall [mm] rel. to the guide center
C A ... TAN(a), a=angle of the lamella w.r.t. guide axis

C  left      
      ZL=(IH+1)*1./BENDER.NLH-0.5+DJ
      XL=ZL*BENDER.FRAME.SIZE(1)-0.5*BENDER.DLH
      AL=(BENDER.W2-BENDER.FRAME.SIZE(1))/BENDER.FRAME.SIZE(3)*ZL
C  right      
      ZR=IH*1./BENDER.NLH-0.5+DJ
      XR=ZR*BENDER.FRAME.SIZE(1)+0.5*BENDER.DLH
      AR=(BENDER.W2-BENDER.FRAME.SIZE(1))/BENDER.FRAME.SIZE(3)*ZR
C  top      
      ZT=(IV+1)*1./BENDER.NLV-0.5
      XT=ZT*BENDER.FRAME.SIZE(2)-0.5*BENDER.DLV
      AT=(BENDER.H2-BENDER.FRAME.SIZE(2))/BENDER.FRAME.SIZE(3)*ZT
C  bottom      
      ZB=IV*1./BENDER.NLV-0.5
      XB=ZB*BENDER.FRAME.SIZE(2)+0.5*BENDER.DLV
      AB=(BENDER.H2-BENDER.FRAME.SIZE(2))/BENDER.FRAME.SIZE(3)*ZB
      
c      if (rep) then
c         write(*,101) 'left:   ',ZL,XL,AL,(R(1)-XL)
c         write(*,101) 'right:  ',ZR,XR,AR,(R(1)-XR)
c         write(*,101) 'top:    ',ZT,XT,AT,(R(2)-XT)
c         write(*,101) 'bottom: ',ZB,XB,AB,(R(2)-XB)
c         write(*,101) 'times:  ', 
c     &  (R(1)-XL-AL*R(3))/(AL*K(3)-K(1)),
c     &  (R(1)-XR-AR*R(3))/(AR*K(3)-K(1)),
c     &  (R(2)-XT-AT*R(3))/(AT*K(3)-K(2)), 
c     &  (R(2)-XB-AB*R(3))/(AB*K(3)-K(2))
c      endif

C***** beginning of the guide tracing cycle
50    CONTINUE
      
      TL=1.D35
      TR=1.D35
      TT=1.D35
      TB=1.D35
C get point (time) of contact with slit walls      
      IF (CONTACT(R(1)-XL,R(3),AL,BENDER.CURV,K(1),K(3),T)) TL=T       
      IF (CONTACT(R(1)-XR,R(3),AR,BENDER.CURV,K(1),K(3),T)) TR=T            
      IF (CONTACT(R(2)-XT,R(3),AT,0.D0,K(2),K(3),T)) TT=T       	 
      IF (CONTACT(R(2)-XB,R(3),AB,0.D0,K(2),K(3),T)) TB=T 
	     	 
C select the closest one      
      DT=MIN(TL,TR,TT,TB)      
      IF(DT.EQ.1.D35) THEN ! no contact, send the neutron outside
         DT=1.5*BENDER.FRAME.SIZE(3)/K(3)
      ENDIF
c         DT=0.
c	 write(*,*) 'WARNING: unexpected event in '//BENDER.FRAME.NAME
c         write(*,101) (R(i),i=1,3),(K(i),i=1,3)
c	 write(*,*) 'IH,IV: ',IH,IV
c	 goto 300
c	 pause 
c      ENDIF	 
	 
C go to a point of reflection             
      DO I=1,3
         R2(I)=R(I)  ! remember old position !
         R(I)=R(I)+K(I)*DT  ! go to the contact point
      END DO   

c      if (rep) write(*,101) 'times:  ',TL,TR,TT,TB
c      if (rep) write(*,101) 'DT, R2  ',DT,R
c      if (rep) pause
      
      
      IF (R(3).GT.BENDER.FRAME.SIZE(3)) THEN  ! no contact, go to guide exit
         DT=(BENDER.FRAME.SIZE(3)-R2(3))/K(3)
	 DO I=1,3
               R(I)=R2(I)+DT*K(I)
         END DO 
	 NEUF.T=NEUF.T+DT/hovm  
c	 CH='E'
D         write(15,101) (R(i),i=1,3),(K(i),i=1,3),NEUI.P*PP,CH 
c       if (BENDER.FRAME.NAME(1:4).EQ.'col2') then  
D       if (rep) then  
D	 write(*,101) dname,(R(i),i=1,3),(K(i),i=1,3),PP,CH
D	 pause
D       endif 
	 GOTO 199  
      ENDIF 
      
      NEUF.T=NEUF.T+DT/hovm 
      
C Q=k*grazing angle             
D      CH='X'       
      IF (DT.EQ.TL) THEN
	  Q=K(1)-(AL+R(3)*BENDER.CURV)*kk
          P=BENDER_REF(0,BENDER,Q,NEUI.S)
	  IF (P.GT.0) K(1)=K(1)-2*Q  
D	  CH='L'          
      ELSE IF (DT.EQ.TR) THEN
	  Q=-K(1)+(AR+R(3)*BENDER.CURV)*kk
          P=BENDER_REF(1,BENDER,Q,NEUI.S)
	  IF (P.GT.0) K(1)=K(1)+2*Q 
D	  CH='R'             
      ELSE IF (DT.EQ.TT) THEN
	  Q=K(2)-AT*kk
          P=BENDER_REF(2,BENDER,Q,NEUI.S)
	  IF (P.GT.0) K(2)=K(2)-2*Q 
D 	  CH='T'
      ELSE IF (DT.EQ.TB) THEN
	  Q=-K(2)+AB*kk
          P=BENDER_REF(3,BENDER,Q,NEUI.S)
	  IF (P.GT.0) K(2)=K(2)+2*Q 
D	  CH='B'             
      ENDIF
      PP=PP*P
      IF (PP.LT.1.D-4) PP=0
      IF (PP.GT.0) then
        DUM=SQRT(K(1)**2+K(2)**2+K(3)**2)
        DO I=1,3
          K(I)=K(I)*KK/DUM
        ENDDO

D        if(rep) write(*,101) dname,(r(i),i=1,3),(k(i),i=1,3),pp,ch 

	GOTO 50
      ENDIF 

D      if(rep) then
D        write(*,101) dname, (r(i),i=1,3),(k(i),i=1,3),pp,ch 
D        write(*,101) dname, q,dt,(r2(i),i=1,3)
D        write(*,*) 'no way !' 
D        pause  
D      endif 

      
D     CLOSE(15)
D     if (PP.GT.0) then
D         write(*,*) dname, 'error' 
D         pause
D     endif

      GOTO 300
     
199   CONTINUE
D     CLOSE(15)
D     if (PP.GT.0) then
D         write(*,*) dname 
D         pause
D     endif

200   CONTINUE
	 
c      if (BENDER.FRAME.NAME.EQ.'guide') then
c      if (ABS(abs(R(2))-BENDER.FRAME.SIZE(2)/2.).LE.BENDER.DLV) THEN       
c      if (BENDER.FRAME.COUNT.GT.40000) THEN       
c         write(*,101) (R(i),i=1,3),(K(i),i=1,3)
c	 pause 
c      endif
c      endif
      IF (BETA0.NE.0) THEN  !  correction for beam deflection
         DELTA0=0.5*BENDER.CURV*BENDER.FRAME.SIZE(3)**2        
         R2(1)=R(1)-R(3)*BETA0+DELTA0
         R2(3)=R(3)+R(1)*BETA0
         R2(2)=R(2)
	 K(1)=K(1)-K(3)*BETA0
         K(3)=K(3)+K(1)*BETA0
         DUM=SQRT(K(1)**2+K(2)**2+K(3)**2)
         DO I=1,3                        ! renormalize k
              K(I)=K(I)*KK/DUM
              R(I)=R2(I)
         END DO                 
      ENDIF
      DO I=1,3
           NEUF.R(I)=R(I)
           NEUF.K(I)=K(I)
      END DO    	   
      NEUF.P=NEUI.P*PP
      NEUF.S=NEUI.S
      BENDER.FRAME.COUNT=BENDER.FRAME.COUNT+1      
      BENDER_GO=.TRUE.
      RETURN

700   DO I=1,3
         NEUF.R(I)=R(I)
         NEUF.K(I)=K(I)
      END DO    	   
      NEUF.P=NEUI.P*P
      NEUF.T=NEUI.T+T
      NEUF.S=NEUI.S
      BENDER.FRAME.COUNT=BENDER.FRAME.COUNT+1      
      BENDER_GO=.TRUE.
      RETURN

300   CONTINUE
      BENDER_GO=.FALSE.
      NEUF.P=0 	
      RETURN
      END        

C--------------------------------------------------------------
      SUBROUTINE SOLLER_GO(BENDER,R,K,P,T)
C GO procedure for a simple collimator (non reflecting, TYP=0)      
C INPUT:   assume R,K at the entry in local coordinates !
C RETURN:  R,K at the exit in local coordinates, P=P*transmission, T=T+passage time
C--------------------------------------------------------------	
      IMPLICIT NONE
      INCLUDE 'structures.inc'
      
      RECORD /BENDER/ BENDER
      REAL*8 R(3),K(3),P,T
      LOGICAL*4 LOG1, BENDER_PASS
      INTEGER*4 IH,IV,IH1,IV1,I

C  check passage through the entry
      LOG1=BENDER_PASS(BENDER,R,IH,IV,0)
      IF (.NOT.LOG1) GOTO 100                  
C  move to the exit
      DO I=1,2
	 R(I)=R(I)+BENDER.FRAME.SIZE(3)/K(3)*K(I)
      ENDDO        
      R(3)=BENDER.FRAME.SIZE(3)
      T=T+R(3)/HOVM/K(3)
C  check passage through the same slit at the exit
      LOG1=(LOG1.AND.BENDER_PASS(BENDER,R,IH1,IV1,1))       
      IF ((.NOT.LOG1).OR.(IH1.NE.IH).OR.(IV1.NE.IV)) GOTO 100
      RETURN
      
C no passage
100   P=0.D0
      END        






C//////////////////  End of definition - BENDER  ///////////////////
