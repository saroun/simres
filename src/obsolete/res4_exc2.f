C---------------------------------------
      SUBROUTINE EXCI(IOCONS,ICOM,Q,omexc,SQOM)
C---------------------------------------      
c
c     The eigenvectors have to be given in C-Convention!     
C
c
      implicit none
c
      integer j,i,jtr,ICOM,IOCONS
c      
      CHARACTER*60 fileph
      real*8 pi,pihalf,zero,z1,z2,z3
      real*8 q(4),omexc(6),tau(3),qphon(3),sqom(6),fqsq(6),bose,xt,temp
      real*8 qphon0(3),sigma1(3),sigma2(3),sigma3(3),
     *       f(6),e(6,3),omega0(6),domdq(6),b1(6),b2(6),b3(6)
      real*8 qe,qsq,qphsq,qphnrm,dq1,dq2,dq3,dqphi,omegaint(6)
c
      common /exptr/ tau,qphon0,sigma1,sigma2,sigma3,
     *         f,e,omega0,domdq,b1,b2,b3,temp
      
      CHARACTER*20 phonname,datname,resname 
      COMMON /FNAMES/ phonname,datname,resname 
      PARAMETER (PI=3.141592653,pihalf=1.570796327,zero=.0)
c
c
c///  Initialization part (ICOM=0):

      IF(ICOM.EQ.0) THEN

2     format(' Phonon parameter file [phonon.par] : ',$)
3     format(a)
4     format(a20)
      write(IOCONS+1,2)       
      read(IOCONS,3) fileph
      write(phonname,4) fileph(1:20)
      open(1,name=fileph,err=50,status='old')
      goto 52
50    open(1,name='phonon.par',err=999,status='old')
      write(phonname,4) 'phonon.par            '
      write(IOCONS+1,*) 'Parameters read from default file phonon.par'
52    read (1,*,err=998)
      read (1,*,err=998) temp
      read (1,*,err=998) tau,qphon0
      read (1,*,err=998) sigma1
      read (1,*,err=998) sigma2
      read (1,*,err=998) sigma3
      do j=1,6
        read (1,*,err=998) f(j),e(j,1),e(j,2),e(j,3),
     *        omega0(j),domdq(j),b1(j),b2(j),b3(j)
      enddo
      close(1)

C///  Supposing orthogonal lattice  !
C///  normalizing sigma1,2,3:

      z1=sigma1(1)**2+sigma1(2)**2+sigma1(3)**2
      z2=sigma2(1)**2+sigma2(2)**2+sigma2(3)**2
      z3=sigma3(1)**2+sigma3(2)**2+sigma3(3)**2
      do i=1,3
         if (z1.ne.0) sigma1(i)=sigma1(i)/sqrt(z1)
         if (z2.ne.0) sigma2(i)=sigma2(i)/sqrt(z2)
         if (z3.ne.0) sigma3(i)=sigma3(i)/sqrt(z3)
      end do         
      
C///  normalizing e1,2,3:

      do j=1,6
      z1=e(j,1)**2+e(j,2)**2+e(j,3)**2
      z2=e(j,1)**2+e(j,2)**2+e(j,3)**2
      z3=e(j,1)**2+e(j,2)**2+e(j,3)**2
      do i=1,3
         if (z1.ne.0) e(j,i)=e(j,i)/sqrt(z1)
         if (z2.ne.0) e(j,i)=e(j,i)/sqrt(z2)
         if (z3.ne.0) e(j,i)=e(j,i)/sqrt(z3)
      end do         
      end do
      
      RETURN
998   write(IOCONS+1,*) 'Cannot read phonon parameters !'
      RETURN
999   write(IOCONS+1,*) 'Cannot open file with phonon parameters !'
      RETURN                
      
      ELSE       ! End of initializtion

C///  Returns six energies and structure factors:
C///  Orthogonal lattice is supposed !

      do 11 i=1,3
11        qphon(i) = q(i)-tau(i) 
      qsq = q(1)**2+q(2)**2+q(3)**2
      qphsq = qphon(1)**2+qphon(2)**2+qphon(3)**2
      qphnrm = sqrt(qphsq)
c
      dq1 = .0
      dq2 = .0
      dq3 = .0
      do i=1,3
           dqphi = qphon(i)-qphon0(i)
         dq1 = dq1+dqphi*sigma1(i)
         dq2 = dq2+dqphi*sigma2(i)
         dq3 = dq3+dqphi*sigma3(i)
      enddo
c
      do j=1,6
        omegaint(j) = omega0(j)+dq1*domdq(j)+
     *                dq1**2*b1(j)+dq2**2*b2(j)+dq3**2*b3(j)
      enddo
c
c *** calculate structure factor
c
      do j=1,6
        qe = q(1)*e(j,1)+q(2)*e(j,2)+q(3)*e(j,3)
        fqsq(j) = (f(j)*qe)**2
      enddo
c 
c *** Bose factor [meV, K] ***
c      
        
      do j=1,6
        jtr = NINT(sign(1.D0,omega0(j)))      
        omegaint(j)=abs(omegaint(j))        
        IF(omegaint(j).LT.0.0001) omegaint(j)=0.0001
        IF(temp.le.0.) temp=0.01
        xt = 11.55*omegaint(j)/temp
        IF (XT.LE.0) THEN
           bose=0.
        ELSE IF (XT.GT.10) THEN
           bose = .5*(JTR+1.)
        ELSE 
           bose =.5*(JTR+1.)+1./(EXP(XT)-1.)
        ENDIF
        if(fqsq(j).eq.0) then
          sqom(j)=0.
        else  
          sqom(j)  = bose*(fqsq(j)/omegaint(j))
        endif  
        omexc(j) = jtr*omegaint(j)
      enddo
      
      ENDIF
          
C
c *** End! ************************
c
      RETURN
      END
c

