!C//////////////////////////////////////////////////////////////////////
!C////  $Id: niti.f,v 1.8 2019/08/15 15:02:09 saroun Exp $
!C////
!C////  R E S T R A X - Simulation of neutron three-axis spectrometers
!C////
!C////     Copyright (C) 1995-2018, All rights reserved
!C////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!C////     Institut Laue Langevin, Grenoble, France
!C////
!C////     Written by:  Jan Saroun
!C////     $Revision: 1.8 $
!C////     $Date: 2019/08/15 15:02:09 $
!C//////////////////////////////////////////////////////////////////////
!C////
!C////  N E S S - Ray-tracing package for  RESTRAX
!C////
!C////  NiTi supermirror data and functions
!C////
!C//////////////////////////////////////////////////////////////////////
      module NITI
      use CONSTANTS
      IMPLICIT NONE
      private
      SAVE
! Ni, Ti material data
      ! scattering cross-sections for Ni(1) and Ti(2)
      real(kind(1.D0)), parameter :: sigc(2) = (/1.208, 0.081/)  ! 1/cm
      real(kind(1.D0)), parameter :: sigi(2) = (/0.475, 0.1631/)  ! 1/cm
      real(kind(1.D0)), parameter :: siga(2) = (/0.228, 0.1921/)  ! 1/cm/A
      real(kind(1.D0)), parameter :: siga_58Ni = 1.024*siga(1)
      real(kind(1.D0)), parameter :: sigi_58Ni = 0.D0
      real(kind(1.D0)), parameter :: sigc_58Ni = 1.96*sigc(1)

      ! scattering length density
      real(kind(1.D0)), parameter :: rho(2) = (/9.413D-6, -1.908D-6/) ! 1/A^2
      real(kind(1.D0)), parameter :: rho_58Ni = 1.398*rho(1) ! 1/A^2
      ! Bragg edge
      real(kind(1.D0)), parameter :: bge(2) = (/4.07, 5.111/) ! A

      logical :: incl_coh = .true.
      logical :: monitor = .false. ! true if we need to monitor process: capture, scattering, ...

      public NITI_PROB, NITI_GETD, NITI_GETSIG, NITI_SET_MODE
      ! public NITI_PROB_DBG

      contains

!---------------------------------------------------------
      subroutine NITI_SET_MODE(adjRef, isMonitor)
! Temporary setting of the reflectivity model.
! Procedures calling NITI_PROB should make this setting first.
! Input:
! adjRef:  true to adjust reflectivity for NiTi coherent scattering:
!   (1) reflectivity is multiplied by attenuation due to coherent scattering
!   (2) coherent scattering is considered as an attenuation process
! isonitor: true if the model calculates probabilities for scattering,
! capture and free transmission
!---------------------------------------------------------
      logical,intent(in) :: adjRef, isMonitor
      incl_coh = adjRef
      monitor = isMonitor
      end subroutine NITI_SET_MODE


!---------------------------------------------------------
      subroutine NITI_GETD(m, d)
! Multilayer thickness, empirical function fitted on SN supermirror data.
! m = Multilayer m-value
! Returns in d:
! total material thickness in um: (1) Ni, (2)  Ti, (3) all
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: m
      REAL(KIND(1.D0)),intent(out) :: d(3)
      REAL(KIND(1.D0)), parameter :: Tot(3) = (/0.41647D0, 0.796D0, 1.953D0/)
      REAL(KIND(1.D0)), parameter :: Ni(3) = (/0.31252D0, 0.42611D0, 1.92135D0/)
      if (m>1.2) then
        d(3) = Tot(1) + Tot(2)*(m-1)**Tot(3)
        d(1) = Ni(1) + Ni(2)*(m-1)**Ni(3)
        d(2) = d(3)-d(1)
      else ! pure Ni
        d=(/0.1D0, 0.0D0, 0.1D0/)
      endif
      end subroutine NITI_GETD

!---------------------------------------------------------
      subroutine NITI_GETSIG(d, sig)
! Multilayer cross-sections averaged over total layer thicknesses
! d = array(3) with thickness in um: (1) Ni, (2) Ti, (3) All
! Returns in sig:
! scattering cross-section in 1/cm: (1) capture, (2) incoherent, (3) coherent
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: d(3)
      REAL(KIND(1.D0)),intent(out) :: sig(3)
      sig(1) = (d(1)*siga(1) + d(2)*siga(2))/d(3)
      sig(2) = (d(1)*sigi(1) + d(2)*sigi(2))/d(3)
      sig(3) = (d(1)*sigc(1) + d(2)*sigc(2))/d(3)
      end subroutine NITI_GETSIG

!---------------------------------------------------------
!      REAL(KIND(1.D0)) function  NITI_COH(d, lambda)
! Coherent scattering contribution in 1/cm
! Very primitive model assuming 1 or 0 below and above Bragg edge
! d = array(3) with thickness in um: (1) all, (2) Ni, (3) Ti
! lambda = wavelength [A]
!---------------------------------------------------------
!      REAL(KIND(1.D0)),intent(in) :: d(3), lambda
!      REAL(KIND(1.D0)) :: res
!      res = 0.D0
!      if (incl_coh) then
!        if (lambda<bge(1)) res = res + d(1)*sigc(1)
!        if (lambda<bge(2)) res = res + d(2)*sigc(2)
!      endif
!      NITI_COH = res/d(3)
!      end function  NITI_COH

!---------------------------------------------------------
      subroutine  NITI_TOT(lambda, pcoh, sigt)
! Total removal cross-section [1/cm]
! lambda = wavelength [A]
! pcoh = coherent scattering factor
! Returns in sigt:
! Total removal cross-section for Ni (1) and Ti (2)
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: lambda, pcoh
      REAL(KIND(1.D0)),intent(out) :: sigt(2)
      integer :: i
      REAL(KIND(1.D0)) :: sc
!1     format('NITI_TOT i, pcoh, lambda, sigi, siga, sigc, sigt: ',7(1x, G12.4))
      sigt=0.D0
      do i=1,2
        sc = 0.D0
        sigt(i) = sigi(i) + lambda*siga(i)
        if (incl_coh) then
          if (lambda<bge(i)) sc = pcoh*sigc(i)
        endif
        sigt(i) = sigt(i) + sc
        !write(*,1) i, int(pcoh), lambda, sigi(i),lambda*siga(i), sc, sigt(i)
      enddo
      end subroutine  NITI_TOT

!---------------------------------------------------------
!      REAL(KIND(1.D0)) function  NITI_AVE(d, lambda, pcoh)
! Total removal cross section averaged over materials in 1/cm
! d = array(3) with thickness in um: (1) all, (2) Ni, (3) Ti
! lambda = wavelength [A]
! pcoh = coherent scattering factor
!---------------------------------------------------------
!      REAL(KIND(1.D0)),intent(in) :: d(3), lambda, pcoh
!      REAL(KIND(1.D0)) :: res, sigt(2)
!      call NITI_TOT(lambda, pcoh, sigt)
!      res = d(1)*sigt(1) + d(2)*sigt(2)
!      NITI_AVE = res/d(3)
!      end function  NITI_AVE


!---------------------------------------------------------
      subroutine  NITI_MU(lambda, qT, pcoh, mu)
! Depth attenuation coefficient [1/cm]
! lambda = wavelength [A]
! pcoh = coherent scattering factor
! qT = 1/2 of the reflection vector in qTNi units
! Returns in mu:
! Total removal cross-section for Ni (1) and Ti (2)
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: lambda, qT, pcoh
      REAL(KIND(1.D0)),intent(out) :: mu(2)
      REAL(KIND(1.D0)), parameter :: PI4=2.D0*TWOPI
      integer :: i
      REAL(KIND(1.D0)) :: sigt(2), q, A, B
1     format(a,6(1x, G12.4))
      call NITI_TOT(lambda, pcoh, sigt)
      do i=1,2
        q = qT*qTNi
        A = q*q - PI4*rho(i)
        B = TWOPI/lambda*sigt(i)*1.D-8 ! convert 1/cm to 1/A
        ! calculate mu = 2*Im(sqrt(A + iB))
        mu(i) = 1.0D8*sqrt(2.D0*(sqrt(A*A + B*B) - A))
        !write(*,1) 'NITI_MU qT, q, PI4, rho: ', qT, q, PI4, rho(i)
        !write(*,1) 'NITI_MU pcoh, qT, A, B, sigt,mu: ', int(pcoh), A, B, sigt(i), mu(i)
      enddo

      end subroutine  NITI_MU

!---------------------------------------------------------
      REAL(KIND(1.D0)) function  NITI_ATT(qT, d, mval, R0, lambda, pcoh, PR)
! Attenuation factor for incident neutrons.
! qT = 1/2 of the reflection vector in qTNi units
! d = array(3) with thickness in um: (1) all, (2) Ni, (3) Ti
! mval = supermirror m-value
! R0 = reflectivity at qT<1
! lambda = wavelength [A]
! pcoh = coherent scattering factor
! PR probability of reflection
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: qT, d(3), mval, R0, lambda, pcoh, PR
      REAL(KIND(1.D0)) :: mu(2), arg, arg1, res, a1, a2
      REAL(KIND(1.D0)), parameter :: delta = 1.D-2
1     format(a,6(1x, G12.4))
      call NITI_MU(lambda, qT, pcoh, mu)
      arg1 = (d(1)*mu(1)+d(2)*mu(2))*1.D-4 ! d is in um, mu is in 1/cm !!
      !write(*,1) 'NITI_ATT d1, d2, m1, mu2, arg1: ', d, mu, arg1
      if (qT<=1.D0) then
        res = 1.D0 ! no transmission at qT<1
      else
        a1 = exp(-arg1)
        arg = 2.D0*arg1*(qT*qT-1.D0)/(mval*mval-1.D0)
        a2 = exp(-arg)
        res = PR*a2 + (1.D0-PR)*a1
      endif
      !write(*,1) 'NITI_ATT qqT, PR, a1, a2, arg2, ATT: ', qT, PR, a1, a2, arg
      !write(*,1) 'NITI_ATT pcoh, ATT: ', int(pcoh), res
      NITI_ATT = res
      end function  NITI_ATT


!---------------------------------------------------------
      subroutine  NITI_PROB(qT, lambda, d, r, mpar, Rc, PROB)
! Calculates probabilities of each process.
! Input:
! d = array(3) with thickness in um: (1) all, (2) Ni, (3) Ti
! lambda = wavelength [A]
! qT = 1/2 of the reflection vector in qTNi units
! r = reflectivity at qT as given in tables (excludes coherent scattering)
! mpar(2) = reflectivity parameters:
!     mpar(1) = m-value
!     mpar(2) = reflectivity at qT<1
! Output:
! Rc = reflectivity at qT corrected for coherent scattering (if incl_coh=true)
! PROB(4) = array of probabilities:
!   (1) capture
!   (2) incoherent scattering
!   (3) coherent scattering
!   (4) free transmision
! The algorithm ensures that Rc + sum(P) = 1
! If > 1, then P(3) is adjusted so that P(4)=0 and Rc + sum(P) = 1 holds
! This can happen if the tabeled reflectivity gives higher values then is
! the calculated 1 - sum of capture and incoheret scattering.
!---------------------------------------------------------
      REAL(KIND(1.D0)),intent(in) :: d(3), lambda, qT, r, mpar(2)
      REAL(KIND(1.D0)),intent(out) :: Rc, PROB(4)
      integer :: i
      logical :: useNiTi
      REAL(KIND(1.D0)) :: P(4), SS(4), r1, NORM, mval, R0, x
      REAL(KIND(1.D0)) :: ATT, ATT0, DATT
1     format(a,6(1x, G12.4))
      mval = mpar(1)
      R0 = mpar(2)
      useNiTi = (mval>1.3)
      P = 0.D0

      ! divide attenuation probability to individual processes:
      ! SS(1) ... capture
      ! SS(2) ... incoherent scattering
      ! SS(3) ... coherent scattering
      SS = 0.D0
      do i=1,2
        SS(1) = SS(1) + siga(i)*d(i)*lambda
        SS(2) = SS(2) + sigi(i)*d(i)
        if (incl_coh) then
          if (lambda<bge(i)) SS(3) = SS(3) + d(i)*sigc(i)
        endif
      enddo
      NORM = SS(1)+SS(2)+SS(3)
      SS = SS/NORM
      ! Calculate total probabilities
      ! PP(1) ... capture
      ! PP(2) ... incoherent scattering
      ! PP(3) ... coherent scattering
      ! PP(4) ... free transmission
      ! Rc ... reflection
      ! ensure SUM(P) + Rc = 1

      ! divide (1-R0) part according to cross-sections:
      P = SS*(1.D0 - R0)
      r1 = r
      ATT = 1.D0
      if ((qT>1.D0).and.useNiTi) then
        ! correct reflectivity for coherent scattering if required
        if (incl_coh.and.(qT>1.D0).and.(r>0.D0)) then
          ATT0 = NITI_ATT(qT, d, mval, R0, lambda, 0.D0, r/R0)
          ATT = NITI_ATT(qT, d, mval, R0, lambda, 1.D0, r/R0)
          DATT = ATT0-ATT ! drop of reflecivity due to coherent scattering
          ! correction of reflectivity for coherent scattering
          r1 = r*(1.D0 - DATT)
        endif
        ! attenuation factor, 1-ATT = removal below r=R0 line
        ATT = NITI_ATT(qT, d, mval, R0, lambda, 1.D0, r1/R0)
        ! add removal probabilities below r = R0:
        P = P + SS*R0*(1.D0 - ATT)
      endif
      ! Assumed free transmission = rest to 1
      P(4) = 1.D0 - r1 - (P(1)+P(2)+P(3))
      ! If necessary, make adjustments so that r1 + P1 + P2 + P3 + x = 1 for with x >=0
      if (P(4) < 0.D0) then
        ! nothing left for transmision
        ! try to match SUM(P)+r1=1 by lowering scattering probability
        P(4)=0.D0
        x = 1.D0 - r1 - (P(1)+P(2))
        if (x>=0.D0) then
          P(3) = x
        else
        ! if it is not enough, do the same with incoherent scattering
          P(3) = 0.D0
          x = 1.D0 - r1 - P(1)
          if (x>=0.D0) then
            P(2) = x
          else
          ! If it still does not help, reduce capture. This should not happen
          ! unless the tabeled reflectivity is too high ...
            P(2) = 0.D0
            x = 1.D0 - r1
            if (x>=0.D0) then
              P(1) = x
            endif
          endif
        endif
      endif
      !! now we are sure that SUM(P) + r1 = 1
      Rc = r1
      PROB = P
      end subroutine  NITI_PROB

      end  module NITI





