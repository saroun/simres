C// Ray-tracing control data
C// $Id: tracing.inc,v 1.1 2007/10/08 08:35:37 saroun Exp $

c global neutron structures
      STRUCTURE  /NEUTRON/
          REAL*8 R(3),K(3),S(3),T,P,K0 ! r,k,spin,time,probability,|k|
      END STRUCTURE
      RECORD /NEUTRON/ NEUI,NEUF,NEUI1,NEUF1
      COMMON /NEUIF/ NEUI,NEUF,NEUI1,NEUF1


c setting of current MC run
      INTEGER*4 NEV     ! target number of events
      INTEGER*4 TRACE_UP ! =1 for up-stream simulation between source and sample
      COMMON /MCSETTING/ NEV,TRACE_UP
