!//////////////////////////////////////////////////////////////////////
!////  $Id: components_io.f,v 1.35 2019/08/15 15:02:08 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.35 $
!////     $Date: 2019/08/15 15:02:08 $
!//////////////////////////////////////////////////////////////////////
!////
!////  - Storage of namespaces for classes (param. names, enum. definitions etc.)
!////  - Storage of instrument input values (as a single array)
!////  - Data exchange between the input array and beamline objects
!////  - Text output of the input parameters
!////
!////  NOTE:
!////  Parameters of beamline components are stored twice:
!////  1) as fields of corresponding objects, converted to working units
!////  2) as a single REAL*8 BPAR_ARRAY(:), which contains only input data
!////     in user units(those used at the input/output)
!////
!////  Subroutines in this module permit to access parameters in both the places
!////  and read/write the data from/to each of them
!////  *_INP subroutines should be used as the only way of assigning values to the components.
!////  This ensures correct unit conversions and some syntax checks
!////
!//////////////////////////////////////////////////////////////////////
      MODULE COMPONENTS_IO
      use CONSTANTS
      use CLASSDEF
      use FIELDDATA
      use COMPONENTS
      use VCTABLE

      implicit none


! All input parameters of the instrument are stored as a REAL array.
! Appart of the data in objects, the values in arrays are in physical units for input
! Unit conversion is done in appropriate subroutines below

    ! array must have enough space for 2 beamlines + 1 sample
      integer, parameter, private :: BPAR_DIM=(2*MAX_FIELDS*BEAMLINE_DIM+MAX_FIELDS)*16
      integer :: BPAR_NP   ! number of parameters
      real(KIND(1.D0)) :: BPAR_ARRAY(BPAR_DIM) ! values
! partition table, which permits to address component parameters as a single array values
! BEAMLINE_INDX(i-1)+1 is the index of the first value for the i-th class of the instrument
      integer :: BPAR_INDX(0:2*BEAMLINE_DIM+2)
      INTEGER,private :: BPAR_INDX_INS,BPAR_INDX_SAM,BPAR_INDX_OPT(2)
      integer,private :: BPAR_INDX_COM(BEAMLINE_DIM,2)

! local char array for full component namespaces, shared by _INP procedures
!      character(2*DIM_ID),private :: SGLOB

      contains

! ********************  PUBLIC procedures ************************************************


C---------------------------------------------------------
      SUBROUTINE BEAM_PRN(IO)
C Print the list of beamline
C INPUT
C    IO      .... output unit
C---------------------------------------------------------
      INTEGER,intent(in) :: IO
      integer :: i
      type(TCLASS) :: OBJ
      character(LEN_ID) :: CLSID,OBJID
      character(LEN_NAME) :: OBJNAME,CNAME

1     FORMAT('   ',a16,' ',a16,' ',a)
11    format(' [',a,'] ')
2     FORMAT(a)

      write(IO,2) 'OPTIONS:'
      call GetClassAttrib(PCLASS(POPTOBJ(TROPT)),CLSID,OBJID,OBJNAME)
      write(CNAME,11) trim(CLSID)
      write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      call GetClassAttrib(PCLASS(POPTOBJ(REPOPT)),CLSID,OBJID,OBJNAME)
      write(CNAME,11) trim(CLSID)
      write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      if (INSOBJ%ICLS.gt.0) then
        write(IO,2) 'INTERFACE:'
        call GetClassAttrib(PCLASS(INSOBJ),CLSID,OBJID,OBJNAME)
        write(CNAME,11) trim(CLSID)
        write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      endIF
      if (SAMOBJ%ICLS.gt.0) then
        write(IO,2) 'SPECIMEN:'
        call GetClassAttrib(PCLASS(SAMOBJ),CLSID,OBJID,OBJNAME)
        write(CNAME,11) trim(CLSID)
        write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      endif
      do i=1,BEAMLINE_NC(1)
        if (i.eq.1) write(IO,2) 'PRIMARY:'
        OBJ=PCLASS(GetComByIndex(i,1))
        call GetClassAttrib(OBJ,CLSID,OBJID,OBJNAME)
        write(CNAME,11) trim(CLSID)
        write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      enddo
      do i=1,BEAMLINE_NC(2)
        if (i.eq.1) write(IO,2) 'SECONDARY:'
        OBJ=PCLASS(GetComByIndex(i,2))
        call GetClassAttrib(OBJ,CLSID,OBJID,OBJNAME)
        write(CNAME,11) trim(CLSID)
        write(IO,1) trim(adjustl(OBJID)),trim(adjustl(CNAME)),trim(OBJNAME)
      enddo
      END SUBROUTINE BEAM_PRN

C-------------------------------------------------------------
      SUBROUTINE BEAM_OUT_APPEND(OBJ,IDX)
C add next CLASS data to BPAR for output
C--------------------------------------------------------------
      type(TCLASS) :: OBJ
      integer,intent(in) :: IDX
      integer :: NA,NP
      TYPE(TFIELD) :: ARG
      TYPE(PCLSDEF) :: CDEF
      TYPE(TFIELDDEF),pointer :: FDEF
      integer :: i,j,k,m,NCLS,PCLS(1:10)
      type(PCLSDEF) :: CLS1
      character(LEN_ID+1) :: prefix
      logical :: dbg
      dbg=.false.
1     format('BEAM_OUT: incorrect No. of parameters in ',a,' NA=',I5,'NC=',I5)
2     format(a,6(1x,G12.5))
      NP=CLASSES_LENGTH(OBJ%ICLS)
      BPAR_INDX(IDX)=BPAR_INDX(IDX-1)+NP
      BPAR_NP=BPAR_INDX(IDX)
      call GetClassParents(OBJ%ICLS,PCLS,NCLS)
      j=NCLS
      NA=0
      do while((j.gt.0))
        call GetClassDef(PCLS(j),CDEF)
        if (associated(CDEF%C)) then
          do i=1,CDEF%C%NF
            FDEF=>CDEF%C%FIELDS(i)%F
            ARG%DEF=>FDEF
            select case (trim(FDEF%TID))
        ! class field
            case(ID_FIELD_CLASS)
              call GetClassDef(trim(FDEF%CID),CLS1)
              if (associated(CLS1%C).and.(CLS1%C%NF.gt.0)) then
                prefix=trim(FDEF%ID)//'.'
                do m=1,CLS1%C%NF
                  if (associated(CLS1%C%FIELDS(m)%F)) then
                    ARG%DEF=>CLS1%C%FIELDS(m)%F
                    if (ALLOC_FIELD(ARG,ARG%DEF%DIM)) then
                      call PCLASS_OUT(OBJ,trim(prefix)//trim(ARG%DEF%ID),ARG,k)
                      select case(trim(ARG%DEF%TID))
                      case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
                        call FIELD2ARRAY(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                      case(ID_FIELD_TABLE)
                        call FIELD2ARRAY(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                        !write(*,*) 'BEAM_OUT_APPEND 1 TABLE2ARRAY k=',k
                      case default
                        k=1
                      end select
                      NA=NA+k
                    endif
                  endif
                enddo
              endif
            ARG%DEF=>FDEF
        ! standard field
            case default
              if (ALLOC_FIELD(ARG,FDEF%DIM)) then
                call PCLASS_OUT(OBJ,FDEF%ID,ARG,k)
                select case(trim(FDEF%TID))
                case(ID_FIELD_FLOAT,ID_FIELD_INT,ID_FIELD_ENUM)
                  call FIELD2ARRAY(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                case(ID_FIELD_TABLE)
                  call FIELD2ARRAY(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                  !write(*,*) 'BEAM_OUT_APPEND 2 TABLE2ARRAY k=',k
                case default
                  k=1
                end select
                NA=NA+k
              endIF
            end select
            call DISPOSE_FIELD(ARG)
          enddo
        endif
        j=j-1
      enddo
      if (NA.ne.NP) then
        write(*,1) trim(OBJ%DEF%ID)//':'//trim(OBJ%ID),NA,NP
      endif
      ! write(*,2) 'BEAM_OUT_APPEND '//trim(OBJ%DEF%ID),IDX,NA,BPAR_NP
      end SUBROUTINE BEAM_OUT_APPEND

C-------------------------------------------------------------
      SUBROUTINE BEAM_INP_APPEND(OBJ,IDX)
C read next CLASS data from BPAR to a component instance
C--------------------------------------------------------------
      type(TCLASS) :: OBJ
      integer,intent(in) :: IDX
      integer :: NA,NP
      TYPE(TFIELD) :: ARG
      TYPE(PCLSDEF) :: CDEF
      TYPE(TFIELDDEF),pointer :: FDEF
      type(PCLSDEF) :: CLS1
      character(LEN_ID+1) :: prefix
      integer :: i,j,k,m,NCLS,PCLS(1:10)
1     format('BEAM_INP: incorrect No. of parameters in ',a,' NA=',I5,'NC=',I5)
2     format(a,6(1x,G12.5))
      NP=CLASSES_LENGTH(OBJ%ICLS)
      call GetClassParents(OBJ%ICLS,PCLS,NCLS)
      j=NCLS
      NA=0
      do while((j.gt.0))
        call GetClassDef(PCLS(j),CDEF)
        if (associated(CDEF%C)) then
          do i=1,CDEF%C%NF
            FDEF=>CDEF%C%FIELDS(i)%F
            ARG%DEF=>FDEF
            select case (trim(FDEF%TID))
        ! class field
            case(ID_FIELD_CLASS)
              call GetClassDef(trim(FDEF%CID),CLS1)
              if (associated(CLS1%C).and.(CLS1%C%NF.gt.0)) then
                prefix=trim(FDEF%ID)//'.'
                do m=1,CLS1%C%NF
                  if (associated(CLS1%C%FIELDS(m)%F)) then
                    ARG%DEF=>CLS1%C%FIELDS(m)%F
                    if (ALLOC_FIELD(ARG,ARG%DEF%DIM)) then
                      call PCLASS_INP(OBJ,trim(prefix)//trim(ARG%DEF%ID),ARG,k)
                      call ARRAY2FIELD(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                      NA=NA+k
                    endif
                  endif
                enddo
              endif
        ! standard field
            case default
              if (ALLOC_FIELD(ARG,FDEF%DIM)) then
                call PCLASS_INP(OBJ,FDEF%ID,ARG,k)
                call ARRAY2FIELD(ARG,0,BPAR_ARRAY(BPAR_INDX(IDX-1)+NA+1),BPAR_DIM,k)
                NA=NA+k
              endif
            end select
            call DISPOSE_FIELD(ARG)
          enddo
        endif
        j=j-1
      enddo
      if (NA.ne.NP) then
        write(*,1) trim(OBJ%DEF%ID)//':'//trim(OBJ%ID),NA,NP
      endif
      ! write(*,2) 'BEAM_INP_APPEND '//trim(OBJ%DEF%ID),IDX,NA,BPAR_NP
      end SUBROUTINE BEAM_INP_APPEND

C-------------------------------------------------------------
      SUBROUTINE BEAM_OUT
C write data from the beamline components to BPAR_ARRAY
C BPAR_INDX and BPAR_NP are updated here
C--------------------------------------------------------------
      integer :: i_inp
      type(TCLASS) :: OBJ
      integer :: i,ib
! clear BPAR array
      BPAR_INDX(0)=0
      BPAR_NP=0
      i_inp=0
      BPAR_INDX_INS=0
      BPAR_INDX_SAM=0
      BPAR_INDX_COM=0
      BPAR_INDX_OPT=0
! OPTIONS
        i_inp=i_inp+1
        BPAR_INDX_OPT(OCLS_TRACING)=BPAR_NP+1
        OBJ=PCLASS(POPTOBJ(TROPT))
        call BEAM_OUT_APPEND(OBJ,i_inp)
        i_inp=i_inp+1
        BPAR_INDX_OPT(OCLS_REPORTS)=BPAR_NP+1
        OBJ=PCLASS(POPTOBJ(REPOPT))
        call BEAM_OUT_APPEND(OBJ,i_inp)
! INTERFACE
      if (INSOBJ%ICLS.gt.0) then
        i_inp=i_inp+1
        BPAR_INDX_INS=BPAR_NP+1
        OBJ=PCLASS(INSOBJ)
        call BEAM_OUT_APPEND(OBJ,i_inp)
      endif
! SAMPLE
      if (SAMOBJ%ICLS.gt.0) then
        i_inp=i_inp+1
        BPAR_INDX_SAM=BPAR_NP+1
        OBJ=PCLASS(SAMOBJ)
        call BEAM_OUT_APPEND(OBJ,i_inp)
      endif
! BEAMLINE
      DO ib=1,2
        do i=1,BEAMLINE_NC(ib)
          i_inp=i_inp+1
          BPAR_INDX_COM(i,ib)=BPAR_NP+1
          OBJ=PCLASS(getComByIndex(i,ib))
          call BEAM_OUT_APPEND(OBJ,i_inp)
        enddo
      enddo
      ! write(*,*) 'BEAM_OUT NPAR=',BPAR_NP,i_inp,BPAR_DIM
      end SUBROUTINE BEAM_OUT

C-------------------------------------------------------------
      SUBROUTINE BEAM_INP
C write data from BPAR_ARRAY to the instrument
C--------------------------------------------------------------
      integer :: i_inp
      type(TCLASS) :: OBJ
      integer :: i,ib
      i_inp=0
! OPTIONS
        i_inp=i_inp+1
        OBJ=PCLASS(POPTOBJ(TROPT))
        call BEAM_INP_APPEND(OBJ,i_inp)
        i_inp=i_inp+1
        OBJ=PCLASS(POPTOBJ(REPOPT))
        call BEAM_INP_APPEND(OBJ,i_inp)
! INTERFACE
      if (INSOBJ%ICLS.gt.0) then
        i_inp=i_inp+1
        OBJ=PCLASS(INSOBJ)
        call BEAM_INP_APPEND(OBJ,i_inp)
      endif
! SAMPLE
      if (SAMOBJ%ICLS.gt.0) then
        i_inp=i_inp+1
        OBJ=PCLASS(SAMOBJ)
        call BEAM_INP_APPEND(OBJ,i_inp)
      endif
! BEAMLINE
      DO ib=1,2
        do i=1,BEAMLINE_NC(ib)
          i_inp=i_inp+1
          OBJ=PCLASS(getComByIndex(i,ib))
          call BEAM_INP_APPEND(OBJ,i)
        enddo
      enddo
      ! write(*,*) 'BEAM_INP i_inp=',BPAR_NP,i_inp,BPAR_DIM
      end SUBROUTINE BEAM_INP

      end MODULE COMPONENTS_IO
