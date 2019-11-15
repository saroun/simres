!//////////////////////////////////////////////////////////////////////
!////  $Id: table_atoms.f,v 1.6 2013/07/12 08:16:51 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2006, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.6 $
!////     $Date: 2013/07/12 08:16:51 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Table with scattering lengths, absortption cross-section etc
!//////////////////////////////////////////////////////////////////////
      MODULE TABLE_ATOMS
      implicit none

      private

      integer,parameter :: MAX_TABLE_ATOMS=100
      TYPE TTAB_ATOMS; SEQUENCE
        character*2 ::  ID               !  atom symbol
        INTEGER ::  Z                    !  proton number
        real(kind(1.D0)) ::  M           !  mass [at. mass units]
        real(kind(1.D0)) ::  BC          !  mean coherent scattering length (bound) [fm]
        real(kind(1.D0)) ::  SIGI        !  incoherent scattering cross-section     [barn]
        real(kind(1.D0)) ::  SIGA        !  absorption cross-section [barn/A]
      end TYPE TTAB_ATOMS

      TYPE(TTAB_ATOMS) :: SC_TABLE(MAX_TABLE_ATOMS)

      interface  GET_ATOM_DATA
        module procedure GET_ATOM_DATA_ID
        module procedure GET_ATOM_DATA_IDX
      end interface

      public GET_ATOM_DATA
      public TABLE_ATOMS_READLN
      public GET_ATOM_INDEX
      public TTAB_ATOMS

      contains


!---------------------------------------------------------
      integer function GET_ATOM_INDEX(ID)
! return index of atom data in SC_TABLE by ID
!---------------------------------------------------------
        character(*),intent(in) :: ID
        integer :: i,ia
        i=0
        ia=0
        do while ((i.lt.MAX_TABLE_ATOMS).and.(ia.eq.0))
          i=i+1
          if (trim(ID).EQ.trim(SC_TABLE(i)%ID)) ia=i
        enddo
        GET_ATOM_INDEX=ia
      end function GET_ATOM_INDEX

!---------------------------------------------------------
      subroutine GET_ATOM_DATA_ID(ID,ATOM)
! return atom data from SC_TABLE by name
!---------------------------------------------------------
        character(*),intent(in) :: ID
        TYPE(TTAB_ATOMS),intent(out) :: ATOM
        call GET_ATOM_DATA_IDX(GET_ATOM_INDEX(ID),ATOM)
      end subroutine GET_ATOM_DATA_ID

!---------------------------------------------------------
      subroutine GET_ATOM_DATA_IDX(IDX,ATOM)
! return atom data from SC_TABLE by index
!---------------------------------------------------------
        INTEGER,intent(in) :: IDX
        TYPE(TTAB_ATOMS),intent(out) :: ATOM
        if (IDX.gt.0) then
          ATOM=SC_TABLE(IDX)
        else
          call GET_EMPTY_DATA(ATOM)
        endif
      end subroutine GET_ATOM_DATA_IDX

!---------------------------------------------------------
      subroutine GET_EMPTY_DATA(ATOM)
! return empty atom data
!---------------------------------------------------------
        TYPE(TTAB_ATOMS),intent(out) :: ATOM
        ATOM%ID='  '
        ATOM%Z=0
        ATOM%BC=0.D0
        ATOM%SIGI=0.D0
        ATOM%SIGA=0.D0
        ATOM%M=0.D0
      end subroutine GET_EMPTY_DATA

!--------------------------------------------------------
      subroutine TABLE_ATOMS_READLN(LINE,error)
! read table row from text line in the format
! Z ID M BC SIGI SIGA
!--------------------------------------------------------
        character(len=*),intent(in) :: LINE
        logical,intent(out) :: error
        integer :: I
        real(kind(1.D0)) :: z
        TYPE(TTAB_ATOMS) :: AT
        logical :: isInteger,isFloat
        character*(256) :: S
        integer :: NA,IS,IL
        error=.false.
        call GET_EMPTY_DATA(AT)
      ! ignore comments and empty lines
        if ((len_trim(LINE).gt.0).and.(LINE(1:1).ne."#")) then
          call READ_STRARRAY(LINE,'|',S,NA)
          if (NA.ge.6) then
            call FINDSTRPAR(S,'|',1,IS,IL)
            if (isInteger(S(IS:IS+IL-1),i)) AT%Z=i
            call FINDSTRPAR(S,'|',2,IS,IL)
            AT%ID=S(IS:IS+IL-1)
            call FINDSTRPAR(S,'|',3,IS,IL)
            if (isFloat(S(IS:IS+IL-1),z)) AT%M=z
            call FINDSTRPAR(S,'|',4,IS,IL)
            if (isFloat(S(IS:IS+IL-1),z)) AT%BC=z
            call FINDSTRPAR(S,'|',5,IS,IL)
            if (isFloat(S(IS:IS+IL-1),z)) AT%SIGI=z
            call FINDSTRPAR(S,'|',6,IS,IL)
            if (isFloat(S(IS:IS+IL-1),z)) AT%SIGA=z
          endif
          error=((len_trim(AT%ID).le.0).or.(AT%Z*AT%BC.eq.0.D0))
          if (.not.error) then
            if (AT%Z.le.MAX_TABLE_ATOMS) SC_TABLE(AT%Z)=AT
          endif
        endif
      end subroutine TABLE_ATOMS_READLN

      END MODULE TABLE_ATOMS
