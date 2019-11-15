!//////////////////////////////////////////////////////////////////////
!////  $Id: messages.f,v 1.2 2016/05/04 14:48:51 saroun Exp $
!////
!////  R E S T R A X - Simulation of neutron three-axis spectrometers
!////
!////     Copyright (C) 1995-2014, All rights reserved
!////     Nuclear Physics Institute AVCR, Rez, Czech Republic
!////     Institut Laue Langevin, Grenoble, France
!////
!////     Written by:  Jan Saroun
!////     $Revision: 1.2 $
!////     $Date: 2016/05/04 14:48:51 $
!//////////////////////////////////////////////////////////////////////
!////
!////  Messages system
!////
!//////////////////////////////////////////////////////////////////////
      MODULE MESSAGES
      USE XMLINFO
      USE IO
      implicit none
      SAVE
      PRIVATE

! dimensions
      integer,PARAMETER :: LEN_MSG=1024  ! maximum message length
      integer,PARAMETER :: LEN_SRC=32    ! default length of input line
      integer,PARAMETER :: MESSAGES_DIM=64  ! maximum number of messages in queue

! message types
      integer,PARAMETER :: m_info=0
      integer,PARAMETER :: m_warn=1
      integer,PARAMETER :: m_error=2
! message priorities
      integer,PARAMETER :: m_low=0
      integer,PARAMETER :: m_high=1

      TYPE TMESSAGE; SEQUENCE
        INTEGER :: TYP
        INTEGER :: PRIORITY
        character(LEN_SRC)   :: SOURCE
        CHARACTER(LEN_MSG) :: MSG
      END TYPE TMESSAGE

! pointer type
      type PMESSAGE
        integer :: IDX  ! instance index
        TYPE(TMESSAGE),pointer :: X
      end type PMESSAGE

! array of messages
      integer :: MESSAGE_NC
      type(PMESSAGE) :: AMSG(1:MESSAGES_DIM)

      PUBLIC ADD_MESSAGE,FREE_MESSAGES,FLUSH_MESSAGES,APPEND_TEXT,INIT_MESSAGES
      PUBLIC m_info,m_warn,m_error,m_low,m_high

      contains

!---------------------------
! PUBLIC METHODS
!---------------------------

!---------------------------------------------------------------
      SUBROUTINE INIT_MESSAGES
! initialize messages
!---------------------------------------------------------------
      integer :: i
      MESSAGE_NC=0
      do i=1,MESSAGES_DIM
        AMSG(i)%IDX=0
      enddo
      end SUBROUTINE INIT_MESSAGES

!---------------------------------------------------------------
      integer function  ADD_MESSAGE(TYP,SRC,MSG,IPR)
! Add message to the queue
! return message ID
!---------------------------------------------------------------
      INTEGER,intent(in) :: TYP,IPR
      CHARACTER*(*),intent(in) :: SRC,MSG
      integer :: i
      i=MSG_CREATE()
      if (i>0) then
        AMSG(i)%X%TYP=TYP
        AMSG(i)%X%PRIORITY=IPR
        AMSG(i)%X%SOURCE=trim(SRC)
        AMSG(i)%X%MSG=trim(MSG)
      endif
      ADD_MESSAGE=i
      ! write(*,*) 'ADD_MESSAGE ',MESSAGE_NC,TYP,' ',trim(SRC),' : ',trim(MSG)
      end function ADD_MESSAGE

!---------------------------------------------------------------
      SUBROUTINE APPEND_TEXT(ID,MSG)
! Append line to the ID-th message with '|' delimiter
!---------------------------------------------------------------
      CHARACTER*(*),intent(in) :: MSG
      integer,intent(in) :: ID
      integer :: i,L,LM
      if (MSG_isValid(ID)) then
        i=ID
        LM=len_trim(AMSG(i)%X%MSG)
        L=len_trim(MSG)
        if (LM==0) then
          AMSG(i)%X%MSG=trim(MSG)
        else if ((L>0).and.(LM+L<LEN_MSG)) then
          AMSG(i)%X%MSG=trim(AMSG(i)%X%MSG)//' | '//trim(MSG)
        endif
        ! write(*,*) 'APPEND_TEXT ',LM+L,' ',trim(MSG)
      endif
      end SUBROUTINE APPEND_TEXT

!---------------------------------------------------------------
      SUBROUTINE FREE_MESSAGES
! dispose messages (deallocate buffer etc.)
!---------------------------------------------------------------
      integer :: i
        do i=1,MESSAGES_DIM
          call MSG_DISPOSE(i)
        enddo
        MESSAGE_NC=0
        ! write(*,*) 'FREE_MESSAGES '
      end SUBROUTINE FREE_MESSAGES

!---------------------------------------------------------------
      SUBROUTINE FLUSH_MESSAGES
! dispose messages (deallocate buffer etc.)
!---------------------------------------------------------------
      integer :: i,L
      ! write(*,*) 'FLUSH_MESSAGES ',MESSAGE_NC
      if (MESSAGE_NC>0) then
        call XML_RSXDUMP(SMES,' ',1)
      endif
      do i=1,MESSAGE_NC
        if (MSG_isValid(i)) then
          L=len_trim(AMSG(i)%X%MSG)
          if (L>0) then
            select case (AMSG(i)%X%TYP)
            case(m_info)
              call MSG_INFO(trim(AMSG(i)%X%MSG),0)
            case(m_warn)
              call MSG_WARN(trim(AMSG(i)%X%MSG),0)
            case(m_error)
              call MSG_ERROR(trim(AMSG(i)%X%SOURCE),trim(AMSG(i)%X%MSG),AMSG(i)%X%PRIORITY,0)
            end select
          endif
        endif
      enddo
      if (MESSAGE_NC>0) then
        call XML_RSXDUMP(SMES,' ',0)
      endif
      call FREE_MESSAGES
      end SUBROUTINE FLUSH_MESSAGES


!---------------------------
! PRIVATE METHODS
!---------------------------

!-------------------------------------------------------------------------
! Create message instance, return instance number
! Memory is allocated on the first free element of AMSG array
!-------------------------------------------------------------------------
      integer function MSG_CREATE()
      integer :: ierr,i
      character(32) :: S
        ierr=1
        i=1
        do while ((i<MESSAGES_DIM).and.(AMSG(i)%IDX>0))
          i=i+1
        enddo
        if (AMSG(i)%IDX.le.0) then
          allocate(AMSG(i)%X,STAT=ierr)
          if(ierr.eq.0) then
            MESSAGE_NC=MESSAGE_NC+1
            AMSG(i)%IDX=i
          endif
        endif
        if (ierr.eq.0) then
          MSG_CREATE=i
        else
          call INT2STR(MESSAGES_DIM,S)
          call MSG_ERROR('MSG_CREATE','Reached maximum number or messages in queue  ('//trim(S)//')',1,1)
          MSG_CREATE=0
        endif
      end function MSG_CREATE

!---------------------------------------------------------
      subroutine MSG_DISPOSE(INST)
!---------------------------------------------------------
        integer,intent(in) :: INST
        if ((INST.gt.0).and.(INST.le.MESSAGES_DIM)) then
          if (associated(AMSG(inst)%X)) then
            DEallocate(AMSG(INST)%X)
            nullify(AMSG(INST)%X)
            MESSAGE_NC=MESSAGE_NC-1
          endif
          AMSG(inst)%IDX=0
        endif
      end subroutine MSG_DISPOSE


!-------------------------------------------------------------
      logical function MSG_isValid(INST)
! check index and association
!-------------------------------------------------------------
      integer,intent(in) :: INST
      MSG_isValid= ((INST.gt.0).and.(INST.le.MESSAGES_DIM).and.associated(AMSG(INST)%X))
      end function MSG_isValid


      end MODULE MESSAGES