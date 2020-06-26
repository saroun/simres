! xmlparse.f90 - Simple, limited XML parser in Fortran
!
! $Id: xmlparse.f,v 1.24 2019/08/16 17:16:27 saroun Exp $
! Adapted for SIMRES from xmlparse.f90,v 1.13 2007/04/03 19:23:20 by Arjen Markus
! no more compatible with other code from the original project
! Main changes:
!   format changed to fixed, 128 columns
!   changed type XML_PARSE
!      - contains tag name
!      - contains info on the type of element (start, end tag, data)
!   rewritten xml_get
!      - only one line of data per call
!      - skips all comments (including multilines)
!   added sml_get_attrib, which processes attributes separately

! -------------------------------------------------------------------
! Arjen Markus
!
! General information:
! The module reads XML files by:
! - Identifying the tag and all attributes and data belonging
!   to the tag.
! - Returning to the calling subprogram to let it take care of
!   the tag, attributes and data.
! - If the tag is actually an ending tag, then this is flagged
!   too.
! - Handling all the data is left to the calling subprogram,
!   the module merely facilitates in the parsing.
!
! Note:
! The module in its current version has a number of limitations:
! - It does not handle escape sequences (like &gt. to signify
!   a ">" sign)
! - It does not handle tags with attributes that are spread
!   over more than one line
! - The maximum length of a line is 1000 characters
! - It may report too many lines of data (empty lines)
! - No DOM support nor support for an object tree
! - It is probably not very robust in detecting malformed XML files
!
! Some questions:
! - What to do with leading blanks?
!
! Update - several ideas:
! - Introduce at least two options (via xml_options):
!   - ignore_whitespace  - remove leading blanks and leading and trailing
!                          empty lines from the PCDATA
!   - no_data_truncation - consider truncation of data (more
!                          attributes or lines of character data than
!                          can be stored) a read error
! - Introduce convenience functions and subroutines:
!   - xml_ok()           - all is well, reading can continue
!   - xml_data_trunc()   - was there truncation of the data?
!   - xml_find_attrib()  - find an attribute by name
!
! Further ideas:
!   - simple checking via a table: parent, tag, id, min, max
!-------------------------------------------------------------------------
      module XMLPARSE
      USE XMLINFO
      implicit none

      integer, parameter :: XML_BUFFER_LENGTH = 1024
      integer, parameter :: XML_TAG_LENGTH = 256
      integer, parameter,PRIVATE :: TABL = 2
   !
   ! Define the data type that holds the parser information
   !
      type XML_PARSE
      integer          :: lun                ! LU-number of the XML-file
      integer          :: level              ! Indentation level (output)
      integer          :: lineno             ! Line in file
      logical          :: ignore_whitespace  ! Ignore leading blanks etc.
      logical          :: no_data_truncation ! Do not allow data truncation
      logical          :: too_many_attribs   ! More attributes than could be stored?
      logical          :: eof                ! End of file?
      logical          :: error              ! Invalid XML file or other error?
      logical          :: comment_tag        ! continuation of a comment tag
      logical          :: hasStart           ! line contains a starting tag
      logical          :: hasEnd             ! line contains an end tag
      logical          :: hasData            ! line contains data
      LOGICAL          :: mustRead           ! file open for read (T) or write (F)
      character(len=XML_TAG_LENGTH) :: tag   ! current tag name
      character(len=XML_BUFFER_LENGTH) :: tagpath   ! current tag path
      character(len=256) :: fname            ! current file name
      end type XML_PARSE

  ! added by JS: info about current path in the document tree
      character(len=10*XML_TAG_LENGTH), private :: tagpath   ! colon separated path
      character(len=XML_TAG_LENGTH) :: parenttag             ! parent tag of the current tag

   !
   ! Global options
   !
      integer, parameter    :: XML_STDOUT       = -1
      integer, private      :: report_lun_      = XML_STDOUT
      logical, private      :: report_errors_   = .true.
      logical, private      :: report_details_  = .false.

   !
   ! Global data (the ampersand must come first)
   !
      character(len=10), dimension(2,3), save, private :: entities =
     & reshape( (/ '&    ', '&amp;',
     &             '>    ', '&gt; ',
     &             '<    ', '&lt; ' /), (/2,3/) )


   ! global constants - xml element type
      integer, parameter    :: exml_starttag=1
      integer, parameter    :: exml_endtag=2
      integer, parameter    :: exml_data=3
      integer, parameter    :: exml_comment=4

   ! Auxiliary routines - private

      private               :: xml_put_open_tag_
      private               :: xml_put_element_
      private               :: xml_put_close_tag_
      private               :: xml_replace_entities_
   !
   ! Interfaces to reporting routines
   !
      private               :: xml_report_details_int_
      private               :: xml_report_details_string_
      private               :: xml_report_errors_int_
      private               :: xml_report_errors_string_

      interface xml_report_details
        module procedure xml_report_details_int_
        module procedure xml_report_details_string_
      end interface
      interface xml_report_errors
        module procedure xml_report_errors_int_
        module procedure xml_report_errors_string_
        module procedure xml_report_errors_extern_
      end interface

      type(XML_PARSE) :: locinfo
      private locinfo
      contains

!-----------------------------------------------------------------------
! xml_report_details_int_ --
!    Routine to write a text with an integer value
! Arguments:
!    text        Text to be written
!    int         Integer value to be added
!
!-----------------------------------------------------------------------
      subroutine xml_report_details_int_( text, int )
      character(len=*), intent(in)     :: text
      integer,          intent(in)     :: int

      if ( report_details_ ) then
      if ( report_lun_ .eq. XML_STDOUT ) then
         write(*,*) trim(text), int
      else
         write(report_lun_,*) trim(text), int
      endif
      endif
      end subroutine xml_report_details_int_

!-----------------------------------------------------------------------
! xml_report_details_string_ --
!    Routine to write a text with a string value
! Arguments:
!    text        Text to be written
!    string      String to be added
!
!-----------------------------------------------------------------------
      subroutine xml_report_details_string_( text, string )
      character(len=*), intent(in)     :: text
      character(len=*), intent(in)     :: string
      if ( report_details_ ) then
        if ( report_lun_ .eq. XML_STDOUT ) then
          write(*,*) trim(text), ' ', trim(string)
        else
          write(report_lun_,*) trim(text), ' ', trim(string)
        endif
      endif
      end subroutine xml_report_details_string_


!-----------------------------------------------------------------------
! xml_report_errors_string_ --
!    Routine to write an error message text with an integer value
! Arguments:
!    text        Text to be written
!    int         Integer value to be added
!    lineno      Line number in the file
!
!-----------------------------------------------------------------------
      subroutine xml_report_errors_int_( text, int, lineno )
      character(len=*),  intent(in)     :: text
      integer,           intent(in)     :: int
      integer, optional, intent(in)     :: lineno

      if ( report_errors_ .or. report_details_ ) then
      if ( report_lun_ .eq. XML_STDOUT ) then
         write(*,*) trim(text), int
         if ( present(lineno) ) then
            write(*,*) '   At or near line', lineno
         endif
      else
         write(report_lun_,*) trim(text), int
         if ( present(lineno) ) then
            write(report_lun_,*) '   At or near line', lineno
         endif
      endif
      endif
      end subroutine xml_report_errors_int_

!-----------------------------------------------------------------------
! xml_report_errors_string_ --
!    Routine to write an error message text with a string value
! Arguments:
!    text        Text to be written
!    string      String to be added
!    lineno      Line number in the file
!
!-----------------------------------------------------------------------
      subroutine xml_report_errors_string_( text, string, lineno )
      character(len=*),  intent(in)     :: text
      character(len=*),  intent(in)     :: string
      integer, optional, intent(in)     :: lineno

      if ( report_errors_ .or. report_details_ ) then
      if ( report_lun_ .eq. XML_STDOUT ) then
         write(*,*) trim(text), ' ', trim(string)
         if ( present(lineno) ) then
            write(*,*) '   At or near line', lineno
         endif
      else
         write(report_lun_,*) trim(text), ' ', trim(string)
         if ( present(lineno) ) then
            write(report_lun_,*) '   At or near line', lineno
         endif
      endif
      endif
      end subroutine xml_report_errors_string_

!-----------------------------------------------------------------------
! xml_report_errors_extern_ --
!    Routine to write an error message text with a string value
! Arguments:
!    info        Structure holding information on the XML-file
!    text        Text to be written
! Note:
!    This routine is meant for use by routines outside
!    this module
!
!-----------------------------------------------------------------------
      subroutine xml_report_errors_extern_( info, text )
      type(XML_PARSE),   intent(in)     :: info
      character(len=*),  intent(in)     :: text

1     format(a, ",line=",G10.4,", tag=",a,", tagpath=",a)
      if ( report_lun_ .eq. XML_STDOUT ) then
      write(*,1) trim(text), info%lineno, trim(info%tag), trim(info%tagpath)
      else
      write(report_lun_,1) trim(text), info%lineno, trim(info%tag), trim(info%tagpath)
      endif
      end subroutine xml_report_errors_extern_

!-----------------------------------------------------------------------
      subroutine clearInfo(info)
!-----------------------------------------------------------------------
      type(XML_PARSE),  intent(out)    :: info
      info%lun = 6
      info%ignore_whitespace  = .false.
      info%no_data_truncation = .false.
      info%too_many_attribs   = .false.
      info%eof                = .false.
      info%error              = .false.
      info%level              = -1
      info%lineno             =  0
      info%comment_tag = .false.
      info%hasStart = .false.
      info%hasEnd = .false.
      info%hasData = .false.
      info%fname='STDIN'
      info%tag=''
      end subroutine clearInfo

!-----------------------------------------------------------------------
! xml_open --
!    Routine to open an XML file for reading or writing
! Arguments:
!    info        Structure holding information on the XML-file
!    fname       Name of the file
!    mustread    The file will be read (.true.) or written (.false.)
!-----------------------------------------------------------------------
      subroutine xml_open( info, fname, mustread )
      character(len=*), intent(in)     :: fname
      logical,          intent(in)     :: mustread
      type(XML_PARSE),  intent(out)    :: info

      integer                          :: i
      logical                          :: opend
      logical                          :: exists
      call clearInfo(info)
      info%mustRead=mustread
      if (trim(fname).ne.'STDIN') then
        info%lun = 10
        do i = 10,99
          inquire( unit = i, opened = opend )
          ! write(*,*) 'xml_open ',i, opend
          if ( .not. opend ) then
            info%lun = i
            inquire( file = fname, exist = exists )
            if ( .not. exists .and. mustread ) then
              call xml_report_errors( 'XML_OPEN: file does not exist:', trim(fname))
              info%lun   = -1
              info%error = .true.
            else
              info%fname=trim(adjustl(fname))
              open( unit = info%lun, file = info%fname, ERR=99 )
              call xml_report_details( 'XML_OPEN: opened file ', info%fname )
              call xml_report_details( 'at LU-number: ', info%lun )
            endif
            EXIT
          endif
        enddo
      else
        if (.not.mustread) call XML_RSXDUMP(info%lun,' ',1)
      endif
    ! write xml header line, not on console
      if ( (.not. info%error) .and. (.not. mustread)) then
        if (trim(fname).ne.'STDIN') write( info%lun, '(a)' ) '<?xml version="1.0" ?>'
      endif
      return

99    call xml_report_errors( 'XML_OPEN: cannot open file ', trim(fname))
      info%lun   = -1
      info%error = .true.
      return
      end subroutine xml_open

!-----------------------------------------------------------------------
! xml_close --
!    Routine to close an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!
!-----------------------------------------------------------------------
      subroutine xml_close(info)
      type(XML_PARSE),  intent(inout)    :: info
      integer :: ierr
      if (trim(info%fname).ne.'STDIN') then
        close( info%lun, iostat=ierr )
   ! Only clean up the LU-number, so that the calling program
   ! can examine the last condition
        if (ierr.le.0) then
          call xml_report_details( 'XML_CLOSE: Closing file with LU-number ', info%lun )
        endif
        info%lun              = -1
      else
        if (.not.info%mustread) call XML_RSXDUMP(info%lun,' ',0)
      endif
      end subroutine xml_close


!-----------------------------------------------------------------------
! Get attributes from a string
! attributes must be a space separated sequence of NAME="VALUE" strings
!-----------------------------------------------------------------------
      subroutine xml_get_attrib(info,line,attribs, no_attribs)
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in)                  :: line
      character(len=*), intent(out), dimension(:,:) :: attribs
      integer,          intent(out)                 :: no_attribs
      integer :: idxat,keq,kfirst,ksecond,kend,kname
      idxat=0
      keq=index(line,'=')
      kend=0
      no_attribs=0
      do while ((keq.gt.kend+1).and.(idxat .lt. size(attribs,2)))
        kname    = kend+1
        kfirst   = index(line(keq+1:), '"' )+keq
        ksecond  = index(line(kfirst+1:), '"' )+kfirst
        if ((kname.ge.keq).or.(kfirst.le.keq).or.(ksecond.le.kfirst)) then
          call xml_report_errors( 'XML_GET - malformed attribute-value pair: ',trim(line(kend+1:)),info%lineno)
          info%error = .true.
          return
        endif
        idxat = idxat + 1
        attribs(1,idxat) = adjustl(line(kname:keq-1))
        if (len_trim(attribs(1,idxat)).le.0) then
          call xml_report_errors( 'XML_GET - no name for an attribute: ',trim(line(kend+1:)),info%lineno)
          info%error = .true.
          return
        endif
        if (kfirst.eq.ksecond-1) then
          attribs(2,idxat) =' '
        else
          attribs(2,idxat) = adjustl(line(kfirst+1:ksecond-1))
        endif
        no_attribs=idxat
        kend=ksecond
        keq=index(line(kend+1:),'=')+kend
        if ((keq.gt.kend+1).and.(idxat.ge.size(attribs,2))) then
          call xml_report_errors( 'XML_GET - more attributes than could be stored: ',trim(line(kend+1:)),info%lineno)
          info%too_many_attribs = .true.
          return
        endif
      endDO
      end subroutine xml_get_attrib

!-----------------------------------------------------------------------
! read a single line from given unit
! check for EOF and skip spaces
!-----------------------------------------------------------------------
      subroutine xml_readline(info,line)
      type(XML_PARSE),  intent(inout)  :: info
      character(len=*), intent(out):: line
      character(len=XML_BUFFER_LENGTH) :: nextline
      integer :: ierr,N

      nextline=' '
      if ( info%lun .lt. 0 ) then
        call xml_report_details( 'XML_READLINE on closed file ', ' ' )
        info%error=.true.
        return
      endif
      info%lineno = info%lineno + 1
      if (info%lun.ne.6) then
        read( info%lun, '(a)', iostat = ierr) nextline
      else
        call RSXREAD(nextline)
        ierr=0
        call SBUFF_LCOUNT(N)
        if (N.le.0) ierr=-1
      endif
      if (ierr.le.0) then
        line = trim(adjustl(nextline))
        info%error=.false.
      else
        line=' '
        info%error=.true.
      endif
      info%eof=(ierr.lt.0)

!      call xml_report_details( 'XML_READLINE: ', trim(nextline))
!      READ(*,*)

      if (info%error) then
        call xml_report_errors( 'XML_READLINE - cannot read from file ',info%lun)
      endif
      end subroutine xml_readline


!-----------------------------------------------------------------------
! xml_get --
!    Routine to get the next bit of information from an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    xmldata     Lines of character data found
!-----------------------------------------------------------------------
      subroutine xml_get(info,line,attribs,no_attribs,xmldata)
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in)                  :: line
      character(len=*), intent(out), dimension(:,:) :: attribs
      integer,          intent(out)                 :: no_attribs
      character(len=*), intent(out)                 :: xmldata

      integer  :: atstart,bracket_len
      integer  :: kend
      integer  :: j,j2,k,ipos
      integer  :: element
      character(len=XML_BUFFER_LENGTH) :: nextline
      character(len=XML_BUFFER_LENGTH) :: atstring

   ! Initialise the output
      no_attribs = 0
      info%error=.false.
      info%too_many_attribs = .false.
      info%hasData  = .false.
      info%hasStart = .false.
      info%hasEnd   = .false.

      nextline=trim(line)

    ! handle xml header line
    ! Assume that <?xml ... ?> appears on a single line!
      if ((.not.info%eof).and.(.not.info%error).and.(info%lineno.eq.1)) then
        k = index(nextline, '<?xml' )
        if (k.ge. 1 ) then
          kend = index( nextline, '?>' )
          if ( kend .le. 0 ) then
            call xml_report_errors( 'XML_GET - error reading file with LU-number: ', info%lun )
            call xml_report_errors( 'Line starting with "<?xml" should end with "?>" in file ', trim(info%fname) )
            info%error = .true.
          endif
          return
        else
          call xml_report_errors( 'Warning: no XML header in file', trim(info%fname) )
        endif
      endif

!      call xml_report_details( 'XML_GET - start: ', nextline )
!      READ(*,*)

    ! possible continuation of a comment from previous line
      if (info%comment_tag) then
        kend=index(nextline,'-->')
        if (kend.gt.0) then
          nextline=nextline(kend+3:)
          info%comment_tag=.false.
        else
          return
        endif
      endif

      do while (len_trim(nextline).gt.0)
    ! determine the type of next xml element
        if (nextline(1:4).eq.'<!--') then
          element=exml_comment
        else if (nextline(1:2).eq.'</') then
          element=exml_endtag
        else if (nextline(1:1).eq.'<') then
          element=exml_starttag
        else
          element=exml_data
        endif

!      write(*,*) 'nestline :',trim(nextline),info%lineno,kend
!      call xml_report_details( 'XML_GET - element found:', element)

    ! get the end position of current element
        info%comment_tag=.false.
        select case(element)
        case(exml_data)
          kend=index(nextline,'<')-1
          if (kend.le.-1) kend=len_trim(nextline)
        case(exml_comment)
          kend=index(nextline,'-->')+2
          if (kend.le.2) then       ! end of the comment not found on the current line
            kend=len_trim(nextline)
            info%comment_tag=.true. ! next line starts as a comment
            return
          endif
        case default
          kend=index(nextline,'>')
          if (kend.lt.3) then
            info%error = .true.
            call xml_report_errors( 'XML_GET - end of tag not found ',info%lineno)
            return
          endif
        end select

      ! handle individual element types
        select case(element)
        case(exml_data)
          if ( info%ignore_whitespace ) then
            xmldata=adjustl(nextline(1:kend))
            info%hasData=(len_trim(xmldata).gt.0)
          else
            xmldata=nextline(1:kend)
            info%hasData=(kend.gt.0)
          endif
        ! replace entities
          do k = 1,size(entities,2)
            j=1
            ipos   = index( xmldata(j:), trim(entities(2,k)) )
            do while (ipos.gt.0)
              ipos   = index( xmldata(j:), trim(entities(2,k)) )
              if ( ipos .gt. 0 ) then
                j     = j + ipos - 1
                j2    = j + len_trim(entities(2,k))
                xmldata(j:) = trim(entities(1,k)) // xmldata(j2:)
                j     = j2
              endif
              ipos = index( xmldata(j:), trim(entities(2,k)) )
            enddo
          enddo
        case(exml_endtag)
          if (info%hasStart.and.(nextline(3:kend-1).ne.trim(info%tag))) then
            info%error = .true.
            call xml_report_errors( 'XML_GET - unpaired end tag ',info%lineno)
            return
          endif
          info%tag=nextline(3:kend-1)
          info%hasEnd=.true.
        case(exml_starttag)
      ! check format: no end tag before start tag
          if (info%hasEnd) then
            info%error = .true.
            call xml_report_errors( 'XML_GET - start tag after end tag ',info%lineno)
            return
          endif
          info%hasStart=.true.
          bracket_len=1
          if (nextline(kend-1:kend).eq.'/>') then ! combined start/end tag, like <BR/>
            info%hasEnd=.true.
            bracket_len=2
          endif
        ! start of the attribute section
          atstart=index(nextline(1:kend),' ')+1
        ! no attributes
          if (atstart.lt.4) then
            info%tag=nextline(2:kend-bracket_len)
            atstring=' '
        ! with attributes
          else
            info%tag=nextline(2:atstart-2)
            atstring=nextline(atstart:kend-bracket_len)
          endif
          if ((index(atstring,'=').le.0).or.(index(atstring,'"').le.0)) then
            atstring=' '
          endif
          call xml_get_attrib(info,atstring,attribs, no_attribs)
        end select
        nextline=adjustl(nextline(kend+1:))
      enddo

      if (report_details_) then
        if (.not.info%error) then
          if (info%hasStart.and.info%hasEnd) then
            call xml_report_details( 'XML_GET ELEM '//trim(info%tag)//' data['//trim(xmldata)//'] attributes: ',no_attribs)
          else if (info%hasStart) then
            call xml_report_details( 'XML_GET START '//trim(info%tag)//', attributes: ', no_attribs)
          else if (info%hasEnd) then
            call xml_report_details( 'XML_GET END ',trim(info%tag))
          else if (info%hasData) then
            call xml_report_details( 'XML_GET DATA ',trim(xmldata))
          endif
        endif
      endif

      end subroutine xml_get

!-----------------------------------------------------------------------
! xml_put --
!    Routine to write a tag with the associated data to an XML file
! Arguments:
!    info        Structure holding information on the XML-file
!    tag         Tag that was encountered
!    endtag      Whether the end of the element was encountered
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    xmldata        Lines of character data found
!    type        Type of action:
!                open - just the opening tag with attributes
!                elem - complete element
!                close - just the closing tag
!
!-----------------------------------------------------------------------
      subroutine xml_put(info, attribs, no_attribs,xmldata,TASK)

      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in), dimension(:,:)  :: attribs
      integer,          intent(in)                  :: no_attribs
      character(len=*), intent(in)                  :: xmldata
      character(len=*)                              :: TASK
      character(len=300), parameter :: indent = ' '
      select case(TASK)
      case('open')
      ! update tagpath and parenttag
        call update_tagpath(info,'open')
        call xml_put_open_tag_(info,attribs, no_attribs)
      case('elem')
        call xml_put_element_(info,attribs, no_attribs,xmldata)
      case('data')
        call xml_put_data_(info,xmldata)
      case('close')
        call xml_put_close_tag_(info)
      ! update tagpath and parenttag
        call update_tagpath(info,'close')
      end select

      end subroutine xml_put

!-----------------------------------------------------------------------
! xml_put_open_tag_ --
!    Routine to write the opening tag with the attributes
! Arguments:
!    info        Structure holding information on the XML-file
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!
!-----------------------------------------------------------------------
      subroutine xml_put_open_tag_(info,attribs,no_attribs)
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in), dimension(:,:)  :: attribs
      integer,          intent(in)                  :: no_attribs
      integer         :: i
      character(len=300), parameter :: indent = ' '

      write( info%lun, '(3a)', advance = 'no' ) indent(1:TABL*info%level), '<', trim(adjustl(info%tag))
      do i=1,no_attribs
        if (len_trim(attribs(1,i)).gt.0) then
          write( info%lun, '(5a)', advance = 'no' ) ' ',trim(attribs(1,i)),'="', trim(attribs(2,i)),'"'
        endif
      enddo
      write( info%lun, '(a)' ) '>'
      info%level = info%level + 1
      end subroutine xml_put_open_tag_

!-----------------------------------------------------------------------
! xml_put_element_ --
!    Routine to write the complete element
! Arguments:
!    info        Structure holding information on the XML-file
!    attribs     List of attribute-value pairs
!    no_attribs  Number of pairs in the list
!    xmldata        Lines of character data found
!
!-----------------------------------------------------------------------
      subroutine xml_put_element_(info,attribs, no_attribs,xmldata)

      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in), dimension(:,:)  :: attribs
      integer,          intent(in)                  :: no_attribs
      character(len=*), intent(in)                  :: xmldata

      logical          :: logic
      integer          :: i, ii

      character(len=300), parameter :: indent = ' '
      logic = .true.
      do ii = 1,no_attribs
          logic = logic .and. (len_trim(attribs(1,ii)).le.0)
      enddo
      logic = logic .and. (len_trim(xmldata).le.0)

      if ( logic ) then
        return
      else
        write( info%lun, '(3a)', advance = 'no' ) indent(1:TABL*info%level), '<', trim(adjustl(info%tag))
        do i = 1,no_attribs
          if (len_trim(attribs(1,i)).gt.0) then
            write( info%lun, '(5a)', advance = 'no' ) ' ',trim(attribs(1,i)),'="', trim(attribs(2,i)),'"'
          endif
        enddo
        if ( no_attribs.gt.0 .and. len_trim(xmldata).le.0 ) then
          write( info%lun, '(a)' ) '/>'
        else
          write( info%lun, '(a)') '>'//trim(adjustl(xmldata))//'</'//trim(adjustl(info%tag))//'>'
        endif
      endif
      end subroutine xml_put_element_

!-----------------------------------------------------------------------
! xml_put_close_tag_ --
!    Routine to write the closing tag
! Arguments:
!    info        Structure holding information on the XML-file
!-----------------------------------------------------------------------
      subroutine xml_put_close_tag_(info)
      type(XML_PARSE),  intent(inout)               :: info
      character(len=300), parameter :: indent = ' '
      info%level=info%level-1
      write( info%lun, '(4a)' ) indent(1:TABL*info%level), '</', trim(adjustl(info%tag)), '>'
      end subroutine xml_put_close_tag_

!-----------------------------------------------------------------------
! xml_put_data_ --
!    Routine to write the closing tag
! Arguments:
!    info        Structure holding information on the XML-file
!-----------------------------------------------------------------------
      subroutine xml_put_data_(info,xmldata)
      type(XML_PARSE),  intent(inout)               :: info
      character(len=*), intent(in)                  :: xmldata
      character(len=300), parameter :: indent = ' '
      if (len_trim(xmldata).gt.0) write( info%lun, '(2a)' ) indent(1:TABL*info%level),trim(adjustl(xmldata))
      end subroutine xml_put_data_

!-----------------------------------------------------------------------
! xml_replace_entities_ --
!    Routine to replace entities such as &gt; by their
!    proper character representation
! Arguments:
!    xmldata        Lines of character data found
!    no_data     (Nett) number of lines of character data
!-----------------------------------------------------------------------
      subroutine xml_replace_entities_( xmldata, no_data )
      character(len=*), intent(inout), dimension(:)    :: xmldata
      integer,          intent(inout)                  :: no_data

      integer :: i
      integer :: j
      integer :: j2
      integer :: k
      integer :: pos
      logical :: found

      do i = 1,no_data
      j = 1
      do
         do k = 1,size(entities,2)
            found = .false.
            pos   = index( xmldata(i)(j:), trim(entities(2,k)) )
            if ( pos .gt. 0 ) then
               found = .true.
               j     = j + pos - 1
               j2    = j + len_trim(entities(2,k))
               xmldata(i)(j:) = trim(entities(1,k)) // xmldata(i)(j2:)
               j     = j2
            endif
         enddo
         if ( .not. found ) exit
      enddo
      enddo

      end subroutine xml_replace_entities_

!-----------------------------------------------------------------------
! xml_options --
!    Routine to handle the parser options
! Arguments:
!    info                Structure holding information on the XML-file
!    ignore_whitespace   Ignore whitespace (leading blanks, empty lines) or not
!    no_data_truncation  Consider truncation of strings an error or not
!    report_lun          LU-number for reporting information
!    report_errors       Write messages about errors or not
!    report_details      Write messages about all kinds of actions or not
!-----------------------------------------------------------------------
      subroutine xml_options( info, ignore_whitespace, no_data_truncation,
     &                    report_lun, report_errors,
     &                    report_details )
      type(XML_PARSE),  intent(inout)               :: info
      logical, intent(in), optional                 :: ignore_whitespace
      logical, intent(in), optional                 :: no_data_truncation

      integer, intent(in), optional                 :: report_lun
      logical, intent(in), optional                 :: report_errors
      logical, intent(in), optional                 :: report_details

      if ( present(ignore_whitespace) ) then
        info%ignore_whitespace = ignore_whitespace
      endif
      if ( present(no_data_truncation) ) then
        info%no_data_truncation = no_data_truncation
      endif
      if ( present(report_lun) ) then
        report_lun_ = report_lun
      endif
      if ( present(report_errors) ) then
        report_errors_ = report_errors
      endif
      if ( present(report_details) ) then
        report_details_ = report_details
      endif
      end subroutine xml_options

!-----------------------------------------------------------------------
! xml_ok --
!    Function that returns whether all was okay or not
! Arguments:
!    info                Structure holding information on the XML-file
! Returns:
!    .true. if there was no error, .false. otherwise
!
!-----------------------------------------------------------------------
      logical function xml_ok( info )
      type(XML_PARSE),  intent(in)               :: info
      xml_ok =  .not. info%error
      end function xml_ok

!-----------------------------------------------------------------------
! xml_error --
!    Function that returns whether there was an error
! Arguments:
!    info                Structure holding information on the XML-file
! Returns:
!    .true. if there was an error, .false. if there was none
!-----------------------------------------------------------------------
      logical function xml_error( info )
      type(XML_PARSE),  intent(in)               :: info
      xml_error = info%error
      end function xml_error


!-----------------------------------------------------------------------
      integer function xml_find_attrib( attribs, no_attribs, name, value )
!-----------------------------------------------------------------------
      character(len=*), dimension(:,:)  :: attribs
      integer                           :: no_attribs
      character(len=*)                  :: name
      character(len=*)                  :: value
      integer :: i
      xml_find_attrib = -1
      do i = 1,no_attribs
      if ( name .eq. attribs(1,i) ) then
         value           = attribs(2,i)
         xml_find_attrib = i
         exit
      endif
      enddo
      end function xml_find_attrib


!-----------------------------------------------------------------------
      subroutine update_tagpath(info,TASK)
! maintain info about tag tree in tagpath
! parent of the current tag is stored in perenttag
!-----------------------------------------------------------------------
      type(XML_PARSE),  intent(inout) :: info
      character(len=*)   :: TASK
      integer IS,IL
      select case(TASK)
      case('open')
  ! get new parent tag = previous tag = 1st item on tagpath
        parenttag=' '
        call FINDSTRPAR(tagpath,':',1,IS,IL)
        if (IL.gt.0) parenttag=tagpath(IS:IS+IL-1)
  ! prepend tagpath with the new tag
        call PREPENDPAR(tagpath,':',trim(info%tag))
      case('close')
  ! remove 1st item from tagpath
        call FINDSTRPAR(tagpath,':',1,IS,IL)
        tagpath=trim(tagpath(IS+IL+1:))
  ! get new parent tag = 2nd item on tagpath
        parenttag=' '
        call FINDSTRPAR(tagpath,':',2,IS,IL)
        if (IL.gt.0) parenttag=tagpath(IS:IS+IL-1)
      end select
      info%tagpath = trim(tagpath)
      end subroutine update_tagpath

!-----------------------------------------------------------------------
! xml_process --
!    Routine to read the XML file as a whole and distribute processing
!    the contents over three user-defined subroutines
! Arguments:
!    info                structure with parsing info
!    attribs             Array for holding the attributes
!    xmldata             Array for holding the character data
!    startfunc           Subroutine to handle the start of elements
!    datafunc            Subroutine to handle the character data
!    endfunc             Subroutine to handle the end of elements
!    error               Indicates if there was an error or not
!-----------------------------------------------------------------------
      subroutine xml_process(startfunc,datafunc,endfunc)
      character(len=128), dimension(2,64) :: attribs
      character(len=XML_BUFFER_LENGTH)    :: xmldata
      character(len=XML_BUFFER_LENGTH)    :: line

      interface
      subroutine startfunc(tag,attribs,error )
        character(len=*) :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical  :: error
      end subroutine
      end interface

      interface
      subroutine datafunc(tag,xmldata,error)
        character(len=*) :: tag
        character(len=*) :: xmldata
        logical          :: error
      end subroutine
      end interface

      interface
      subroutine endfunc( tag, error )
        character(len=*) :: tag
        logical  :: error
      end subroutine
      end interface

      integer  :: noattribs

! these are defaults
! change any of them by calling xml_options or directly in locinfo
!      locinfo%lun = 10
!      locinfo%ignore_whitespace  = .false.
!      locinfo%no_data_truncation = .false.
!      locinfo%too_many_attribs   = .false.
!      locinfo%eof                = .false.
!      locinfo%error              = .false.
!      locinfo%level              = -1
!      locinfo%lineno             =  0
!      locinfo%comment_tag = .false.
!      locinfo%hasStart = .false.
!      locinfo%hasEnd = .false.
!      locinfo%hasData = .false.
!      report_lun_      = XML_STDOUT
!      report_errors_   = .true.
!      report_details_  = .false.

!      call xml_options( locinfo,report_details=.true. )

  ! store the current path in document tree as a sequence of tags separated by colon (top at the end)
      tagpath=' '
      parenttag=' '
      noattribs=0
      do while (xml_ok(locinfo).and.(.not. locinfo%eof))
        call xml_readline(locinfo,line)
      !  if (locinfo%lun.eq.6) write(*,*) 'xml_process ',locinfo%error,trim(line)
        if (locinfo%error) exit
        xmldata=' '
        call xml_get(locinfo,line,attribs,noattribs,xmldata)
  ! start tag
        if (xml_ok(locinfo).and.locinfo%hasStart) then
      ! update tagpath and parenttag
          call update_tagpath(locinfo,'open')
  !        if (trim(locinfo%tag).eq.'COMMAND') write(*,*) 'xml_process ',trim(line)
      ! call handler
          call startfunc(locinfo%tag, attribs(:,1:noattribs), locinfo%error)
        endif
  ! data
        if (xml_ok(locinfo).and.locinfo%hasData) then
          call datafunc(locinfo%tag, xmldata, locinfo%error)
        endif
  ! end tag
        if (xml_ok(locinfo).and.locinfo%hasEnd) then
          call endfunc(locinfo%tag, locinfo%error)
      ! get new current tag = former parenttag
          locinfo%tag=trim(parenttag)
      ! update tagpath and parenttag
          call update_tagpath(locinfo,'close')
        endif
    !    if (trim(locinfo%fname).eq.'STDIN') then
    !      write(*,*) 'locinfo: line=',locinfo%lineno,' tag=',trim(locinfo%tag),' err=',locinfo%error,' eof=',locinfo%eof
    !    endif
      !  if (len_trim(parenttag).le.0) exit
      enddo

      end subroutine xml_process



!-----------------------------------------------------------------------
! xml_suprocess --
!    Process a tag group in XML file
!    Use to provide specific handlers
! Usage: call xml_suprocess from the main startfunc of xml_process
! Should work recurrently (not tested)
!-----------------------------------------------------------------------
      subroutine xml_subprocess(startfunc,datafunc,endfunc)
      character(len=128), dimension(2,64) :: attribs
      character(len=XML_BUFFER_LENGTH)    :: xmldata
      character(len=XML_BUFFER_LENGTH)    :: line
      character(len=XML_TAG_LENGTH) :: root  ! root tag name
      logical :: stop_tag
      interface
      subroutine startfunc(tag,attribs,error )
        character(len=*) :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical  :: error
      end subroutine
      end interface

      interface
      subroutine datafunc(tag,xmldata,error)
        character(len=*) :: tag
        character(len=*) :: xmldata
        logical          :: error
      end subroutine
      end interface

      interface
      subroutine endfunc( tag, error )
        character(len=*) :: tag
        logical  :: error
      end subroutine
      end interface
      integer  :: noattribs

      noattribs=0
      ! xml_subprocess is called by start function in xml_process
      ! therefore locinfo already contains information to be handled by subprocess start function
      root=trim(locinfo%tag)
      call startfunc(locinfo%tag, attribs(:,1:noattribs), locinfo%error)
      stop_tag=.false.
      do while ((.not.stop_tag).and.xml_ok(locinfo).and.(.not. locinfo%eof))
        call xml_readline(locinfo,line)
        if (locinfo%error) exit
         xmldata=' '
        call xml_get(locinfo,line,attribs,noattribs,xmldata)
  ! start tag
        if (xml_ok(locinfo).and.locinfo%hasStart) then
      ! update tagpath and parenttag
          call update_tagpath(locinfo,'open')
      ! call handler
          call startfunc(locinfo%tag, attribs(:,1:noattribs), locinfo%error)
        endif
  ! data
        if (xml_ok(locinfo).and.locinfo%hasData) then
          call datafunc(locinfo%tag, xmldata, locinfo%error)
        endif
  ! end tag
        if (xml_ok(locinfo).and.locinfo%hasEnd) then
          stop_tag=(root.eq.trim(locinfo%tag))
          call endfunc(locinfo%tag, locinfo%error)
      ! get new current tag = former parenttag
          locinfo%tag=trim(parenttag)
      ! update tagpath and parenttag
          call update_tagpath(locinfo,'close')
        endif
      enddo

      end subroutine xml_subprocess

!-----------------------------------------------------------------------
      subroutine xml_process_file(filename,startfunc,datafunc,endfunc,error)
!-----------------------------------------------------------------------
      character(len=*),intent(in)       :: filename
      logical,intent(out)               :: error
      logical :: std  = .false.
 !     type(XML_PARSE) :: info

      interface
      subroutine startfunc(tag,attribs,error )
        character(len=*) :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical  :: error
      end subroutine
      end interface

      interface
      subroutine datafunc(tag,xmldata,error)
        character(len=*) :: tag
        character(len=*) :: xmldata
        logical          :: error
      end subroutine
      end interface

      interface
      subroutine endfunc( tag, error )
        character(len=*) :: tag
        logical  :: error
      end subroutine
      end interface

      call xml_options(locinfo, report_errors=.true., report_details=.false.)
      call xml_open(locinfo,filename,.true.)
      call xml_process(startfunc,datafunc,endfunc)
      if (locinfo%error) then
        write(*,*) 'XML parsing of ',trim(filename),': failed'
        call xml_report_errors(locinfo,'XML parsing stoped')
      ! else
      !  write(*,*) 'XML parsing of ',trim(filename),': O.K.'
      endif
      error=locinfo%error
      call xml_close(locinfo)
      end subroutine xml_process_file

!-----------------------------------------------------------------------
      subroutine PARSEINFO
!-----------------------------------------------------------------------
        write(*,*) 'PARSEINFO: tag=',trim(locinfo%tag),' line=',locinfo%lineno
      end subroutine PARSEINFO


      end module XMLPARSE
