      subroutine test_xmlparse
      use xmlparse
      integer                                :: lunrep
      logical                                :: error
      external startfunc, datafunc, endfunc

      lunrep=6

!      write(*,*) 'test_xmlparse'
!      read(*,*)

      call xml_process_file('test.xml', startfunc, datafunc, endfunc, error )

      READ(*,*)
      end subroutine test_xmlparse

      subroutine startfunc( tag, attribs, error )
        character(len=*)                  :: tag
        character(len=*), dimension(:,:)  :: attribs
        logical                           :: error
        integer :: i
        write(*,*) trim(tag),'  start ',error
        write(*,*) '  attributes: ',size(attribs,2)
        do i=1,size(attribs,2)
          write(*,*) '  ',trim(attribs(1,i)),'=',trim(attribs(2,i))
        enddo
      end subroutine startfunc

      subroutine datafunc( tag, xmldata, error )
        character(len=*)                  :: tag
        character(len=*)                  :: xmldata
        logical                           :: error
        write(*,*) '  xmldata: ',trim(tag),':',trim(xmldata),error
      end subroutine datafunc

      subroutine endfunc( tag, error )
         character(len=*)                  :: tag
         logical                           :: error
         write(*,*) trim(tag),'  end ',error
      end subroutine endfunc
