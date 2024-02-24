c
      subroutine cont
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
      use b2mod_boundary_namelist
      implicit none
      character*256 tmpstring, titlestring
      real (kind=R8) :: wrkx(mpg%nCv), wrky(mpg%nCv)
      integer is

      if(b2plot_debug.ge.10) write(*,*) 'Plotting array: nCv, nz = ',
     1 ndim(1,iwrk),ndim(2,iwrk)
      wrkx(:) = geo%cvX(:)
      wrky(:) = geo%cvY(:)
      do is=0,ndim(2,iwrk)-1
        tmpstring=trim(graphlabel(iwrk))//textwk(is,iwrk)
        call customise_label(tmpstring,titlestring)
        if(typeofplot .ge. 2 .and. mpg%isClassicalGrid .eq. 1) then
          if(b2plot_debug.ge.7) write(*,*) 'computational domain'
          call init_page(titlestring,
     1     values(1),values(2),values(3),values(4),1,.false.,.true.,
     3     .false.)
          call plotcontour(mpg%nCv,is,
     2     xx(1,1),yy(1,1),
     3     work(1,is,iwrk),
     4     values(9),values(10),.false.,linlog)
          call end_page
        endif
        if(mod(typeofplot,2) .eq. 1) then
          if(b2plot_debug.ge.7) write(*,*) 'physical domain'
          call init_page(titlestring,
     1     values(5),values(6),values(7),values(8),1,.true.,.true.,
     2     .true.)
          call plotcontour(mpg%nCv,is,
     2     wrkx,wrky,
     3     work(1,is,iwrk),
     4     values(9),values(10),.true.,linlog)
          if(arrowl .gt. 0.0_R8) call plot_arrows(is)
          call end_page
        endif
      enddo
      return
      end subroutine cont
