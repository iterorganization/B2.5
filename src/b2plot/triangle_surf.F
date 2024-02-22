c
      subroutine triangle_surf()
      use b2mod_types
      use b2us_data
      use b2mod_b2plot
      use b2mod_b2plot_debug
      use b2mod_b2plot_values
#ifdef B25_EIRENE
      use eirmod_ctrig
      implicit none
      character*256 tmpstring,titlestring
      integer is

      if(b2plot_debug.ge.10)
     > write(*,*) 'Plotting Eirene array: ntrii, nz = ',
     1 nedim(1,iewrk),nedim(2,iewrk)

      do is=0,nedim(2,iewrk)-1
        tmpstring=trim(eirene_label(iewrk))//textewk(is,iewrk)
        call customise_label(tmpstring,titlestring)
        if(b2plot_debug.ge.7) write(*,*) 'physical domain'
        call init_page(titlestring,
     1   values(5),values(6),values(7),values(8),1,.true.,.true.,
     2   .true.)
        call plottripatch(mpg%nCv,ntrii,nrknot,is,xtrian,ytrian,necke,
     3   eirene_stack(:,is,iewrk),
     4   values(38),values(39),linlog)
        if(earrowl .gt. 0.0_R8) call plot_trig_arrows(is)
        call end_page
      enddo
#endif
      return
      end subroutine triangle_surf
