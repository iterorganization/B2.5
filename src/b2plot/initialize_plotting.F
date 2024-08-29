      subroutine initialise_plotting
      use b2mod_types
      use b2mod_b2plot_debug
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
#include <graphics_vars.h>

c     the zoom common
      integer winid
      integer wsnshwin,wsnX11
      common /zoom_com/ winid,wsnshwin,wsnX11

#ifdef NCAR4
      character*256 devenv
      integer iw
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
*     ------------------------------------------------------------------
#ifdef PLPLOT
      call color
      if(device.ne.'xwin') call plsfnam('b2plot.ps')
      call plstart(device,nx_page,ny_page)
c      call plinit()
#else
      ncol=12
#ifdef NCAR4
      devenv = 'x11 ps'
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('B2PLOT_DEV', 0, devenv, lenval, ierror)
#else
      call get_environment_variable('B2PLOT_DEV',
     . status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('B2PLOT_DEV',
     . value=devenv)
#endif
#endif
      if(b2plot_debug.ge.8) write(*,*) 'devenv: ',trim(devenv)
      nwork=0
      call gopks(6,0)
      if(index(devenv,'x11').ne.0.or.
     1   index(devenv,'X11').ne.0) then
        if(b2plot_debug.ge.12) write(*,*) 'X11'
        nwork=nwork+1
#ifdef WINDOWS
        call gopwk(nwork,0,1)      ! X11          !srv 07.05.01
#else
        call gopwk(nwork,0,8)      ! X11
#endif
        wsnX11=nwork
      else
        wsnX11=0
      endif
      if(index(devenv,'shwin').ne.0) then
        if(b2plot_debug.ge.12) write(*,*) 'shwin'
        nwork=nwork+1
        call gopwk(nwork,winid,7)  ! shwin
        wsnshwin=nwork
      else
        wsnshwin=0
      endif
      if(index(devenv,'ps').ne.0.or.
     1   index(devenv,'PS').ne.0) then
        if(b2plot_debug.ge.12) write(*,*) 'ps'
        nwork=nwork+1
        call ngsetc('ME','b2plot.ps')
        if(b2plot_debug.ge.13) write(*,*) 'ngsetc'
        if(xwidth.lt.1) then       ! portrait
          call ngseti('lx',50)
          call ngseti('ly',50)
          call ngseti('ux',792)
          call ngseti('uy',792)
#ifdef WINDOWS
          call gopwk(nwork,2,2700) ! PS       !srv 07.05.01
#else
          call gopwk(nwork,0,20)   ! PS
#endif
        elseif(yheight.lt.1) then  ! landscape
          call ngseti('lx',50)
          call ngseti('ly',50)
          call ngseti('ux',792)
          call ngseti('uy',792)
#ifdef WINDOWS
          call gopwk(nwork,2,2701) ! PS       !srv 07.05.01
#else
          call gopwk(nwork,0,26)   ! PS
#endif
        else
#ifdef WINDOWS
          call gopwk(nwork,2,2700) ! PS      !srv 07.05.01
#else
          call gopwk(nwork,0,20)   ! PS
#endif
        endif
      endif
      if(index(devenv,'cgm').ne.0.or.
     1   index(devenv,'CGM').ne.0.or.
     2   nwork.eq.0) then
        if(b2plot_debug.ge.12) write(*,*) 'cgm'
        nwork=nwork+1
        call gopwk(nwork,2,1)      ! CGM
      endif
      do iw=1,nwork
        call gacwk(iw)
      enddo
      if(b2plot_debug.ge.8) write(*,*) 'nwork ',nwork
#else
      call opngks
#endif
      call new_color_map
#endif
      initialised_plotting=.true.
      return
      end subroutine initialise_plotting
c
      subroutine initialise_plotting_flags(time)
      use b2mod_types
      implicit none
      real (kind=R8), intent(in) :: time
#include <graphics_vars.h>
c
      extralabel=' '
      globalheader=' '
      xlabel=' '
      ylabel=' '
      runtime=real(time,R4)
      glsz=1.0_R4
      lbsz=1.0_R4
      ndec=10
      initialised_plotting=.false.
      want_aspect=.true.
      lzero=.false.
      label_same_lines=.true.
      color_map_loaded=.false.
      color_file=' '
      device='psc'
      bold=.false.
      return
      end subroutine initialise_plotting_flags

      subroutine set_color_map(sval)
      use b2mod_types
      implicit none
      character*(*), intent(in) :: sval
#include <graphics_vars.h>
      logical streql
      external streql
c
      color_file=sval
      if (streql(color_file,' ')) then
        write(*,*) 'Setting default color table'
      else
        write(*,*) 'Setting color table ', trim(color_file)
      endif
      color_map_loaded=.false.
      if (nwork.gt.0) call new_color_map
      return
      end subroutine set_color_map

      subroutine set_fill_option(iopt_fill)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt_fill
#include <graphics_vars.h>
c
      fill_option=iopt_fill
      return
      end subroutine set_fill_option

      subroutine initialise_plotting_0
      use b2mod_types
      implicit none
#include <graphics_vars.h>
c
      dirty_page=.false.
      new_page=.true.
      same_page=0
      noax=1
      write(*,*) 'same_page set to ',same_page
      return
      end subroutine initialise_plotting_0

      subroutine set_glsz(rval)
      use b2mod_types
      implicit none
      real (kind=R8), intent(in) :: rval
#include <graphics_vars.h>
c
      glsz=real(rval,R4)
      return
      end subroutine set_glsz

      subroutine set_lbsz(rval)
      use b2mod_types
      implicit none
      real (kind=R8), intent(in) :: rval
#include <graphics_vars.h>
c
      lbsz=real(rval,R4)
      return
      end subroutine set_lbsz

      subroutine set_ndec(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      ndec=iopt
      return
      end subroutine set_ndec

      subroutine get_ndec(iopt)
      use b2mod_types
      implicit none
      integer, intent(out) :: iopt
#include <graphics_vars.h>
c
      iopt=ndec
      return
      end subroutine get_ndec

      subroutine set_noax
      use b2mod_types
      implicit none
#include <graphics_vars.h>
c
      noax=-noax
      return
      end subroutine set_noax

      subroutine swap_aspect
      use b2mod_types
      implicit none
#include <graphics_vars.h>
c
      want_aspect=.not.want_aspect
      write(*,*) 'want_aspect switched to ',want_aspect
      return
      end subroutine swap_aspect

      subroutine set_zero
      use b2mod_types
      implicit none
#include <graphics_vars.h>
c
      lzero=.not.lzero
      write(*,*) 'lzero switched to ',lzero
      return
      end subroutine set_zero

      subroutine toggle_label_lines
      use b2mod_types
      implicit none
#include <graphics_vars.h>
c
      label_same_lines=.not.label_same_lines
      return
      end subroutine toggle_label_lines

      subroutine set_device(devnam)
      use b2mod_types
      implicit none
      character*256, intent(in) :: devnam
#include <graphics_vars.h>
c
      device=devnam
      return
      end subroutine set_device

      subroutine set_line_font_size(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      lfsz=iopt
      return
      end subroutine set_line_font_size

      subroutine set_axes_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      zcol=iopt
      return
      end subroutine set_axes_color

      subroutine set_header_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      hcol=iopt
      return
      end subroutine set_header_color

      subroutine set_title_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      gcol=iopt
      return
      end subroutine set_title_color

      subroutine set_frame_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      fcol=iopt
      return
      end subroutine set_frame_color

      subroutine set_barlabel_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      bcol=iopt
      return
      end subroutine set_barlabel_color

      subroutine set_extralabel_color(iopt)
      use b2mod_types
      implicit none
      integer, intent(in) :: iopt
#include <graphics_vars.h>
c
      xcol=iopt
      return
      end subroutine set_extralabel_color

      subroutine new_color_map
      use b2mod_types
      use b2mod_b2plot_debug
      implicit none
#include <graphics_vars.h>
      real red, green, blue
      integer iw
      logical open_file, streql
      external open_file, streql

      if(color_map_loaded) return
      if(.not.streql(color_file,' ')) then
        write(*,*) 'Using ',trim(color_file)
        if(.not.open_file(99,color_file,
     1   'data.local/palettes|scripts.local/palettes|scripts/palettes'))
     >   then
          call color
          goto 50
        endif
        ncol=0
   10   read(99,*,end=40,err=30) red,green,blue
        do iw=1,nwork
          call gscr(iw,ncol,red,green,blue)
          if(b2plot_debug.ge.13.and.iw.eq.1)
     >     write(*,'(a,i2,a,1p,3(f7.4,1x))')
     >     ' Color ',ncol,': ',red,green,blue
        enddo
        ncol=ncol+1
        goto 10
   30   write(*,*) 'Error in reading ',trim(color_file),
     1   ' on line ',ncol+1
        call color
        goto 50
   40   close(99)
        ncol=ncol-2
        if(ncol.le.0) then
          write(*,*) 'Too few colors in ',trim(color_file)
          call color
          goto 50
        endif
      else
C Set up color table
        call color
c      call color2(ncol)
      endif
   50 continue
      color_map_loaded=.true.
      return
      end subroutine new_color_map
