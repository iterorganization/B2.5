C ================= User-defined colors =========================
      SUBROUTINE COLOR
C ===============================================================
      use b2mod_types
      use b2mod_b2plot_debug
      integer iw
#include <graphics_vars.h>
#ifdef PLPLOT
      call plscmap0n(18)
      ncol=17-2
      call plcolset(0,255,255,255)
      call plcolset(1,0,0,0)
      call plcolset(2,204,204,255)
      call plcolset(3,153,0,204)
      call plcolset(4,102,89,204)
      call plcolset(5,63,114,242)
      call plcolset(6,0,191,255)
      call plcolset(7,0,229,255)
      call plcolset(8,51,204,51)
      call plcolset(9,51,255,127)
      call plcolset(10,127,255,0)
      call plcolset(11,178,255,51)
      call plcolset(12,255,255,0)
      call plcolset(13,255,216,0)
      call plcolset(14,255,165,0)
      call plcolset(15,255,0,51)
      call plcolset(16,229,63,0)
      call plcolset(17,216,114,204)
#else
      do iw=1,nwork
C
C     BACKGROUND COLOR
C White
        if(b2plot_debug.ge.13) write(*,*) 'Color 0: White'
        CALL GSCR(iw,  0, 1.0, 1.0, 1.0)
C
C     FOREGROUND COLORS
C Black
        if(b2plot_debug.ge.13) write(*,*) 'Color 1: Black'
        CALL GSCR(iw,  1, 0., 0., 0.)
        ncol=17-2
C Orchid
        if(b2plot_debug.ge.13) write(*,*) 'Color 17: Orchid'
        CALL GSCR(iw, 17, 0.85, 0.45, 0.8)
C Red
        if(b2plot_debug.ge.13) write(*,*) 'Color 16: Red'
        CALL GSCR(iw, 16, 0.9, 0.25, 0.0)
C OrangeRed
        if(b2plot_debug.ge.13) write(*,*) 'Color 15: OrangeRed'
        CALL GSCR(iw, 15, 1.0, 0.0, 0.2)
C Orange
        if(b2plot_debug.ge.13) write(*,*) 'Color 14: Orange'
        CALL GSCR(iw, 14, 1.0, 0.65, 0.0)
C Golden
        if(b2plot_debug.ge.13) write(*,*) 'Color 13: Golden'
        CALL GSCR(iw, 13, 1.0, 0.85, 0.0)
C Yellow
        if(b2plot_debug.ge.13) write(*,*) 'Color 12: Yellow'
        CALL GSCR(iw, 12, 1.0, 1.0, 0.0)
C GreenYellow
        if(b2plot_debug.ge.13) write(*,*) 'Color 11: GreenYellow'
        CALL GSCR(iw, 11, 0.7, 1.0, 0.2)
C Chartreuse
        if(b2plot_debug.ge.13) write(*,*) 'Color 10: Chartreuse'
        CALL GSCR(iw, 10, 0.5, 1.0, 0.0)
C Celeste
        if(b2plot_debug.ge.13) write(*,*) 'Color 9: Celeste'
        CALL GSCR(iw,  9, 0.2, 1.0, 0.5)
C Green
        if(b2plot_debug.ge.13) write(*,*) 'Color 8: Green'
        CALL GSCR(iw,  8, 0.2, 0.8, 0.2)
C Aqua
        if(b2plot_debug.ge.13) write(*,*) 'Color 7: Aqua'
        CALL GSCR(iw,  7, 0.0, 0.9, 1.0)
C DeepSkyBlue
        if(b2plot_debug.ge.13) write(*,*) 'Color 6: DeepSkyBlue'
        CALL GSCR(iw,  6, 0.0, 0.75, 1.0)
C RoyalBlue
        if(b2plot_debug.ge.13) write(*,*) 'Color 5: RoyalBlue'
        CALL GSCR(iw,  5, 0.25, 0.45, 0.95)
C SlateBlue
        if(b2plot_debug.ge.13) write(*,*) 'Color 4: SlateBlue'
        CALL GSCR(iw,  4, 0.4, 0.35, 0.8)
C DarkViolet
        if(b2plot_debug.ge.13) write(*,*) 'Color 3: DarkViolet'
        CALL GSCR(iw,  3, 0.6, 0.0, 0.8)
C Lavender
        if(b2plot_debug.ge.13) write(*,*) 'Color 2: Lavender'
        CALL GSCR(iw,  2, 0.8, 0.8, 1.0)
C Done.
      enddo
#endif
      RETURN
      END
c
      integer function no_available_colors(ityp)
      use b2mod_types
      integer ityp
#include <graphics_vars.h>
      if(partition_colors.ne.0) then
      if(ityp.eq.0) then
        no_available_colors=NCOL*partition_colors
      else
        no_available_colors=NCOL*(1-partition_colors)
      endif
      else
      no_available_colors=NCOL
      endif
      end
c
      function choose_color(icol,ityp)
      use b2mod_types
      integer icol, ityp
#include <graphics_vars.h>
      integer choose_color
      if(partition_colors.ne.0) then
      if(ityp.eq.0) then
        if(fill_option.eq.0) then
          choose_color=(icol)*NCOL*partition_colors/(ncont+3)+2
        else
          choose_color=MOD(icol,int(NCOL*partition_colors))+2
        endif
      else
        if(fill_option.eq.0) then
          choose_color=(icol)*NCOL*(1-partition_colors)/(ncont+3)+2+
     1     (NCOL*partition_colors)
        else
          choose_color=MOD(icol,int(NCOL*(1-partition_colors)))+2+
     1     (NCOL*partition_colors)
        endif
      endif
      else
      if(fill_option.eq.0) then
        choose_color=(icol)*NCOL/(ncont+3)+2
      else
        choose_color=MOD(icol,NCOL)+2
      endif
      endif
      return
      end
c
      subroutine set_partition_colors(val)
      use b2mod_types
#include <graphics_vars.h>
      real (kind=R8) :: val
      if(val.lt.0.0.or.val.gt.1.0) then
      write(*,*) 'set_partition_colors: parameter out of range'
      partition_colors=1.0
      else
      partition_colors=val
      endif
      return
      end

