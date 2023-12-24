      subroutine b2pgeom(lun)
c
c  version : 07.09.2000 16:07
c=======================================================================
c***  Produces tables with coordinates of grid cell corners along the
c***  separatrix and along the magnetic surfaces
c***  for subsequent plotting with standalone software
c***
c*** Output files:
c***
c***  gm_sep.dat first separatrix
c***  gm_sep_i.dat inner second separatrix for DDN
c***  gm_sep_o.dat outer second separatrix for DDN
c***  gm_SOL.dat magnetic surfaces in first SOL
c***  gm_SOL_i.dat magnetic surfaces in inner second SOL for DDN
c***  gm_SOL_o.dat magnetic surfaces in outer second SOL for DDN
c***  gm_PFR.dat magnetic surfaces in lower PFR
c***  gm_PFR_u.dat magnetic surfaces in upper PFR for DN
c***  gm_core.dat magnetic surfaces in the core part of the grid
c***
c*** For CDN, no gm_sep.dat or gm_SOL.dat is created.
c=======================================================================
      use b2us_data
      use b2mod_diag
      use b2mod_b2plot
      implicit none
      integer lun

c*** Local variables

      integer i,j,l,i1,i2
      character hlp_frmr*100, tab, ofile*12
c=======================================================================

      tab=char(9)

      if(mpg%nXpt.eq.1) then !{                      === single null ===
#ifdef WG_TODO
        k=jsep

        ofile='gm_sep.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# Separatrix shape'
        write(lun,'(2a)',err=998) '# ',trim(label)
        write(lun,'(3a)',err=998) '# R',tab,'Z'
        do i=0,nx !{
          write(lun,'(1p,e15.7,a,e15.7)',err=998) crx(i,k,2),tab,
     ,                                                        cry(i,k,2)
        end do !}
        close(lun)

        ofile='gm_SOL.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# SOL magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        write(hlp_frmr,'(a,i3,a)')
     ,    '(a,4x,',ny-k,'(a1,i2.2,13x,a1,i2.2,13x))'
        write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=k,ny-1)
        write(hlp_frmr,'(a,i3,a)') '(1p,',2*(ny-k),'(e14.7,2x))'
        do i=0,nx !{
          write(lun,hlp_frmr,err=998) (crx(i,j,2),cry(i,j,2),j=k,ny-1)
        end do !}
        close(lun)

        ofile='gm_PFR.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# PFR magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        write(hlp_frmr,'(a,i3,a)')
     ,   '(a,4x,',k+1,'(a1,i2.2,13x,a1,i2.2,13x))'
        write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=0,k)
        write(hlp_frmr,'(a,i3,a)') '(1p,',2*(k+1),'(e14.7,2x))'
        i1=0
        i2=leftcut(1)-1
        do l=1,2 !{
          do i=i1,i2 !{
            write(lun,hlp_frmr,err=998) (crx(i,j,2),cry(i,j,2),j=0,k)
          end do !}
          write(lun,hlp_frmr,err=998) (crx(i2,j,3),cry(i2,j,3),j=0,k)
          i1=rightcut(1)
          i2=nx-1
        end do !}
        close(lun)

        ofile='gm_core.dat'
        open(lun,file=ofile,err=999)
        rewind(lun)
        write(lun,' (a)',err=998) '# Core magnetic surfaces'
        write(lun,'(2a)',err=998) '# ',trim(label)
        write(hlp_frmr,'(a,i3,a)')
     ,    '(a,4x,',k+1,'(a1,i2.2,13x,a1,i2.2,13x))'
        write(lun,hlp_frmr,err=998) '#',('R',i,'Z',i,i=0,k)
        write(hlp_frmr,'(a,i3,a)') '(1p,',2*(k+1),'(e14.7,2x))'
        i1=leftcut(1)
        i2=rightcut(1)-1
        do i=i1,i2 !{
          write(lun,hlp_frmr,err=998) (crx(i,j,2),cry(i,j,2),j=0,k)
        end do !}
        write(lun,hlp_frmr,err=998) (crx(i2,j,3),cry(i2,j,3),j=0,k)

        close (lun)
#endif
      else if(mpg%nXpt.eq.2) then !}{ === double null ===

        write(0,*) 'b2pgeom: double null configuration yet to be ',
     ,                                                    'implemented'
        return

      else !}{
        write(0,*) 'b2pgeom: wrong value of nXpt (must be 1 or 2) ',
     ,             mpg%nXpt
        return
      end if !}

      return
c=======================================================================
  998 write(0,*) 'b2pgeom: error writing to the file ',ofile
      return
c=======================================================================
  999 write(0,*) 'b2pgeom: unable to open file ',ofile
      return
c=======================================================================
      end subroutine b2pgeom
