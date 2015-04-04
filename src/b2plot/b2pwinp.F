      subroutine b2pwinp
c
c  version : 07.11.97 13:19
c
c======================================================================
c*** Data input from dg files.
c*** MUST be compiled in single precision!!!
c======================================================================
      use b2mod_types
      use b2mod_b2plot_wall_loading
      implicit none
#     include "impcon.inc"
      integer lun
      parameter (lun=96)
      integer i
      logical open_file
      character*16 ufile
      external b2pwinm, open_file
c======================================================================
c*** This routine works once
c
      if(l2pwinp) return
c======================================================================
      l2pwinp=.true.
c
c*** Set the default values
c
      nwldprt=0
      do i=1,nwldprtx
        nwldlis(i)=0
        lwldstrt(i)=0
      end do
      lppp=0
      jedge(1:4)=0
      add_rad=0.0_R8
      nadsrc=0
c
c*** Set up the input routine
c
      immed=.false.
      mimp=lun
      ufile='param.dg'
c      open(lun,file=ufile,status='old',err=998)
      if(.not.open_file(lun,ufile,'/modules/Uinp/')) then
        l2pwinp=.false.
        write(*,*) 'b2pwinp: unable to open param.dg'
        return
      end if
c
c*** ... read the data from the DG output file ...
c
      call import(b2pwinm)
      close(lun)
c
c*** ... and prepare the target and shadowing lists
c
      call b2pwlprp
c
      return
c======================================================================
 998  write(*,*) 'Error opening the file ',ufile
      stop
c======================================================================
      end subroutine b2pwinp
