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
      use b2mod_subsys
      implicit none
#     include "impcon.inc"
      integer lun
      parameter (lun=96)
      integer i
      logical open_file, file_ok
      character*24 ufile
      external b2pwinm, open_file, xerrab
c======================================================================
c*** This routine works once
c
      if(l2pwinp.or..not.shadow) return
c======================================================================
      call subini ('b2pwinp')
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
      jedge(1:2*ntrgsx)=0
      add_rad=0.0_R8
      nadsrc=0
      do i=1,ntrgsx
        ltrgedge(1,i)=0
        ltrgedge(2,i)=0
        mtrgedge(1,i)=0
        mtrgedge(2,i)=0
      end do
      ntrgsa=0
      nskipel=0
      lnumboff=0
      lvdg20=.false.
c
c*** Set up the input routine
c
      lena=8
      lfin=72
      jmess=2
      immed=.false.
      mimp=lun
      keylen=8
      fixkey=.false.
      finish = 'finish'
      promp(1) = repeat(' ',8)
      promp(2) = repeat(' ',8)
      ufile='param.dg'
      call find_file(ufile,file_ok)
      if(.not.open_file(lun,ufile,'modules/Uinp')) then
        l2pwinp=.false.
        shadow=.false.
        write(*,*) 'b2pwinp: unable to open param.dg'
        call subend ()
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
      call b2pwlprp_nodiff
c
      call subend ()
      return
c======================================================================
      end subroutine b2pwinp
