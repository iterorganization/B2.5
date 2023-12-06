      subroutine maskfun(nCv, compval, exval, toll, maskf)
c     This routine ####
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_types
      use b2mod_math
      use b2mod_indirect
      use b2mod_constants
      use b2mod_subsys

      implicit none

*   ..input arguments (unchanged on exit)
      integer, intent(in) :: nCv
      real (kind=R8), intent(in) :: compval(nCv), exval(nCv), toll
c     Output variables
      integer, intent(inout):: maskf(nCv)

c     Internal variables
      integer :: iCv
      real (kind=R8) :: rel(nCv)
      character chns*3


      call subini ('maskfun')

      if(SIZE(compval).ne.SIZE(exval)) then
        call xerrab('Dimensional mismatch between array sizes of'//
     &            'values to compare w.r.t exact values')
      endif

      call ifill (nCv, 0, maskf, 1)
      call sfill (nCv, 0.0_R8, rel, 1)

      do iCv = 1, nCv
        rel(iCv) = abs((compval(iCv) - exval(iCv))/exval(iCv))
        if (rel(iCv).le.toll) then
          maskf(iCv) = 1
        endif
      enddo

*   ..Check whether all the elements of mask function are zero and print warning
      if (ALL( maskf==0 )) then
        write(*,*) 'Warning: maskf is zero everywhere'
      endif

*   ..return
      call subend ()
      return
      end subroutine maskfun
