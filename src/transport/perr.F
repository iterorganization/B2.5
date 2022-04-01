



      subroutine perr(n)
c--------------------------------------------------------------------
c     this subroutine does the error handling
c--------------------------------------------------------------------
      use b2mod_subsys
      implicit none

      integer ::
     *        n,ifile, iwarn1, iwarn2, iwarn3, iwarn4, iwarn5,
     +        iwarn6, iwarn7
      logical nlwarn

      data ifile / 0 /
      data iwarn1, iwarn2, iwarn3, iwarn4, iwarn5, iwarn6, iwarn7
     +    /     0,      0,      0,      0,      0,      0,      0/

      call subini ('perr')
      nlwarn = .false.
      if (n.le.10) then
        goto 1
      else
c       to switch off all warnings activate the line below
        if (.not.nlwarn) then
          call subend ()
          return
        end if
        if (ifile.eq.0) then
          open(11,file = 'neowarn.txt')
          ifile = 1
        endif
        goto 100
      endif
      call subend ()
      return

c-------------------------------------------------------------------
c     the error messages
c-------------------------------------------------------------------
    1 continue
      select case (n)
      case (1)
        write(*,*)'the parameters nsm and/or ncm are not consistently'
        write(*,*)'set throughout the code. they appear in different'
        write(*,*)'subroutines but must all have the same value'
        call xerrab ('Inconsistent parameters nsm and/or ncm')
      case (2)
        write(*,*)'sorry, you cannot run the code for rho = 0.'
        call xerrab ('Rho = 0 not allowed')
      case (3)
        write(*,*)'the parameter nmaxgr is not consistently'
        write(*,*)'set throughout the code. they appear in different'
        write(*,*)'subroutines but must have the same value'
        call xerrab ('Inconsistent parameter nmaxgr')
      case (4)
        write(*,*)'unable to reach convergence in viscos'
        call xerrab ('No convergence in viscos')
      case (5)
        write(*,*)'unable to reach convergence in the routine geom'
        call xerrab ('No convergence in geom')
      end select

c-------------------------------------------------------------------
c     the warning messages
c-------------------------------------------------------------------
  100 continue
      select case (n)
      case (11)
        if (iwarn1.eq.0) then
          write(11,*)'warning'
          write(11,*)'you have used nreg.ne.0. therefore a certain',
     +               ' collision'
          write(11,*)'regime will be forced. the results might be ',
     +               'wrong'
          write(11,*)
          iwarn1 = 1
        endif

      case (12)
        if (iwarn2.eq.0) then
          write(11,*)'warning'
          write(11,*)'you have used nleg.ne.3. the number of ',
     +               'legendre harmonics'
          write(11,*)'is not equal 3. this can make the ',
     +            'results inaccurate'
          write(11,*)
          iwarn2 = 1
        endif

      case (13)
        if (iwarn3.eq.0) then
          write(11,*)'warning'
          write(11,*)'you have used nenergy.ne.1 the energy ',
     +               'scattering in'
          write(11,*)'the calculation of the viscosity is not ',
     +               'taken into account'
          write(11,*)
          iwarn3 = 1
        endif

      case (14)
        if (iwarn4.eq.0) then
          write(11,*)'warning'
          write(11,*)'you have used ncofc.ne.1. ion-electron ',
     +               'collisions are'
          write(11,*)'not properly accounted for'
          write(11,*)
          iwarn4 = 1
        endif

      case (15)
        if (iwarn5.eq.0) then
          write(11,*)'warning'
          write(11,*)'one of the sigmas does not have the value ',
     +               '1. and one'
          write(11,*)'of the couplings in the pfirsch schlueter ',
     +               'calculation is'
          write(11,*)'therefore taken into account'
          write(11,*)
          iwarn5 = 1
        endif

      case (16)
        if (iwarn6.eq.0) then
          write(11,*)'warning'
          write(11,*)'both neogeo and neofrc are set true. however',
     +               ', the viscosity'
          write(11,*)'coefficients contain geometric dependencies',
     +               ' and must be'
          write(11,*)'newly calculated'
          write(11,*)
          iwarn6 = 1
        endif

      case (17)
        if (iwarn7.eq.0) then
          write(11,*)'warning'
          write(11,*)'the routine neoart has been called with ',
     +               'different rho/isel'
          write(11,*)'more than nkeep times. the values of the ',
     +               'geometry dependent'
          write(11,*)'quantities can no longer be kept in core'
          write(11,*)
          iwarn7 = 1
        endif

      end select
      call subend ()
      return

      end subroutine perr
