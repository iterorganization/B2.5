     



      subroutine perr(n)
c--------------------------------------------------------------------
c     this subroutine does the error handling
c--------------------------------------------------------------------

      implicit none
 
      integer ::
     *        n,ifile, iwarn1, iwarn2, iwarn3, iwarn4, iwarn5, 
     +        iwarn6, iwarn7
      logical nlwarn
      
      data ifile / 0 /
      data iwarn1, iwarn2, iwarn3, iwarn4, iwarn5, iwarn6, iwarn7
     +    /     0,      0,      0,      0,      0,      0,      0/
      
      nlwarn = .false. 
      if (n.le.10) then
        goto (1,2,3,4,5) n
      else
c       to switch off all warnings activate the line below
        if (.not.nlwarn) return
        if (ifile.eq.0) then
          open(11,file = 'neowarn.txt')
          ifile = 1
        endif
        goto (100,101,102,103,104,105,106)n-10
      endif
      return

c-------------------------------------------------------------------
c     the error messages
c-------------------------------------------------------------------
 1    continue
      write(*,*)'the parameters nsm and/or ncm are not consistently'
      write(*,*)'set throughout the code. they appear in different'
      write(*,*)'subroutines but must all have the same value'
      stop

 2    continue
      write(*,*)'sorry, you can not run the code for rho = 0.'
      stop

 3    continue 
      write(*,*)'the parameter nmaxgr is not consistently'
      write(*,*)'set throughout the code. they appear in different'
      write(*,*)'subroutines but must have the same value'
      stop

 4    continue
      write(*,*)'unable to reach convergence in viscos'
      stop

 5    continue
      write(*,*)'unable to reach convergence in the routine geom'
      stop
      
c-------------------------------------------------------------------
c     the warning messages
c-------------------------------------------------------------------
 100  continue
      if (iwarn1.eq.0) then
        write(11,*)'warning'
        write(11,*)'you have used nreg.ne.0. therefore a certain',
     +            ' collision'
        write(11,*)'regime will be forced. the results might be ',
     +            'wrong'
        write(11,*)
        iwarn1 = 1
      endif 
      return

 101  continue
      if (iwarn2.eq.0) then
        write(11,*)'warning'
        write(11,*)'you have used nleg.ne.3. the number of leg',
     +            'endre har-'
        write(11,*)'monics is not equal 3. this can make the ',
     +            'results inaccurate'
        write(11,*)
        iwarn2 = 1
      endif
      return

 102  continue
      if (iwarn3.eq.0) then
        write(11,*)'warning'
        write(11,*)'you have used nenergy.ne.1 the energy scat',
     +            'tering in'
        write(11,*)'the calculation of the viscosity is not ta',
     +            'ken into account'
        write(11,*)
        iwarn3 = 1
      endif
      return

 103  continue
      if (iwarn4.eq.0) then
        write(11,*)'warning'
        write(11,*)'you have used ncofc.ne.1. ion-electron col',
     +            'lisions are'
        write(11,*)'not properly accounted for'
        write(11,*)
        iwarn4 = 1
      endif
      return

 104  continue
      if (iwarn5.eq.0) then
        write(11,*)'warning'
        write(11,*)'one of the sigmas does not have the value ',
     +             '1. and one'
        write(11,*)'of the couplings in the pfirsch schlueter ',
     +             'calculation is'
        write(11,*)'therefore taken into account'
        write(11,*)
        iwarn5 = 1
      endif
      return

 105  continue
      if (iwarn6.eq.0) then
        write(11,*)'warning'
        write(11,*)'both neogeo and neofrc are set true. however',
     +             ', the viscosity'
        write(11,*)'coefficients contain geometric dependencies',
     +             ' and must be'
        write(11,*)'newly calculated'
        write(11,*)
        iwarn6 = 1
      endif
      return
 
 106  continue
      if (iwarn7.eq.0) then
        write(11,*)'warning'
        write(11,*)'the routine neoart has been called with dif',
     +             'ferent rho/isel'
        write(11,*)'more than nkeep times. the values of the ',
     +             'geometry dependent'
        write(11,*)'quantities can no longer be kept in core'
        write(11,*)
        iwarn7 = 1
      endif
      return
  
      end
