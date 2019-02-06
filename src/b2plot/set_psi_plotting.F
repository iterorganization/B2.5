c
      subroutine set_psi_plotting(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        plot_vs_psi=.not.plot_vs_psi
        write (*,*) 'Plotting versus psi option changed'
      else
        plot_vs_psi=.false.
        write (*,*) 'Not plotting versus psi'
      endif
      return
      end
