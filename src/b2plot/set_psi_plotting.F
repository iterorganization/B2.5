c
      subroutine set_norm_psi_plotting(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        plot_vs_norm_psi=.not.plot_vs_norm_psi
        write (*,*) 'Plotting versus normalised psi option changed'
      else
        plot_vs_norm_psi=.false.
        write (*,*) 'Not plotting versus normalised psi'
      endif
      if (plot_vs_norm_psi) then
        if (plot_vs_full_psi) call set_full_psi_plotting(.false.)
        if (plot_vs_root_psi) call set_root_psi_plotting(.false.)
        if (plot_vs_rmaj) call set_rmaj_plotting(.false.)
      end if
      return
      end subroutine set_norm_psi_plotting
c
      subroutine set_full_psi_plotting(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        plot_vs_full_psi=.not.plot_vs_full_psi
        write (*,*) 'Plotting versus full psi option changed'
      else
        plot_vs_full_psi=.false.
        write (*,*) 'Not plotting versus full psi'
      endif
      if (plot_vs_full_psi) then
        if (plot_vs_norm_psi) call set_norm_psi_plotting(.false.)
        if (plot_vs_root_psi) call set_root_psi_plotting(.false.)
        if (plot_vs_rmaj) call set_rmaj_plotting(.false.)
      end if
      return
      end subroutine set_full_psi_plotting
c
      subroutine set_root_psi_plotting(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        plot_vs_root_psi=.not.plot_vs_root_psi
        write (*,*) 'Plotting versus root psi option changed'
      else
        plot_vs_root_psi=.false.
        write (*,*) 'Not plotting versus root psi'
      endif
      if (plot_vs_root_psi) then
        if (plot_vs_norm_psi) call set_norm_psi_plotting(.false.)
        if (plot_vs_full_psi) call set_full_psi_plotting(.false.)
        if (plot_vs_rmaj) call set_rmaj_plotting(.false.)
      end if
      return
      end subroutine set_root_psi_plotting
c
      subroutine set_rmaj_plotting(onflag)
      use b2mod_types
      use b2mod_b2plot
      logical onflag
      if(onflag) then
        plot_vs_rmaj=.not.plot_vs_rmaj
        write (*,*) 'Plotting versus major radius option changed'
      else
        plot_vs_root_psi=.false.
        write (*,*) 'Not plotting versus major radius'
      endif
      if (plot_vs_rmaj) then
        if (plot_vs_norm_psi) call set_norm_psi_plotting(.false.)
        if (plot_vs_full_psi) call set_full_psi_plotting(.false.)
        if (plot_vs_root_psi) call set_root_psi_plotting(.false.)
      end if
      return
      end subroutine set_rmaj_plotting
