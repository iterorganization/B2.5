      module b2mod_exp_file_names
	 character(32) :: fnm,
     &                  fnm_power='pw'
	 logical(1) ::    l_power=.false.

	 character(32) :: fnm_ne_outer_mdpl='ne_outer_mdpl',
     &                  fnm_te_outer_mdpl='te_outer_mdpl', 
     &                  fnm_ti_outer_mdpl='ti_outer_mdpl'
	 logical(1) ::    l_ne_outer_mdpl=.false.,
     &                  l_te_outer_mdpl=.false., 
     &                  l_ti_outer_mdpl=.false.

	 character(32) :: fnm_ne_inner_lower='ne_inner_lower',
     &                  fnm_te_inner_lower='te_inner_lower', 
     &                  fnm_ti_inner_lower='ti_inner_lower',
     &                  fnm_js_inner_lower='js_inner_lower'
	 logical(1) ::    l_ne_inner_lower=.false.,
     &                  l_te_inner_lower=.false., 
     &                  l_ti_inner_lower=.false.,
     &                  l_js_inner_lower=.false.

	 character(32) :: fnm_ne_inner_upper='ne_inner_upper',
     &                  fnm_te_inner_upper='te_inner_upper', 
     &                  fnm_ti_inner_upper='ti_inner_upper',
     &                  fnm_js_inner_upper='js_inner_upper'
	 logical(1) ::    l_ne_inner_upper=.false.,
     &                  l_te_inner_upper=.false., 
     &                  l_ti_inner_upper=.false.,
     &                  l_js_inner_upper=.false.

	 character(32) :: fnm_ne_outer_lower='ne_outer_lower',
     &                  fnm_te_outer_lower='te_outer_lower', 
     &                  fnm_ti_outer_lower='ti_outer_lower',
     &                  fnm_js_outer_lower='js_outer_lower'
	 logical(1) ::    l_ne_outer_lower=.false.,
     &                  l_te_outer_lower=.false., 
     &                  l_ti_outer_lower=.false.,
     &                  l_js_outer_lower=.false.

	 character(32) :: fnm_ne_outer_upper='ne_outer_upper',
     &                  fnm_te_outer_upper='te_outer_upper', 
     &                  fnm_ti_outer_upper='ti_outer_upper',
     &                  fnm_js_outer_upper='js_outer_upper'
	 logical(1) ::    l_ne_outer_upper=.false.,
     &                  l_te_outer_upper=.false., 
     &                  l_ti_outer_upper=.false.,
     &                  l_js_outer_upper=.false.


      end module b2mod_exp_file_names
