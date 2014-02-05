      module b2mod_optim
       use b2mod_types

	 integer :: nfunc, istyle
       integer, parameter :: nx = 116, ny =21, ns = 2
       real (kind=R8) :: power_mod
*   .. for inner lower
	 real (kind=R8) :: 
     &      cy_ne_inner_lower(-1:ny), ne_inner_lower_mod(-1:ny), 
     &       b_ne_inner_lower_mod(-1:ny), c_ne_inner_lower_mod(-1:ny),
     &       d_ne_inner_lower_mod(-1:ny),
     &      cy_te_inner_lower(-1:ny), te_inner_lower_mod(-1:ny),
     &       b_te_inner_lower_mod(-1:ny), c_te_inner_lower_mod(-1:ny),
     &       d_te_inner_lower_mod(-1:ny),
     &      cy_ti_inner_lower(-1:ny), ti_inner_lower_mod(-1:ny),
     &       b_ti_inner_lower_mod(-1:ny), c_ti_inner_lower_mod(-1:ny),
     &       d_ti_inner_lower_mod(-1:ny),
     &      cy_js_inner_lower(-1:ny), js_inner_lower_mod(-1:ny),
     &       b_js_inner_lower_mod(-1:ny), c_js_inner_lower_mod(-1:ny),
     &       d_js_inner_lower_mod(-1:ny)

*   .. for inner outer
	 real (kind=R8) :: 
     &      cy_ne_inner_upper(-1:ny), ne_inner_upper_mod(-1:ny), 
     &       b_ne_inner_upper_mod(-1:ny), c_ne_inner_upper_mod(-1:ny),
     &       d_ne_inner_upper_mod(-1:ny),
     &      cy_te_inner_upper(-1:ny), te_inner_upper_mod(-1:ny),
     &       b_te_inner_upper_mod(-1:ny), c_te_inner_upper_mod(-1:ny),
     &       d_te_inner_upper_mod(-1:ny),
     &      cy_ti_inner_upper(-1:ny), ti_inner_upper_mod(-1:ny),
     &       b_ti_inner_upper_mod(-1:ny), c_ti_inner_upper_mod(-1:ny),
     &       d_ti_inner_upper_mod(-1:ny),
     &      cy_js_inner_upper(-1:ny), js_inner_upper_mod(-1:ny),
     &       b_js_inner_upper_mod(-1:ny), c_js_inner_upper_mod(-1:ny),
     &       d_js_inner_upper_mod(-1:ny)

*   .. for outer midplane
	 real (kind=R8) :: 
     &      cy_ne_outer_mdpl(-1:ny), ne_outer_mdpl_mod(-1:ny), 
     &       b_ne_outer_mdpl_mod(-1:ny), c_ne_outer_mdpl_mod(-1:ny),
     &       d_ne_outer_mdpl_mod(-1:ny),
     &      cy_te_outer_mdpl(-1:ny), te_outer_mdpl_mod(-1:ny),
     &       b_te_outer_mdpl_mod(-1:ny), c_te_outer_mdpl_mod(-1:ny),
     &       d_te_outer_mdpl_mod(-1:ny),
     &      cy_ti_outer_mdpl(-1:ny), ti_outer_mdpl_mod(-1:ny),
     &       b_ti_outer_mdpl_mod(-1:ny), c_ti_outer_mdpl_mod(-1:ny),
     &       d_ti_outer_mdpl_mod(-1:ny)

*   .. for outer lower
	 real (kind=R8) :: 
     &      cy_ne_outer_lower(-1:ny), ne_outer_lower_mod(-1:ny), 
     &       b_ne_outer_lower_mod(-1:ny), c_ne_outer_lower_mod(-1:ny),
     &       d_ne_outer_lower_mod(-1:ny),
     &      cy_te_outer_lower(-1:ny), te_outer_lower_mod(-1:ny),
     &       b_te_outer_lower_mod(-1:ny), c_te_outer_lower_mod(-1:ny),
     &       d_te_outer_lower_mod(-1:ny),
     &      cy_ti_outer_lower(-1:ny), ti_outer_lower_mod(-1:ny),
     &       b_ti_outer_lower_mod(-1:ny), c_ti_outer_lower_mod(-1:ny),
     &       d_ti_outer_lower_mod(-1:ny),
     &      cy_js_outer_lower(-1:ny), js_outer_lower_mod(-1:ny),
     &       b_js_outer_lower_mod(-1:ny), c_js_outer_lower_mod(-1:ny),
     &       d_js_outer_lower_mod(-1:ny)

*   .. for inner outer
	 real (kind=R8) :: 
     &      cy_ne_outer_upper(-1:ny), ne_outer_upper_mod(-1:ny), 
     &       b_ne_outer_upper_mod(-1:ny), c_ne_outer_upper_mod(-1:ny),
     &       d_ne_outer_upper_mod(-1:ny),
     &      cy_te_outer_upper(-1:ny), te_outer_upper_mod(-1:ny),
     &       b_te_outer_upper_mod(-1:ny), c_te_outer_upper_mod(-1:ny),
     &       d_te_outer_upper_mod(-1:ny),
     &      cy_ti_outer_upper(-1:ny), ti_outer_upper_mod(-1:ny),
     &       b_ti_outer_upper_mod(-1:ny), c_ti_outer_upper_mod(-1:ny),
     &       d_ti_outer_upper_mod(-1:ny),
     &      cy_js_outer_upper(-1:ny), js_outer_upper_mod(-1:ny),
     &       b_js_outer_upper_mod(-1:ny), c_js_outer_upper_mod(-1:ny),
     &       d_js_outer_upper_mod(-1:ny)


      end module b2mod_optim
