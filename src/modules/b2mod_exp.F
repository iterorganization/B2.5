      module b2mod_exp
       use b2mod_types
*   .. for inner lower
       integer :: 
     &           ndat_ne_inner_lower_exp, ndat_te_inner_lower_exp,
     &           ndat_ti_inner_lower_exp, ndat_js_inner_lower_exp

       real (kind=R8), allocatable, save :: 
     &   cy_ne_inner_lower_exp(:), ne_inner_lower_exp(:), 
     &   cy_te_inner_lower_exp(:), te_inner_lower_exp(:),
     &   cy_ti_inner_lower_exp(:), ti_inner_lower_exp(:),
     &   cy_js_inner_lower_exp(:), js_inner_lower_exp(:)

*   .. for inner upper
       integer :: 
     &           ndat_ne_inner_upper_exp, ndat_te_inner_upper_exp,
     &           ndat_ti_inner_upper_exp, ndat_js_inner_upper_exp

       real (kind=R8), allocatable, save :: 
     &   cy_ne_inner_upper_exp(:), ne_inner_upper_exp(:), 
     &   cy_te_inner_upper_exp(:), te_inner_upper_exp(:),
     &   cy_ti_inner_upper_exp(:), ti_inner_upper_exp(:),
     &   cy_js_inner_upper_exp(:), js_inner_upper_exp(:)

*   .. for outer midplane
       integer :: 
     &           ndat_ne_outer_mdpl_exp, ndat_te_outer_mdpl_exp,
     &           ndat_ti_outer_mdpl_exp

       real (kind=R8), allocatable, save :: 
     &   cy_ne_outer_mdpl_exp(:), ne_outer_mdpl_exp(:), 
     &   cy_te_outer_mdpl_exp(:), te_outer_mdpl_exp(:),
     &   cy_ti_outer_mdpl_exp(:), ti_outer_mdpl_exp(:)

*   .. for outer lower
       integer :: 
     &           ndat_ne_outer_lower_exp, ndat_te_outer_lower_exp,
     &           ndat_ti_outer_lower_exp, ndat_js_outer_lower_exp

       real (kind=R8), allocatable, save :: 
     &   cy_ne_outer_lower_exp(:), ne_outer_lower_exp(:), 
     &   cy_te_outer_lower_exp(:), te_outer_lower_exp(:),
     &   cy_ti_outer_lower_exp(:), ti_outer_lower_exp(:),
     &   cy_js_outer_lower_exp(:), js_outer_lower_exp(:)

*   .. for outer upper
       integer :: 
     &           ndat_ne_outer_upper_exp, ndat_te_outer_upper_exp,
     &           ndat_ti_outer_upper_exp, ndat_js_outer_upper_exp

       real (kind=R8), allocatable, save :: 
     &   cy_ne_outer_upper_exp(:), ne_outer_upper_exp(:), 
     &   cy_te_outer_upper_exp(:), te_outer_upper_exp(:),
     &   cy_ti_outer_upper_exp(:), ti_outer_upper_exp(:),
     &   cy_js_outer_upper_exp(:), js_outer_upper_exp(:)
       real (kind=R8) :: power_exp

      end module b2mod_exp
