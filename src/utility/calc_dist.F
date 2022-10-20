      subroutine calc_dist(mpg,geo,clist,nc,iref,ds)
      ! Description
      !============
      ! Subroutine to compute relative distance between a set of cell 
      ! centers (Euclidean). The distance is computed between each 
      ! set of consecutive cells (i.e. x(iCv+1) - x(iCv) in 1D), but no
      ! check is made whether the cells are ordened in a 'logical' way. 

      ! Initialize
      !===========      
      ! Modules
      use b2mod_types
      use b2us_geo
      use b2us_map
      use b2mod_subsys
      implicit none
      ! Input variables
      type(geometry), intent(in) :: geo 
      type(mapping), intent(in) :: mpg 
      integer, intent(in) :: nc, iref, clist(1:nc)

      ! Output variables
      real (kind=R8), intent(out) :: ds(1:nc)
    
      ! Intermediate variables
      real (kind=R8) :: dsref
      integer :: iCv

      ! External functions
      external xertst

      ! Intrinsic functions
      intrinsic sqrt
    
      ! Main program
      !=============
      ! Call subini
      call subini ('output_ds_us')
      ! csc added safety for possible reference not within list number
      call xertst(iref.gt.1 .and. iref.le.nc,
     & 'output_ds_us: iref out of range')

      ! Compute the initial length as the distance of the first cell to 
      ! the first boundary face
      ds(1) = 0

      ! Loop over number of cells
      do iCv = 2,nc
        ! Compute cumulative distance
        ds(iCv) = ds(iCv-1) 
     &    + sqrt( (geo%cvX(clist(iCv)) - geo%cvX(clist(iCv-1)) )**2 
     &    + (geo%cvY(clist(iCv)) - geo%cvY(clist(iCv-1)) )**2  )
      enddo

      ! Compute reference point
      dsref = (ds(iref) + ds(iref-1))/2.0_R8

      ! Subtract reference length
      ds = ds - dsref
    
      ! Call subend and return
      call subend ()
      return
      end subroutine calc_dist
