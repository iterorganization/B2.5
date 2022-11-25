      subroutine calc_dist_f(mpg,geo,flist,nf,iref,ds)
      ! Description
      !============
      ! Subroutine to compute distance along a set of consecutive faces.
      ! No check is made whether the cells are ordened in a 'logical' way. 

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
      integer, intent(in) :: nf, iref, flist(1:nf)

      ! Output variables
      real (kind=R8), intent(out) :: ds(1:nf)
    
      ! Intermediate variables
      real (kind=R8) :: dsref
      integer :: iFc

      ! External functions
      external xertst

      ! Intrinsic functions
      intrinsic sqrt
    
      ! Main program
      !=============
      ! Call subini
      call subini ('calc_dist_f')
      ! Safety for possible reference not within list number
      call xertst(iref.ge.0 .and. iref.le.nf,
     & 'calc_dist_f: iref out of range')

      ! Compute the initial length as half the length along the first face
      ds(1) = 0.5_R8*geo%fcHt(flist(1))

      ! Loop over number of faces
      do iFc = 2,nf
        ! Compute cumulative distance
        ds(iFc) = ds(iFc-1) + 
     &    0.5_R8*(geo%fcHt(flist(iFc)) + geo%fcHt(flist(iFc-1)))
      enddo

      ! If reference face provided, shift length computation
      if (iref.gt.1) then
        ! Compute reference point
        dsref = (ds(iref) + ds(iref-1))/2.0_R8

        ! Subtract reference length
        ds = ds - dsref
      endif

      ! Call subend and return
      call subend ()
      return
      end subroutine calc_dist_f
