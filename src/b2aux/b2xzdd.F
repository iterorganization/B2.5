      subroutine b2xzdd (nCv, nd, mpg, src)
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv, nd
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   src(nCv, 0:nd)
*     ------------------------------------------------------------------
*      B2XZDD zeroes out the array elements corresponding to dead cells
*      and core boundary cells.
*     ------------------------------------------------------------------
      integer id, iCv
      integer ncall, zero_dead_and_core
      save ncall, zero_dead_and_core
*   ..procedures
      external ipgeti
*   ..initialisation
      data ncall/0/, zero_dead_and_core/1/
*     ------------------------------------------------------------------

      if (ncall.eq.0) then
	call ipgeti('b2xzdd_zero_dead_and_core', zero_dead_and_core)
      endif
      do iCv = 1, nCv
	 if(mpg%cvReg(iCv).eq.0.and.zero_dead_and_core.ge.1) then
	   do id = 0, nd 
	     src(iCv,id) = 0.0_R8
	   enddo
	 elseif(iCv.gt.mpg%nCi .and. mpg%cvOnClosedSurface(iCv) .and.
     1	        zero_dead_and_core.ge.2) then
	   do id = 0, nd
	     src(iCv,id) = 0.0_R8
	   enddo
	 endif
      enddo
      ncall=ncall+1
      return
*     ------------------------------------------------------------------
      end
