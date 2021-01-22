      subroutine b2xzdd (nCv, nd, switch, mpg, src)
      use b2mod_switches
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv, nd
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   src(nCv, 0:nd)
*     ------------------------------------------------------------------
*      B2XZDD zeroes out the array elements corresponding to dead cells
*      and core boundary cells.
*     ------------------------------------------------------------------
      integer id, iCv
*     ------------------------------------------------------------------
      do iCv = 1, nCv
	 if(mpg%cvReg(iCv).eq.0.and.switch%zero_dead_and_core.ge.1) then
	   do id = 0, nd 
	     src(iCv,id) = 0.0_R8
	   enddo
	 elseif(iCv.gt.mpg%nCi .and. mpg%cvOnClosedSurface(iCv) .and.
     1	        switch%zero_dead_and_core.ge.2) then
	   do id = 0, nd
	     src(iCv,id) = 0.0_R8
	   enddo
	 endif
      enddo
      return
*     ------------------------------------------------------------------
      end
