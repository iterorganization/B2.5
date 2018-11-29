      subroutine b2xbzb (nCv, mpg, bb, vol, bzb)              !srv 13.10.06 23.07.10 {
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   bb(nCv,0:3), vol(nCv), bzb(nCv)
*     ------------------------------------------------------------------
*     B2XBZB calculates Bz/B**2 - Bz/<B**2>.
*     ------------------------------------------------------------------
      integer iCv
      real (kind=R8) ::
     *   bbsum, volsum
*     ------------------------------------------------------------------
*$$$  call subini ('b2xbzb')

      bbsum  = 0.0e0_R8
      volsum = 0.0e0_R8
      do iCv = 1, nCv
        if(mpg%cvOnClosedSurface(iCv)) then      ! for CORE
          bbsum  = bbsum  + bb(iCv,3)**2*vol(iCv)
          volsum = volsum + vol(iCv)
        endif
      enddo
      if(volsum.gt.0.0_R8) then
        bzb = bbsum/volsum
        bzb = bb(:,2)/bb(:,3)**2 - bb(:,2)/bzb
      else
        bzb = 0.0_R8
      endif

*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                   !srv 13.10.06 23.07.10 }
