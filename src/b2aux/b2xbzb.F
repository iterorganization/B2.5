      subroutine b2xbzb (nCv, nFc, nVx, mpg, bb, bbf, bbv, vol,
     .  bzb, bzbf, bzbv)                                                !srv 13.10.06 23.07.10 {
      use b2mod_types
      use b2us_map
      implicit none
      integer nCv, nFc, nVx
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   bb(nCv,0:3), bbf(nFc,0:3), bbv(nVx,0:3), vol(nCv),
     *   bzb(nCv), bzbf(nFc), bzbv(nVx)
*     ------------------------------------------------------------------
*     B2XBZB calculates Bz/B**2 - Bz/<B**2>.
*     ------------------------------------------------------------------
      integer iCv
      real (kind=R8) ::
     *   bbsum, volsum, bzb0
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
        bzb0 = bbsum/volsum
        bzb  = bb(:,2)/bb(:,3)**2 - bb(:,2)/bzb0
        bzbf = bbf(:,2)/bbf(:,3)**2 - bbf(:,2)/bzb0
        bzbv = bbv(:,2)/bbv(:,3)**2 - bbv(:,2)/bzb0
      else
        bzb  = 0.0_R8
        bzbf = 0.0_R8
        bzbv = 0.0_R8
      endif

*$$$  call subend ()
      return
*     ------------------------------------------------------------------
      end                                                               !srv 13.10.06 23.07.10 }
