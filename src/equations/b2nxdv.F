      subroutine b2nxdv (nCv, isb, mpg, rob, smqb, smqdu, smb)
      use b2mod_types
      use b2us_map
!WG_TODO      use b2mod_balance !djm Jan2017
!WG_TODO     & , only : b2nxdv_smo0to3, balance_netcdf
      implicit none
      integer nCv, isb
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   rob(nCv), smqb(nCv,0:3), 
     &   smqdu(nCv), smb(nCv,0:3)
!WG_TODO      real (kind=R8) :: smbtmp(nCv,0:3) !djm Jan2017
*     ------------------------------------------------------------------
*     B2NXDV computes the correction to the parallel momentum source
*     that is required to make the total friction force cancel.
*     ------------------------------------------------------------------
      integer ncall, iCv
      integer, save :: style = 1                                       !srv 26.10.08
      save ncall
      external subini, subend, xertst
*   ..initialisation
      data ncall/0/
*     ------------------------------------------------------------------
      call subini ('b2nxdv')
*   ..set internal parameters on first call
      if (ncall.eq.0) then
       call ipgeti ('b2nxdv_style', style)                            !srv 26.10.08
      endif
*   ..test nCv
      call xertst (0.le.nCv, 'faulty argument nCv')

!WG_TODOcdjm Jan2017 Keep source at the start to calculate total change for balance routines
!WG_TODO      smbtmp = smb

*   ..adjust smb
      if(style.eq.0) then
        smb(:,0) = smb(:,0) - smqb(:,1)*smqdu - smqb(:,3)*rob*smqdu
      else
        do iCv = 1, mpg%nCi
          smb(iCv,0) = smb(iCv,0)-smqb(iCv,1)*smqdu(iCv)-
     &                 smqb(iCv,3)*rob(iCv)*smqdu(iCv)
        enddo
      endif

cdjm Store the total change in linearised source for balance
!WG_TODO      if (balance_netcdf.ne.0) then
!WG_TODO        b2nxdv_smo0to3(-1:nx,-1:ny,0:3,isb) = smb-smbtmp
!WG_TODO      endif

      ncall = ncall+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
