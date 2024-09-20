      subroutine b2nxdv (nCv, isb, switch, mpg, rob, smqb, smqdu, smb)
      use b2mod_types
      use b2mod_switches
      use b2us_map
#ifndef FIXED_POINT
      use b2mod_balance !djm Jan2017
     & , only : b2nxdv_smo0to3, balance_netcdf
#endif
      use b2mod_subsys
      implicit none
      integer nCv, isb
      type (switches), intent (in) :: switch
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *   rob(nCv), smqb(nCv,0:3), 
     &   smqdu(nCv), smb(nCv,0:3)
      real (kind=R8) :: smbtmp(nCv,0:3) !djm Jan2017
*     ------------------------------------------------------------------
*     B2NXDV computes the correction to the parallel momentum source
*     that is required to make the total friction force cancel.
*     ------------------------------------------------------------------
      integer iCv
      external xertst
*   ..initialisation
*     ------------------------------------------------------------------
      call subini ('b2nxdv')
*   ..test nCv
      call xertst (0.lt.nCv, 'faulty argument nCv')

cdjm Jan2017 Keep source at the start to calculate total change for balance routines
      smbtmp = smb

*   ..adjust smb
      if(switch%b2nxdv_style.eq.0) then
        smb(:,0) = smb(:,0) - smqb(:,1)*smqdu - smqb(:,3)*rob*smqdu
      else
        do iCv = 1, mpg%nCi
          smb(iCv,0) = smb(iCv,0)-smqb(iCv,1)*smqdu(iCv)-
     &                 smqb(iCv,3)*rob(iCv)*smqdu(iCv)
        enddo
      endif

#ifndef FIXED_POINT
cdjm Store the total change in linearised source for balance
      if (balance_netcdf.ne.0) then
        b2nxdv_smo0to3(1:nCv,0:3,isb) = smb-smbtmp
      endif
#endif

      call subend ()
      return
*     ------------------------------------------------------------------
      end subroutine b2nxdv
