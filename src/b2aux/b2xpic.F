      subroutine b2xpic (nCv, ns, qe, rza, bb, na, ua, zeff,       !srv 13.01.17
     &                   st_ext, fch_pi_c)
      use b2mod_types
      use b2us_plasma
      implicit none
      integer nCv, ns
      type (B2StateExt), intent (in) :: st_ext
      real (kind=R8) ::
     *  qe, rza(nCv,0:ns-1), na(nCv,0:ns-1), ua(nCv,0:ns-1),
     *  bb(nCv,0:3), zeff(nCv)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fch_pi_c(nCv,0:1)
*     ------------------------------------------------------------------
*     B2XPIC computes the poloidal projection of the ion parallel
*     current on cell centers.
*     Attention! There is no multiplication by vol/hx. 
*     It will be done later in B2XPIF.
*     ------------------------------------------------------------------
      integer is
      integer, save :: ncall = 0
      external sfill
*     ------------------------------------------------------------------
      fch_pi_c = 0.0_R8
      do is = 0, ns-1
       fch_pi_c(:,0) = fch_pi_c(:,0) +
     &  rza(:,is)*(rza(:,is)-zeff(:))*na(:,is)*ua(:,is)
      enddo
      do is = 0, st_ext%ns-1
       fch_pi_c(:,0) = fch_pi_c(:,0) +
     &  st_ext%za(:,is)*(st_ext%za(:,is)-zeff(:))*
     &  st_ext%na(:,is)*st_ext%ua(:,is)
      enddo
      fch_pi_c(:,0) = fch_pi_c(:,0)*qe*(bb(:,0)/bb(:,3))/zeff(:)

      ncall = ncall+1
      return
*     ------------------------------------------------------------------
      end subroutine b2xpic
