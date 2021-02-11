      subroutine b2nxdp (nCv, mpg, nregionv, solvereg, qe, te, nep, ne, 
     &  delpo)
      use b2mod_types
      use b2us_map
! csc The following are not necessary for computation but are needed
!     for adjoint AD to avoid side-effect variables
      use b2mod_ad
     & , only : ncall_b2nxdp
      implicit none
      integer nCv, nregionv
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     *  qe, te(nCv), nep(nCv), ne(nCv), delpo(nCv)
      logical solvereg(0:nregionv)
*     ------------------------------------------------------------------
*     B2NXDP computes the correction to the electric potential following
*     an update to the electron density.
*     ------------------------------------------------------------------
      integer iCv, ireg
      external subini, subend, xertst
      external b2xvsg
*   ..initialisation
*     ------------------------------------------------------------------
      call subini ('b2nxdp')
      call xertst (0.le.nCv, 'faulty argument nCv')
      call xertst (0.0_R8.lt.qe, 'faulty argument qe')
*   ..extensive tests on first few calls
      if (ncall_b2nxdp.lt.3) then
       call b2xvsg (nCv, te, 1, 'te', '.gt.')
       call b2xvsg (nCv, nep, 1, 'nep', '.gt.')
       call b2xvsg (nCv, ne, 1, 'ne', '.gt.')
      endif
*   ..compute delpo
      do iCv = 1, nCv
        ireg = mpg%cvReg(iCv)
        if (solvereg(ireg)) then
         delpo(iCv) = (te(iCv)/qe)*
     &    (ne(iCv)-nep(iCv))*2.0_R8/(ne(iCv)+nep(iCv))
        else
         delpo(iCv) = 0.0_R8
        endif
      enddo
      ncall_b2nxdp = ncall_b2nxdp+1
      call subend ()
      return
*     ------------------------------------------------------------------
      end
