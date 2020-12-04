      subroutine test_couple &
       & (nxdd_b2,nydd_b2,nstra_b2,nfl_b2,natm_b2,nmol_b2,nion_b2,nlimps_b2)
#ifdef B25_EIRENE
      use eirmod_parmmod
#endif
      use b2mod_subsys
      implicit none
      integer :: &
       & nxdd_b2, nydd_b2, nstra_b2, nfl_b2, natm_b2, nmol_b2, nion_b2, nlimps_b2
      logical error

      call subini ('test_couple')
      error=.false.
#ifdef B25_EIRENE
      if(nxdd_b2.ne.ndx) then
        write(*,*) 'NXDD(B2) <> NDX(EIRENE) ',nxdd_b2,ndx
        error=.true.
      endif
      if(nydd_b2.ne.ndy) then
        write(*,*) 'NYDD(B2) <> NDY(EIRENE) ',nydd_b2,ndy
        error=.true.
      endif
      if(nstra_b2.ne.nstra) then
        write(*,*) 'NSTRA(B2) <> NSTRA(EIRENE) ',nstra_b2,nstra
        error=.true.
      endif
      if(nfl_b2.ne.nfl) then
        write(*,*) 'NFL(B2) <> NFL(EIRENE) ',nfl_b2,nfl
        error=.true.
      endif
      if(natm_b2.ne.natm) then
        write(*,*) 'NATM(B2) <> NATM(EIRENE) ',natm_b2,natm
        error=.true.
      endif
      if(nmol_b2.ne.nmol) then
        write(*,*) 'NMOL(B2) <> NMOL(EIRENE) ',nmol_b2,nmol
        error=.true.
      endif
      if(nion_b2.ne.nion) then
        write(*,*) 'NION(B2) <> NION(EIRENE) ',nion_b2,nion
        error=.true.
      endif
      if(nlimps_b2.ne.nlimps) then
        write(*,*) 'NLIMPS(B2) <> NLIMPS(EIRENE) ',nlimps_b2,nlimps
        error=.true.
      endif
#endif
      if(error) then
        write(*,*) 'ERROR IN COUPLING PARAMETERS'
        write(*,*) 'PLEASE CHECK THAT BOTH EIRENE AND B2.5 ARE COMPILED'
        write(*,*) 'WITH THE SAME PARAMETERS FROM THE DIMENSIONS.F FILE'
        write(*,*) 'COUPLING TEST FAILED'
      else      
        write(*,*) 'COUPLING TEST PASSED'
      endif
      call xertst(.not.error, 'Coupling test failed!')
      call subend()
      return
      end subroutine test_couple



!!!Local Variables:
!!! mode: f90
!!! End:
