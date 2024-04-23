      subroutine check_boundary_labels(m, nb, bstart, bend, docore)
      ! Description
      !============
      ! Subroutine to check that all boundary faces are assigned to
      ! a boundary condition or recycling strata.
      ! To do so, loop over all guard cells and get the face linked to
      ! each of them. Only one face per guard cell is possible!
      ! Then check that such face has a label different than zero and
      ! that is included in the lists bstart to bend
      ! Initialize
      !===========
      ! Modules
      use b2mod_types
      use b2us_map
      use b2mod_subsys
      implicit none
      ! Input variables
      type(mapping), intent(in) :: m
      integer, intent(in) :: nb, bstart(1:nb), bend(1:nb)
      logical, intent(in) :: docore

      ! Intermediate variables
      integer :: iCv, bnd_fc, ib, ll
      logical :: assigned

      ! External functions
      external xerrab

      ! Main computation
      !=============
      ! Call subini
      call subini ('check_boundary_labels')

      ! Loop over all guard cells
      do iCv = m%nCi+1, m%nCv
        ! get the only face linked to the guard cell
        bnd_fc = m%cvFc(m%cvFcP(iCv,1))
        ! if face is towards core, for recycling we don't care about it
        if (m%cvOnClosedSurface(m%fcCv(bnd_fc,1)) .or.
     &      m%cvOnClosedSurface(m%fcCv(bnd_fc,2)) .and.
     &      .not. docore) cycle
        if (bnd_fc.eq.0) then
          write(*,*) 'guard cell iCv=',iCv
          call xerrab('GUARD CELL HAS NO FACE FACING THE PLASMA!')
        endif
        if (m%fcLbl(bnd_fc).eq.0) then
          write(*,*) 'boundary face iFc=',bnd_fc
          call xerrab('BOUNDARY FACE HAS NO LABEL!')
        endif
        ! check that label is assigned somewhere
        assigned = .false.
        do ib = 1, nb
          do ll = bstart(ib), bend(ib)
            if (m%fcLbl(bnd_fc).eq.ll) then
              assigned = .true.
              exit
            endif
          enddo !ll
          if (assigned) exit
        enddo ! ib
        if (.not. assigned) then
          write(*,*) 'ifc, label', bnd_fc, m%fcLbl(bnd_fc)
          call xerrab('boundary face label not assigned to any BC/RC!')
        endif
      enddo !ic

      ! Call subend and return
      call subend ()
      return
      end subroutine check_boundary_labels
