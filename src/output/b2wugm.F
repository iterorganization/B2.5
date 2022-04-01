*-----------------------------------------------------------------------
*.specification

      subroutine b2wugm (nput, nx, ny)
      use b2mod_types
      use b2mod_geo
      use b2mod_indirect
      use b2mod_version
      use b2mod_b2cmgs
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer nput, nx, ny
*   ..output arguments (unspecified on entry)
*     (none)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WUGM writes to an un*formatted file the geometry and magnetic
*     field.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine writes to unit (nput) the specification of a B2
*     geometry and magnetic field. Unit nput must be connected to an
*     un*formatted file.
*     The data consists of the coordinates of the four corners of all
*     the cells in the mesh, the magnetic flux and field functions at
*     those positions, the components of the magnetic field at cell
*     centers, a collection of metric coefficients, and a specification
*     of singular regions in the geometry. See the code to discern the
*     format of the output file.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nput - integer, input.
*     nput specifies a unit number for un*formatted output.
*
*     nx, ny - integer, input.
*     See the calling routine for description.
*
*     crx, cry, fpsi, ffbz - real*8 array, input.
*     bb, vol, hx, hy, qz, qc, gs, pbs - real*8 array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmgs/ - input.
*     /b2cmgs/ specifies the singular regions in the geometry.
*     See the calling routine for description.
*
*
*  8. further comments
*
*     See routine b2cdca for explanation of 'un*formatted file'.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer n2, idum(0:9), ix, iy, ixl, ixr, iyb, iyt, errors
      integer nxiso1, nxiso2, nyiso1, nyiso2
      logical done
*   ..initialisation
      data nxiso1/-2/, nxiso2/-2/, nyiso1/-2/, nyiso2/-2/
*   ..procedures
      external xertst, cfwure, cfwuin, ipgeti
      external b2xvsg, b2xvgm
      real (kind=R8) ::
     *   distval

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wugm')
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nput, nx, ny
      call xertst (1.le.nput, 'faulty argument nput')
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
*   ..test sign of bb
      call b2xvsg (n2, bb(-1,-1,3), 1, 'bb3', '.gt.')
*   ..test geometry
      call b2xvgm (nx, ny, vol, hx, hy, qz, qc, qs, gs, pbs,
     &  cell_width, cell_height)
*   ..obtain isolated region data
      call ipgeti('b2agdr_nxiso1', nxiso1)
      call ipgeti('b2agdr_nxiso2', nxiso2)
      call ipgeti('b2agdr_nyiso1', nyiso1)
      call ipgeti('b2agdr_nyiso2', nyiso2)
      if (nxiso1.ne.-2) then
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso1.eq.-2) nyiso1 = -1
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nxiso2.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nyiso1.eq.-2) nyiso1 = -1
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nyiso1.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso2.eq.-2) nyiso2 = ny
      endif
      if (nyiso2.ne.-2) then
        if (nxiso1.eq.-2) nxiso1 = -1
        if (nxiso2.eq.-2) nxiso2 = nx
        if (nyiso1.eq.-2) nyiso1 = -1
      endif

cdpc
      write(*,*) 'Checking CRX'
      errors=0
      do iy=-1,ny
c check right neighbours
        do ix=-1,nx-1
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ixr=rightix(ix,iy)
          if(.not.isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            write(*,*) 'No right neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. ix.eq.nxiso2 .and.
     &     iy.ge.nyiso1 .and. iy.le.nyiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. ixr.eq.ix+1) then
            distval=
     1       abs_diff(crx(ix,iy,1),crx(ixr,rightiy(ix,iy),0))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1         'Error with right neighbour of ',
     2          ix,iy,' (lower point), ','dist=',distval
            endif
            distval=
     1       abs_diff(crx(ix,iy,3),crx(ixr,rightiy(ix,iy),2))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1         'Error with right neighbour of ',
     2          ix,iy,' (upper point), ','dist=',distval
            endif
          endif
        enddo
c check left neighbours
        do ix=0,nx
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ixl=leftix(ix,iy)
          if(.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
            write(*,*) 'No left neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. ix.eq.nxiso1 .and.
     &     iy.ge.nyiso1 .and. iy.le.nyiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. ixl.eq.ix-1) then
            distval=
     1       abs_diff(crx(ix,iy,0),crx(ixl,leftiy(ix,iy),1))
            if(distval.ne.0) then
            if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1         'Error with left neighbour of ',ix,iy,' (lower point), ',
     2         'dist=',distval
            endif
            distval=
     1       abs_diff(crx(ix,iy,2),crx(ixl,leftiy(ix,iy),3))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1         'Error with left neighbour of ',ix,iy,' (upper point), ',
     2         'dist=',distval
            endif
          endif
        enddo
      enddo
      do ix=-1,nx
c check top neighbours
        do iy=-1,ny-1
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          iyt=topiy(ix,iy)
          if(.not.isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
            write(*,*) 'No top neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. iy.eq.nyiso2 .and.
     &     ix.ge.nxiso1 .and. ix.le.nxiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. iyt.eq.iy+1) then
            distval=
     1       abs_diff(crx(ix,iy,2),crx(topix(ix,iy),iyt,0))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with top neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(crx(ix,iy,3),crx(topix(ix,iy),iyt,1))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with top neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
c check bottom neighbours
        do iy=0,ny
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          iyb=bottomiy(ix,iy)
          if(.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)))
     >     then
            write(*,*) 'No bottom neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. iy.eq.nyiso1 .and.
     &     ix.ge.nxiso1 .and. ix.le.nxiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. iyb.eq.iy-1) then
            distval=
     1       abs_diff(crx(ix,iy,0),crx(bottomix(ix,iy),iyb,2))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with bottom neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(crx(ix,iy,1),crx(bottomix(ix,iy),iyb,3))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with bottom neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
      write(*,*) 'Checking CRY'
      do iy=-1,ny
c check right neighbours
        do ix=-1,nx-1
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ixr=rightix(ix,iy)
          if(.not.isInDomain(nx,ny,rightix(ix,iy),rightiy(ix,iy))) then
            write(*,*) 'No right neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. ix.eq.nxiso2 .and.
     &     iy.ge.nyiso1 .and. iy.le.nyiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. ixr.eq.ix+1) then
            distval=
     1       abs_diff(cry(ix,iy,1),cry(ixr,rightiy(ix,iy),0))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with right neighbour of ',ix,iy,' (lower point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(cry(ix,iy,3),cry(ixr,rightiy(ix,iy),2))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with right neighbour of ',ix,iy,' (upper point), ',
     2       'dist=',distval
           endif
          endif
        enddo
c check left neighbours
        do ix=0,nx
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          ixl=leftix(ix,iy)
          if(.not.isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
            write(*,*) 'No left neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. ix.eq.nxiso1 .and.
     &     iy.ge.nyiso1 .and. iy.le.nyiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. ixl.eq.ix-1) then
            distval=
     1       abs_diff(cry(ix,iy,0),cry(ixl,leftiy(ix,iy),1))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with left neighbour of ',ix,iy,' (lower point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(cry(ix,iy,2),cry(ixl,leftiy(ix,iy),3))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with left neighbour of ',ix,iy,' (upper point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
      do ix=-1,nx
c check top neighbours
        do iy=-1,ny-1
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          iyt=topiy(ix,iy)
          if(.not.isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
            write(*,*) 'No top neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. iy.eq.nyiso2 .and.
     &     ix.ge.nxiso1 .and. ix.le.nxiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. iyt.eq.iy+1) then
            distval=
     1       abs_diff(cry(ix,iy,2),cry(topix(ix,iy),iyt,0))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with top neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(cry(ix,iy,3),cry(topix(ix,iy),iyt,1))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with top neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
c check bottom neighbours
        do iy=0,ny
          if(isUnusedCell(cflags(ix,iy,CELLFLAG_TYPE))) cycle
          iyb=bottomiy(ix,iy)
          if(.not.isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy)))
     >     then
            write(*,*) 'No bottom neighbour for ',ix,iy
          elseif(region(ix,iy,0).eq.0 .and. iy.eq.nyiso1 .and.
     &     ix.ge.nxiso1 .and. ix.le.nxiso2) then
            write(*,*) 'Testing skipped for isolated region edge at ',
     .      ix,iy
          elseif(periodic_bc.eq.0 .or. iyb.eq.iy-1) then
            distval=
     1       abs_diff(cry(ix,iy,0),cry(bottomix(ix,iy),iyb,2))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with bottom neighbour of ',ix,iy,' (left point), ',
     2       'dist=',distval
            endif
            distval=
     1       abs_diff(cry(ix,iy,1),cry(bottomix(ix,iy),iyb,3))
            if(distval.ne.0) then
              if(distval.gt.1e-6) errors=errors+1
              write(*,*)
     1       'Error with bottom neighbour of ',ix,iy,' (right point), ',
     2       'dist=',distval
            endif
          endif
        enddo
      enddo
cxpb
      done = .false.
      do ix = -1, nx
        do iy = -1, ny
          if(region(ix,iy,0).eq.0) then
            if (.not.done) then
              write(*,*) 'Checking neighbours of dead cells'
              done = .true.
            endif
            if(isInDomain(nx,ny,leftix(ix,iy),leftiy(ix,iy))) then
              if (region(leftix(ix,iy),leftiy(ix,iy),0).ne.0) then
                errors=errors+1
                write(*,*) 'Left neighbour of dead cell ',ix,iy,
     $                     ' alive!'
                write(*,*) 'ixl,iyl',leftix(ix,iy),leftiy(ix,iy)
              endif
            endif
            if(isInDomain(nx,ny,rightix(ix,iy),leftiy(ix,iy))) then
              if (region(rightix(ix,iy),rightiy(ix,iy),0).ne.0) then
                errors=errors+1
                write(*,*) 'Right neighbour of dead cell ',ix,iy,
     $                     ' alive!'
                write(*,*) 'ixr,iyr',rightix(ix,iy),rightiy(ix,iy)
              endif
            endif
            if(isInDomain(nx,ny,bottomix(ix,iy),bottomiy(ix,iy))) then
              if (region(bottomix(ix,iy),bottomiy(ix,iy),0).ne.0) then
                errors=errors+1
                write(*,*) 'Bottom neighbour of dead cell ',ix,iy,
     .                     ' alive!'
                write(*,*) 'ixb,iyb',bottomix(ix,iy),bottomiy(ix,iy)
              endif
            endif
            if(isInDomain(nx,ny,topix(ix,iy),topiy(ix,iy))) then
              if (region(topix(ix,iy),topiy(ix,iy),0).ne.0) then
                errors=errors+1
                write(*,*) 'Top neighbour of dead cell ',ix,iy,
     $                     ' alive!'
                write(*,*) 'ixt,iyt',topix(ix,iy),topiy(ix,iy)
              endif
            endif
          endif
        enddo
      enddo
cxpb
      if(errors.eq.0) then
        write(*,*) 'No errors found in geometry'
      elseif(errors.eq.1) then
        write(*,*) errors,' error found in the geometry'
      else
        write(*,*) errors,' errors found in the geometry'
      endif
      call xertst (errors.eq.0, 'Errors found in geometry')
cdpc-end

* ..save geometry
*   ..save coordinates
      call cfwure (nput, n2*4, crx, 'crx')
      call cfwure (nput, n2*4, cry, 'cry')
*   ..save flux functions
      call cfwure (nput, n2*4, fpsi, 'fpsi')
      call cfwure (nput, n2*4, ffbz, 'ffbz')
*   ..save magnetic field
      call cfwure (nput, n2*4, bb, 'bb')
*   ..save vol, hx, hy
      call cfwure (nput, n2, vol, 'vol')
      call cfwure (nput, n2, hx, 'hx')
      call cfwure (nput, n2, hy, 'hy')
*   ..save qz, qc, gs
      call cfwure (nput, n2*2, qz, 'qz')
cxpb
      if(newversion.ge.'03.001.000') then
        call cfwure (nput, n2*2, qc, 'qc')
        call cfwure (nput, n2*2, qs, 'qs')
      elseif(newversion.ge.'01.001.022') then
        call cfwure (nput, n2, qc, 'qc')
      endif
cxpb-end
      call cfwure (nput, n2*3, gs, 'gs')
*   ..save /b2cmgs/
      idum(0) = nlreg
      call cfwuin (nput, 1, idum, 'nlreg')
      call cfwuin (nput, nlreg, nlxlo, 'nlxlo')
      call cfwuin (nput, nlreg, nlxhi, 'nlxhi')
      call cfwuin (nput, nlreg, nlylo, 'nlylo')
      call cfwuin (nput, nlreg, nlyhi, 'nlyhi')
      call cfwuin (nput, nlreg, nlloc, 'nlloc')
cdpc
      idum(0)=nncut
      call cfwuin (nput, 1, idum, 'nncut')
      call cfwuin (nput, nncut, leftcut, 'leftcut')
      call cfwuin (nput, nncut, rightcut, 'rightcut')
      call cfwuin (nput, nncut, topcut, 'topcut')
      call cfwuin (nput, nncut, bottomcut, 'bottomcut')
      call cfwuin (nput, n2, leftix, 'leftix')
      call cfwuin (nput, n2, rightix, 'rightix')
      call cfwuin (nput, n2, topix, 'topix')
      call cfwuin (nput, n2, bottomix, 'bottomix')
      call cfwuin (nput, n2, leftiy, 'leftiy')
      call cfwuin (nput, n2, rightiy, 'rightiy')
      call cfwuin (nput, n2, topiy, 'topiy')
      call cfwuin (nput, n2, bottomiy, 'bottomiy')
      call cfwuin (nput, n2*3, region, 'region')
      call cfwuin (nput, 3, nnreg, 'nnreg')
      call cfwuin (nput, n2*2, resignore, 'resignore')
cdpc-end
cxpb
      if(newversion.ge.'01.001.000') then
        idum(0)=periodic_bc
        call cfwuin (nput, 1, idum, 'periodic_bc')
      endif
cxpb-end
      if(newversion.ge.'01.001.063') then          !srv 15.02.11
        call cfwure (nput, n2*2, pbs, 'pbs')
      endif
      if(newversion.ge.'03.000.006') then
        idum(0)=redef_pbs
        call cfwuin (nput, 1, idum, 'redef_pbs')
      endif
* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wugm

      end subroutine b2wugm
