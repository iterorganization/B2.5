*-----------------------------------------------------------------------
*.specification

      subroutine b2wgpj_gnuplot (nprt, nx, ny, ns, 
     &  na, ua, te, ti, po, fna, fhe, fhi, fch, 
     &  ne, fne, kinrgy, gs, qc, ni)
      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_subsys
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif
      implicit none
*   ..input arguments (unchanged on exit)
      integer nprt, nx, ny, ns
      real (kind=R8) ::
     *  na(-1:nx,-1:ny,0:ns-1), ua(-1:nx,-1:ny,0:ns-1),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), po(-1:nx,-1:ny),
     &  fna(-1:nx,-1:ny,0:1,0:1,0:ns-1), fhe(-1:nx,-1:ny,0:1,0:1),
     &  fhi(-1:nx,-1:ny,0:1,0:1), fch(-1:nx,-1:ny,0:1,0:1),
     &  ne(-1:nx,-1:ny), fne(-1:nx,-1:ny,0:1,0:1),
     &  kinrgy(-1:nx,-1:ny,0:ns-1), ni(-1:nx,-1:ny,0:1),
     &  gs(-1:nx,-1:ny, 0:2), qc(-1:nx,-1:ny,0:1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2WGPJ provides graphical output of the plasma state on a uniform
*     cartesian mesh.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine produces plots of the plasma state that results from
*     a B2 calculation. Arrays na, ua, te, ti, and po provide the state
*     of the plasma. Arrays fna, fhe, fhi, and fch provide the particle
*     flux, electron heat flux, atom heat flux, and the electric
*     current. This routine does not require any metric information.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nprt - integer, input.
*     nprt specifies a unit number for formatted output.
*
*     nx, ny, ns - integer, input.
*     See the calling routine for description.
*
*     na, ua, te, ti, po, fna, fhe, fhi, fch - real*8 array, input.
*     ne, fne - real*8 array, input.
*     See the calling routine for description.
*
*     (parameters in common:)
*
*     /b2cmpa/ - input.
*     /b2cmpa/ specifies basic physics parameters.
*     See the calling routine for description.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer nmax, ix, iy, islice
      parameter (nmax=300)
      character*4 want,slice
*   ..procedures
      intrinsic max
#ifndef NAGF95
#ifndef NAGFOR
      integer iargc
#else
      integer lenval, ierror
#endif
#endif
      external xertst

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2wgpj_gnuplot')
*   ..test nprt
      call xertst (1.le.nprt, 'faulty argument nprt')
*   ..test nx, ny, ns
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
*   ..test if sufficient workspace
      call xertst (nx.le.nmax.and.ny.le.nmax,
     &  'faulty parameter nmax')

* ..make two-dimensional plots
      islice=1
#ifndef NAGFOR
      if(iargc().le.0) then
#else
      if(command_argument_count().le.0) then
#endif
         want='te  '
      else
#ifndef NAGFOR
         call getarg(1,want)
         if(iargc().ge.2) then
            call getarg(2,slice)
#else
         call get_command_argument(1,want,lenval,ierror)
         if(command_argument_count().ge.2) then
            call get_command_argument(2,slice,lenval,ierror)
#endif
            read(slice,*,err=5) islice
            islice=min(max(islice,0),ns-1)
         endif
      endif
    5 continue
      write(*,'(a,a,i4)') '#',want,islice
      if(want.eq.'te  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,te(ix,iy)/ev
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fhex') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fhe(ix,iy,0,0)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fhey') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fhe(ix,iy,1,1)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'ti  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,ti(ix,iy)/ev
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fhix') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fhi(ix,iy,0,0)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fhiy') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fhi(ix,iy,1,1)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'ne  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,ne(ix,iy)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fnex') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fne(ix,iy,0,0)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fney') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fne(ix,iy,1,1)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'po  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,po(ix,iy)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fchx') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fch(ix,iy,0,0)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fchy') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,fch(ix,iy,1,1)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'na  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,na(ix,iy,islice)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fnax') then
         do iy=-1,ny
            do ix=-1,nx
csw 17feb2012 fixed             write(*,*) ix,iy,fna(ix,iy,islice,0)
               write(*,*) ix,iy,fna(ix,iy,0,0,islice)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'fnay') then
         do iy=-1,ny
            do ix=-1,nx
csw 17feb2012 fixed              write(*,*) ix,iy,fna(ix,iy,islice,1)
               write(*,*) ix,iy,fna(ix,iy,1,1,islice)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'ua  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,ua(ix,iy,islice)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'kinr') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,kinrgy(ix,iy,islice)
            enddo
            write(*,'()')
         enddo
csw 30apr2014 added
      elseif(want.eq.'sx') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,gs(ix,iy,0)
            enddo
            write(*,'()')
         enddo
      elseif(want.eq.'sxpe') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,gs(ix,iy,0)*qc(ix,iy,0)
            enddo
            write(*,'()')
         enddo
csw 17sep2014 added
      elseif(want.eq.'ni  ') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,ni(ix,iy,0)
            enddo
            write(*,'()')
         enddo
csw 30sep2014 added
      elseif(want.eq.'sy') then
         do iy=-1,ny
            do ix=-1,nx
               write(*,*) ix,iy,gs(ix,iy,1)
            enddo
            write(*,'()')
         enddo
      else
         write(*,*) 'Have not coded ',want,' yet'
      endif

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2wgpj_gnuplot

      end subroutine b2wgpj_gnuplot
