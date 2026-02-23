      subroutine b2urbp(nx, ny, ns, rza, pbs, po, te, ti, ne, na,
     &                  fch, sch, nmb, conv)
      use b2mod_types
      use b2mod_indirect
      use b2mod_constants
      use b2mod_b2cmpa
      use b2mod_boundary_namelist
      use b2mod_geo, only : crdy_c_cell
      use b2mod_math
      implicit none
*   ..input arguments (unchanged on exit)
      integer :: nx, ny, ns
      real (kind=R8) :: po(-1:nx,-1:ny),
     &  te(-1:nx,-1:ny), ti(-1:nx,-1:ny), ne(-1:nx,-1:ny)
      real (kind=R8) :: pbs(-1:nx,-1:ny,0:1), fch(-1:nx,-1:ny,0:1),
     &  na(-1:nx,-1:ny,0:ns-1), sch(-1:nx,-1:ny,0:3)
      real (kind=R8) :: rza(-1:nx,-1:ny,0:ns-1)
      logical :: nmb
*     if nmb = .true. then output the cell number
*                otherwise output the radial coordinate
*   ..output arguments (unspecified on entry)
      logical :: conv
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2URBP checks the convergence of the electric potential at
*     the targets using relative (EPSR) and/or absolute (EPSA) criteria.
*     For diagnostic purposes, B2URBP calculates and outputs
*     alternative measures of residuals at the targets if the b2urbp_iout
*     switch is nonzero.
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer :: ix, iy, is, iw
      integer, parameter :: maxd = 14
      real (kind=R8) :: pz(-1:nx,-1:ny), rz(-1:nx,-1:ny)
      real (kind=R8) :: s1, vte, cs, t0, fchi, fche, fch_real_bc,
     &     wdim(maxd)
      real (kind=R8) :: fch_code, seec,
     &     res, dif_real_bc, res_bc, fch_sch
      real (kind=R8),save :: epsr=1.0e-2_R8, epsa=1.0e-2_R8
      character(7), save :: folder
      integer, save :: iappend = 0
      integer, save :: ncall = 0, iout = 0

      conv=.true.
      if (ncall.eq.0) then
        call ipgeti('b2wdat_append', iappend)
!!!   ..try to create a file in optional output folder
        open(70,file='output/__1__.dat',err=10)
!!!   ..folder exists: future output will be done in folder
        folder='output/'
        write(70,*) 'Test file to check on output folder'
        close(70)
        go to 20
!!!   ..folder does not exist: future output will be done in current folder
   10   folder=''
   20   continue
        call ipgetr('b2urbp_epsr', epsr)
        call ipgetr('b2urbp_epsa', epsa)
        call ipgeti('b2urbp_iout', iout)
      endif

      call b2xppz (nx, ny, ns, ne, na, te, ti, pz)
      call b2xprz (nx, ny, ns, mp, am, na, rz)

      if(iout.ne.0) then
        if(LSN) then
          call b2wdat_open(70,'res_at_inner_plate_bcpot=11')
        else
          call b2wdat_open(70,'res_at_outer_plate_bcpot=11')
        endif
      endif

      do iy = -1, ny
        iw = 0
        do ix = -1,nx
          if(leftix(ix,iy).eq.-2) then
            if(iout.ne.0) then
              s1 = abs(pbs(rightix(ix,iy),rightiy(ix,iy),0))
              vte = sqrt(te(ix,iy)/me)
              CS = SQRT(PZ(ix,iy)/RZ(ix,iy))
              t0 = 0.0_R8
              do is = 0, ns-1
                t0 = t0+rza(ix,iy,is)*cs*s1*na(ix,iy,is)
              enddo
              fchi = qe*t0
              fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*ne(ix,iy)*vte*s1*
     &               exp(-qe*po(ix,iy)/te(ix,iy))
              fch_real_bc = (1.0_R8-gammae)*fche-fchi
              res = fch_real_bc - fch(rightix(ix,iy),rightiy(ix,iy),0)
            endif

            fch_sch = sch(ix,iy,0)+sch(ix,iy,1)*po(ix,iy)+
     &                sch(ix,iy,2)*ne(ix,iy)+
     &                sch(ix,iy,3)*po(ix,iy)*ne(ix,iy)
            res_bc = fch_sch - fch(rightix(ix,iy),rightiy(ix,iy),0)
            if(abs(res_bc).ge.
     .         abs(fch(rightix(ix,iy),rightiy(ix,iy),0))*epsr+epsa)
     >       conv=.false.

            if(iout.ne.0) then
              seec = gammae
              t0 = max(fchi,(1.0_R8-seec)*fche)*qe/te(ix,iy)
              fch_code = (1.0_R8-seec)*fche-fchi
              dif_real_bc = fch_real_bc-fch_sch
              if(nmb) then
                wdim(iw+1) = iy
              else
                wdim(iw+1) = crdy_c_cell(ix,iy)
              endif
              wdim(iw+2) = fch_real_bc
              wdim(iw+3) = fch(rightix(ix,iy),rightiy(ix,iy),0)
              wdim(iw+4) = res
              wdim(iw+5) = fch_code
              wdim(iw+6) = res_bc !  dif
              wdim(iw+7) = dif_real_bc

              iw=iw+7
              call xertst(iw.le.maxd,
     .         'index exceeded maxd=14 in b2urbp.F')
            endif
          endif
        enddo
        if(iout.ne.0) write(70,'(2(1x,7es14.6))') (wdim(ix),ix=1,iw)
      enddo
      if(iout.ne.0) close(70)

      if(iout.ne.0) then
        if(LSN) then
          call b2wdat_open(70,'res_at_outer_plate_bcpot=11')
        else
          call b2wdat_open(70,'res_at_inner_plate_bcpot=11')
        endif
      endif
      do iy = -1, ny
        iw = 0
        do ix = -1, nx
          if(rightix(ix,iy).eq.nx+1) then
            if(iout.ne.0) then
              s1 = abs(pbs(ix,iy,0))
              vte = sqrt(te(ix,iy)/me)
              CS = SQRT(PZ(ix,iy)/RZ(ix,iy))
              t0 = 0.0_R8
              do is = 0, ns-1
                t0 = t0+rza(ix,iy,is)*cs*s1*na(ix,iy,is)
              enddo
              fchi = qe*t0
              fche = sqrt(1.0_R8/(2.0_R8*pi))*qe*ne(ix,iy)*vte*s1*
     &         exp(-qe*po(ix,iy)/te(ix,iy))
              fch_real_bc = -((1.0_R8-gammae)*fche-fchi)
              res = fch_real_bc - fch(ix,iy,0)
            endif
            fch_sch = -(sch(ix,iy,0)+sch(ix,iy,1)*po(ix,iy)+
     &                  sch(ix,iy,2)*ne(ix,iy)+
     &                  sch(ix,iy,3)*po(ix,iy)*ne(ix,iy))
            res_bc = fch_sch - fch(ix,iy,0)
            if(abs(res_bc).ge.abs(fch(ix,iy,0))*epsr+epsa)
     >       conv=.false.

            if(iout.ne.0) then
              seec = gammae
              t0 = max(fchi,(1.0_R8-seec)*fche)*qe/te(ix,iy)
              fch_code = -((1.0_R8-seec)*fche-fchi)
              dif_real_bc = fch_real_bc - fch_sch
              if(nmb) then
                wdim(iw+1) = iy
              else
                wdim(iw+1) = crdy_c_cell(ix,iy)
              endif
              wdim(iw+2) = fch_real_bc
              wdim(iw+3) = fch(ix,iy,0)
              wdim(iw+4) = res
              wdim(iw+5) = fch_code
              wdim(iw+6) = res_bc ! dif
              wdim(iw+7) = dif_real_bc
              iw=iw+7
              call xertst(iw.le.maxd,
     &         'index exceeded maxd=14 in b2urbp.F')
            endif
          endif
        enddo
        if(iout.ne.0) write(70,'(2(1x,7es14.6))') (wdim(ix),ix=1,iw)
      enddo
      if(iout.ne.0) close(70)

      ncall = ncall+1
      return

      contains

      subroutine b2wdat_open(nwrite,first_name)
      integer, intent(in) :: nwrite
      character*(*), intent(in) :: first_name
      character*256 :: filename

      filename=trim(folder)//trim(first_name)//'.dat'
      if(iappend.eq.0) then
        open (nwrite,file=trim(filename))
      else
        open (nwrite,file=trim(filename),position='append')
      end if
      return
      end subroutine b2wdat_open

      end subroutine b2urbp
