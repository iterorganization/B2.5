      module b2mod_par_opt
      use b2mod_types
      use b2mod_user_namelist
     & , only : inf_opt, nsigma, sigma,
     & prior_type, prior_par, prior_range,
     & nsigmx, cftype, nvmx
      implicit none
      integer, save :: npar_opt = 0
      real (kind=R8), allocatable, save :: par_opt0(:), par_opt(:),
     &       par_opt_phys(:)
      logical, save :: flag_optim = .false.
      integer nnvar, nncon, nnjac
      double precision :: x0(nvmx), xl(nvmx), xu(nvmx), gl(nvmx),
     &  gu(nvmx), jcol(nvmx*nvmx), jrow(nvmx*nvmx), jj(nvmx*nvmx),
     &  par_rescale(nvmx)
! sc some variables that allows to change optimization parameters, read from b2.optimization.parameters
      character*256, save :: filename, limited_memory_update_type,
     &  hessian_approximation 
      double precision, save :: cpu_opt=0.0, tol_opt=1.0e-7
      integer, save :: maxiter=100, partype(nvmx), nsigma_opt = 0 
      integer, save :: paris(nvmx),parib(nvmx)
      logical, save :: sigma_opt(nsigmx)
      logical, save :: spatial_dep(nvmx)
      integer, save :: spatial_points(nvmx) = 0
      integer :: idd, ipp
      logical :: file_ok
      real (kind=R8) :: rangeold(nvmx,2), parold(nvmx,2)
      integer :: typeold(nvmx)
      private optimization
      namelist /optimization/ nnvar, nncon, nnjac, x0, xl, xu, gl, gu, 
     & jcol, jrow, jj, maxiter, cpu_opt, tol_opt, hessian_approximation,
     & limited_memory_update_type, partype, sigma_opt, spatial_dep,
     & spatial_points, par_rescale, paris, parib

      contains

      subroutine read_b2mod_par_opt(ncon, nele_jac, ns, nbc)
      use b2mod_types
      implicit none
      integer ncon, nele_jac, ns, nbc, ii, isigma

      filename='b2.optimization.parameters'
      write (*,*) 'OPTIM: max number of readable parameters '//
     & 'and constraints:',nvmx
      write (*,*) 'OPTIM: infinity bound set to:',inf_opt
      x0 = inf_opt*10.0
      xl = -inf_opt*10.0
      xu = inf_opt*10.0
      gl = -inf_opt*10.0
      gu = inf_opt*10.0
      par_rescale = 1.0
      NELE_JAC = 0
      nnjac = 0
      jcol = 0
      jrow = 0
      jj = 0.0
      partype = -1
      sigma_opt = .true.
      spatial_dep = .false.
      spatial_points = 0
      hessian_approximation = 'limited-memory'
      limited_memory_update_type = 'b2fgs'
! FIXME set if (read_file)
      call find_file(filename,file_ok)
      if(file_ok) then
        open(99,file=filename)
        write(*,*) 'Reading namelist optimization from ',trim(filename)
#ifdef USE_MPI
        read(99,optimization)
#else
        read(99,optimization)
#endif
        close(99)
        write(*,*) 'Read namelist optimization from ',trim(filename)
      else
        call xerrab ('No '//trim(filename))
      endif
      write (*,*) 'OPTIM: read parameters'
      write(*,optimization)
      call xertst(nnvar.gt.0.and.nnvar.le.nvmx,
     & 'b2mod_par_opt: nnvar<=0 or nnvar>nvmx')
      call xertst(nncon.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnjac.le.nvmx,'b2mod_par_opt: increase size of nvmx')
      call xertst(nnvar.ge.0,'b2mod_par_opt: nnvar<0')
      call xertst(nncon.ge.0,'b2mod_par_opt: nncon<0')
      call xertst(nnjac.ge.0,'b2mod_par_opt: nnjac<0')
      call xertst(all(par_rescale.gt.0.D0),
     & 'b2mod_par_opt: par_rescale<=0')
!     partype refers to the type of design variable:
!     0  - sigma, standard deviation on parameters for Bayesian inference. They must be the last parameters in the vector!
!     1  - Dna, density driven particle diffusion coefficient
!     2  - Dpa, pressure driven particle diffusion coefficient
!     3  - X_i, ion heat diffusion coefficient
!     4  - X_e, electron heat diffusion coefficient
!     5  - vlax, x component of the anomalous velocity
!     6  - vlay, y component of the anomalous velocity FIXME not yet done
!     7  - vsa, anomalous viscosity
!     8  - sig, field-driven current radial conductivity
!     9  - alf, temperature-driven current radial conductivity
!     10 - enepar(parib(ipar), 1), with parib(ipar) identifies boundary
!     11 - enipar(parib(ipar), 1), with parib(ipar) identifies boundary
!     12 - conpar(paris(ipar), parib(ipar), 1), with paris(ipar) identifies species, parib(ipar) identifies boundary
!     13 - mompar(paris(ipar), parib(ipar), 1), with paris(ipar) identifies species, parib(ipar) identifies boundary
!     14 - potpar(parib(ipar), 1), with parib(ipar) identifies boundary
!     15 - enkpar(parib(ipar), 1), k boundary condition, parib(ipar) identifies boundary
!     16 - b2recyc(paris(ipar),parib(ipar)), with paris(ipar) identifies species, parib(ipar) identifies stratum
!     17 - b2tqna_ballooning
!     18 - b2tqna_ballooning_rescale
!     19 - keps_cd, k-model coefficient in Dna_ExB ~ Cd*sqrt(k) 
!     20 - keps_heat, multiplier Ce for X_e = Ce*Dna_ExB 
!     21 - keps_heat_i, multiplier Ci for X_i = Ci*Dna_ExB 
!     22 - keps_sig
!     23 - keps_alf
!     24 - keps_visc
!     25 - keps_dkt
!     26 - keps_dzt
!     27 - b2sikt_fac_diss
!     28 - b2sikt_fac_diss_core
!     29 - b2sikt_fac_sheath, multiplier Cl for Sk_diss = Cl*k*cs/Lconn (valid in SOL & PFR)
!     30 - b2sikt_fac_sheath_core, multiplier Cl for Sk_diss = Cl*k*cs/Lconn (valid in core)
!     31 - keps_shear
!     32 - b2tfhi_fsigkt
!     33 - b2sikt_fac_vis_RS
!     34 - b2tfhi_fflokt
!     35 - b2tfhi_fconkt
!     36 - b2tfhi_fflozt
!     37 - b2tfhi_fconzt
!     38 - b2tfhi_fkt_hie
!     39 - b2tfhe_vis_kt
      call xertst(any(partype(1:nnvar) .ge. 0) .and.
     &  any(partype(1:nnvar) .le. 39),
     &  'b2mod_par_opt: partype<0 or partype>39')
      call xertst(any(spatial_points(1:nnvar) .ge. 0),
     &  'b2mod_par_opt: spatial_points<0')
      call xertst(sum(spatial_points(1:nnvar)) .le. nvmx,
     &  'b2mod_par_opt: sum(spatial_points)<=nvmx, increase nvmx')
      do ii=1,nnvar
       if (partype(ii) .eq. 1 .or. partype(ii) .eq. 2 .or.
     &     partype(ii) .eq. 3 .or. partype(ii) .eq. 5 .or.
     &     partype(ii) .eq. 7 .or. partype(ii) .eq. 12 .or.
     &     partype(ii) .eq. 13 .or. partype(ii) .eq. 16) 
     >  call xertst((paris(ii) .ge. 0) .and. (paris(ii) .le. ns-1),
     &  'b2mod_par_opt: paris<0 or >ns-1')
       if (partype(ii) .ge. 10 .and. partype(ii) .le. 16)
     >  call xertst((parib(ii) .ge. 1) .and. (parib(ii) .le. nBc),
     &   'b2mod_par_opt: parib<0 or >nBc')
      enddo
      if (nsigma .gt. 0 .and. any(sigma_opt(1:nsigma))) then
!      calculate first how many sigmas are being optimized
       nsigma_opt = 0
       do ii = 1, nsigma
        if (sigma_opt(ii)) nsigma_opt = nsigma_opt + 1
       end do
!      check if sigmas are at the end of the parameter vector, if not, issue error
       isigma = nnvar - nsigma_opt + 1 ! position of first sigma in partype
       do ii = 1, nsigma_opt
        call xertst(partype(ii+isigma-1) .eq. 0,'b2mod_par_opt: '//
     &   'partype 0 variables MUST be at the end of the vector')
       end do
      endif
!     Suppose MAP estimation is turned on only when cftype(1)=6
      if (cftype(1) .eq. 6) then
!      prior_type refers to the prior type for the design variables:
!      0 - Uniform distribution.
!      1 - Uninformative, proper, Gaussian prior.
!      2 - Jeffrey's 1/sigma (only for sigma!!)
       call xertst(any(prior_type(1:nnvar) .ge. 0).and.
     &  any(prior_type(1:nnvar) .le. 2),
     &  'b2mod_par_opt: prior_type<0 or prior_type>2')
!      prior_par refers to the prior parameters:
!      For prior_type = 1, prior_par(ii,1) = mu, prior_par(ii,2) = sigma
!      For prior_type = 2, not used?
       do ii = 1, nnvar
        if (prior_type(ii) .eq. 1) 
     &   call xertst(prior_par(ii,2) .gt. 0.0_R8,
     &    'b2mod_par_opt: sigma prior <=0')
        call xertst(prior_range(ii,1).lt.prior_range(ii,2),
     &   'b2mod_par_opt: prior_range(.,1)>prior_range(.,2)')
       end do
      end if
      if (nnjac .gt. 0) then
       call xertst(any(jcol(1:nnjac) .gt. 0),'b2mod_par_opt: jcol<=0')
       call xertst(any(jrow(1:nnjac) .gt. 0),'b2mod_par_opt: jrow<=0')
       call xertst(maxval(jcol(1:nnjac)) .le. nncon,
     & 'b2mod_par_opt: jcol>ncon')
       call xertst(maxval(jrow(1:nnjac)) .le. nnvar,
     & 'b2mod_par_opt: jrow>nvar')
      endif

!     adjust some arrays when there are radially dependent coefficients
!     x0, xl, xu are supposed to have the correct dimension already, i.e. one for each spatial point
!     priors parameters are not
      npar_opt = 0
      idd = 1
      typeold = prior_type
      parold = prior_par
      rangeold = prior_range
      do ipp = 1, nnvar + (nsigma - nsigma_opt)
       prior_type(idd) = typeold(ipp)
       prior_par(idd,1) = parold(ipp,1)
       prior_par(idd,2) = parold(ipp,2)
       prior_range(idd,1) = rangeold(ipp,1)
       prior_range(idd,2) = rangeold(ipp,2)
       if (spatial_dep(ipp)) then
        prior_type(idd:idd+spatial_points(ipp)-1) = typeold(ipp)
        prior_par(idd:idd+spatial_points(ipp)-1,1) = parold(ipp,1)
        prior_par(idd:idd+spatial_points(ipp)-1,2) = parold(ipp,2)
        prior_range(idd:idd+spatial_points(ipp)-1,1) = rangeold(ipp,1)
        prior_range(idd:idd+spatial_points(ipp)-1,2) = rangeold(ipp,2)
        npar_opt = npar_opt + spatial_points(ipp) - 1 
        idd = idd + spatial_points(ipp) - 1 !next index
       end if
       npar_opt = npar_opt + 1
       idd = idd + 1 !next index
      enddo
      npar_opt = npar_opt - (nsigma - nsigma_opt)
      isigma = npar_opt - nsigma_opt + 1 ! position of first sigma in x0
      do ii = 1, nsigma
!      only assign initial value from x0 if that sigma is being optimized
       if (sigma_opt(ii)) then
        sigma(ii) = dble(x0(isigma))
        isigma = isigma + 1 !next sigma in inital vector x0
       endif
      end do
      ncon = nncon
      nele_jac = nnjac
      write (*,*) 'NPAR_OPT = ',npar_opt
      write (*,*) 'NCON = ',ncon
      write (*,*) 'NELE_JAC = ',NELE_JAC

      return
      end subroutine read_b2mod_par_opt

      end module b2mod_par_opt
