      subroutine map_and_interpolate_cf(n1, n2, nomp, b2psi, b2rr,
     &  psi_omp, b2data, cfrr, b2dataOnCf, mapToOmp, ompind, shift)
      use b2mod_types
      use b2mod_subsys
      implicit none
      integer, intent (in) :: n1, n2, nomp, ompind(1:2)
      real(kind=R8), intent (in) :: b2psi(1:n1), psi_omp(1:nomp),
     &  b2data(1:n1), b2rr(1:n1), cfrr(1:n2), shift
      real(kind=R8), intent (inout) :: b2dataOnCf(1:n2)
      logical, intent(in) :: mapToOmp
      ! local variables
      real(kind=R8) :: b2datapsi(1:nomp)
      integer :: i1, i2, nn
      ! This routine interpolates data for the cost function calculation
      ! from the B2.5 CVs (or FCs) onto spatial points where the
      ! cost function is defined.
      ! If the cost function is mapped to the OMP, then this mapping is done
      ! as: first interpolate B2.5 data from the cost function domain onto
      ! the OMP CVs using psi; second from OMP CVs onto cost function spatial
      ! points.  

      call subini('map_and_interpolate_cf')

      if (mapToOmp) then
        ! Mapping might need excluding some CVs on the OMP using ompind
        i1 = ompind(1)
        i2 = ompind(2)
        nn = i2 - i1 + 1
        ! Now interpolate b2data from CF domain (b2psi) to OMP domain (psi_omp)
        b2datapsi = 0.0_R8
        call interp1d(n1, nn, b2psi, psi_omp(i1:i2), b2data,
     &    b2datapsi(1:nn), 0.0_R8)
        ! now interpolate data from b2rr on OMP onto cfdata
        ! (note that b2rr has already been adjusted using ompind)
        call interp1d(nn, n2, b2rr(1:nn), cfrr, b2datapsi(1:nn),
     &    b2dataOnCf, shift)
      else
        ! now interpolate data from b2rr onto cfdata
        call interp1d(n1, n2, b2rr, cfrr, b2data, b2dataOnCf, shift)
      endif

      call subend ()
      return

      end subroutine map_and_interpolate_cf

