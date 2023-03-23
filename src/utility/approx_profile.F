      subroutine get_dim(fnm, ndat)                  !srv 24.04.09 {
      use b2mod_types
      implicit none
      integer :: ndat
      character*(*) fnm
      real (kind=R8) :: t0, t1
*   ..read data file in order to determine ndat
      open(70,file=trim(fnm))
      ndat = 0
      do
        read(70,*,end=10) t0, t1
        ndat = ndat + 1
      enddo
   10 close(70)

      end subroutine get_dim                               !srv 29.04.09 }

      subroutine read_data(fnm, ndat, cy, fy)        !srv 19.07.10 {
      use b2mod_types
      implicit none
      integer ndat
      character*(*) fnm
      real (kind=R8) :: cy(1:ndat), fy(1:ndat)
      integer iy
      open(70,file=trim(fnm))
      do iy = 1, ndat
        read(70,*) cy(iy), fy(iy)
      enddo
      close(70)
      end subroutine read_data                       !srv 19.07.10 }

      subroutine write_data(fnm, ndat, cy, fy)       !srv 19.07.10 {
      use b2mod_types
      implicit none
      integer ndat
      character*(*) fnm
      real (kind=R8) :: cy(1:ndat), fy(1:ndat)
      integer iy
      open(70,file=trim(fnm))
      do iy=1,ndat
        write(70,"(1x,2(es14.6))") cy(iy), fy(iy)
      enddo
      close(70)
      end subroutine write_data                        !srv 19.07.10 }

      subroutine approx_profile(fnm, ndat, r, f, b, c, d)  !srv 01.03.07 {
      use b2mod_types
      implicit none
      character*(*) fnm
      integer :: ndat
!      real (kind=R8) :: coord(-1:n), fun(-1:n)

      real (kind=R8) :: r(1:ndat), f(1:ndat), b(1:ndat), c(1:ndat),
     &                  d(1:ndat)
      integer :: i, iy

      open(70,file=trim(fnm))
      do iy = 1, ndat
        read(70,*) r(iy), f(iy)
      enddo
      close(70)
      do iy = 1, ndat-1
        call xertst (r(iy).lt. r(iy+1),
     1   'Sequence of coordinates in file '//fnm//' is not monotonic')
      enddo

!      call dcscon (ndat, r, f, ibreak, break, cscoef)
!      write(*,*) ' ibreak = ', ibreak
!      ideriv = 0
!      nintv  = ndat-1
!      call dcs1gd (ideriv, ndat, r, nintv, break, cscoef, pf)
!      do iy = 1,ndat
!        pf(iy-2) = dcsval(r(iy), nintv, break, cscoef)
!        write(*,*) r(iy), pf(iy-2)
!      enddo

!      CALL DSPLINE(NDAT, R, F, B, C, D)
*   ..the simplest linear interpolation
      do i = 1, ndat-1
        b(i) = (f(i+1)-f(i))/(r(i+1)-r(i))
      enddo
      b(ndat) = b(ndat-1)
      c = 0.0e0_R8
      d = 0.0e0_R8

      end subroutine approx_profile                                     !srv 01.03.07 }

      subroutine approx_profile_mod(ndat, r, f, b, c, d)                !srv 19.07.10 {
      use b2mod_types
      implicit none
      integer :: ndat
!       real (kind=R8) :: coord(-1:n), fun(-1:n)

      real (kind=R8) :: r(1:ndat), f(1:ndat), b(1:ndat), c(1:ndat),
     &                  d(1:ndat)
      integer :: i

!      integer :: iy
!
!      call dcscon (ndat, r, f, ibreak, break, cscoef)
!      write(*,*) ' ibreak = ', ibreak
!      ideriv = 0
!      nintv  = ndat-1
!      call dcs1gd (ideriv, ndat, r, nintv, break, cscoef, pf)
!      do iy = 1,ndat
!        pf(iy-2) = dcsval(r(iy), nintv, break, cscoef)
!        write(*,*) r(iy), pf(iy-2)
!      enddo

!      CALL DSPLINE(NDAT, R, F, B, C, D)
*   ..the simplest linear interpolation
      do i = 1, ndat-1
         b(i) = (f(i+1)-f(i))/(r(i+1)-r(i))
      enddo
      b(ndat) = b(ndat-1)
      c = 0.0e0_R8
      d = 0.0e0_R8

      end subroutine approx_profile_mod                                 !srv 19.07.10 }

