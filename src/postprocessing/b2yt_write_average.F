
      subroutine b2yt_write_average (iun, nx1, ny1, ns1,
     &    naver, ntotdt, lblmn, zamin1, zamax1, zn1, am1,
     &    na1_mean, te1_mean, ti1_mean, ua1_mean, po1_mean,
     &    e_na1, e_te1, e_ti1, e_ua1, e_po1,
     &    sna1_mean, smo1_mean, smr1_mean, smd1_mean,
     &    she1_mean, shi1_mean,
     &    e_sna1, e_smo1, e_smr1, e_smd1, e_she1, e_shi1)

      use b2mod_types
      use b2mod_version

      implicit none

      integer, intent(in) :: iun, nx1, ny1, ns1, naver, ntotdt
      character*120 :: lblmn

      real (kind=R8), intent(in) ::
     .  zamin1(0:ns1-1), zamax1(0:ns1-1), zn1(0:ns1-1), am1(0:ns1-1),
     .  na1_mean(-1:nx1,-1:ny1,0:ns1-1), te1_mean(-1:nx1,-1:ny1),
     .  ti1_mean(-1:nx1,-1:ny1), ua1_mean(-1:nx1,-1:ny1,0:ns1-1),
     .  po1_mean(-1:nx1,-1:ny1),
     .  e_na1(-1:nx1,-1:ny1,0:ns1-1), e_te1(-1:nx1,-1:ny1),
     .  e_ti1(-1:nx1,-1:ny1), e_ua1(-1:nx1,-1:ny1,0:ns1-1),
     .  e_po1(-1:nx1,-1:ny1),
     .  sna1_mean(-1:nx1,-1:ny1,0:ns1-1),
     .  smo1_mean(-1:nx1,-1:ny1,0:ns1-1),
     .  smr1_mean(-1:nx1,-1:ny1,0:ns1-1),
     .  smd1_mean(-1:nx1,-1:ny1,0:ns1-1),
     .  she1_mean(-1:nx1,-1:ny1),
     .  shi1_mean(-1:nx1,-1:ny1),
     .  e_sna1(-1:nx1,-1:ny1,0:ns1-1),
     .  e_smo1(-1:nx1,-1:ny1,0:ns1-1),
     .  e_smr1(-1:nx1,-1:ny1,0:ns1-1),
     .  e_smd1(-1:nx1,-1:ny1,0:ns1-1),
     .  e_she1(-1:nx1,-1:ny1),
     .  e_shi1(-1:nx1,-1:ny1)

      integer :: n2, idum(0:2)

      n2 = (nx1+2)*(ny1+2)

*  write 2d averaged plasma profiles

      open (unit=iun,file='b2favere.2',access='SEQUENTIAL',
     .          form='FORMATTED')

*   ..write header
      idum(0) = nx1
      idum(1) = ny1
      idum(2) = ns1
      call cfverw (iun, newversion)
      call cfwuin (iun, 3, idum, 'nx,ny,ns')
      call cfwuch (iun, 120, lblmn, 'label')
      call b2wuzd (iun, newversion, ns1, zamin1, zamax1, zn1, am1)

*   ..write average state and sources
      idum(0) = naver
      call cfwuin (iun, 1, idum, 'naver')
      idum(0) = ntotdt
      call cfwuin (iun, 1, idum, 'ntotdt')

      call cfwure (iun, n2*ns1, na1_mean, 'na_mean')
      call cfwure (iun, n2,     te1_mean, 'te_mean')
      call cfwure (iun, n2,     ti1_mean, 'ti_mean')
      call cfwure (iun, n2,     po1_mean, 'po_mean')
      call cfwure (iun, n2*ns1, ua1_mean, 'ua_mean')

      call cfwure (iun, n2*ns1, sna1_mean, 'sna_mean')
      call cfwure (iun, n2,     she1_mean, 'she_mean')
      call cfwure (iun, n2,     shi1_mean, 'shi_mean')
      call cfwure (iun, n2*ns1, smo1_mean, 'smo_mean')
      call cfwure (iun, n2*ns1, smr1_mean, 'smr_mean')
      call cfwure (iun, n2*ns1, smd1_mean, 'smd_mean')

*  write 2d squared averaged plasma profiles
      call cfwure (iun, n2*ns1, e_na1, 'e_na')
      call cfwure (iun, n2,     e_te1, 'e_te')
      call cfwure (iun, n2,     e_ti1, 'e_ti')
      call cfwure (iun, n2,     e_po1, 'e_po')
      call cfwure (iun, n2*ns1, e_ua1, 'e_ua')

      call cfwure (iun, n2*ns1, e_sna1, 'e_sna')
      call cfwure (iun, n2,     e_she1, 'e_she')
      call cfwure (iun, n2,     e_shi1, 'e_shi')
      call cfwure (iun, n2*ns1, e_smo1, 'e_smo')
      call cfwure (iun, n2*ns1, e_smr1, 'e_smr')
      call cfwure (iun, n2*ns1, e_smd1, 'e_smd')

      close (iun)
      write (*,*) ' file b2faveri.2 written '

      return
      end subroutine b2yt_write_average
