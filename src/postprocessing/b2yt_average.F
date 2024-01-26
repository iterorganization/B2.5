*-----------------------------------------------------------------------
*.specification

      subroutine b2yt_average (filename,
     &  nx, ny, ns, nx1, ny1, ns1, nncutmax,
     &  non_commensurate, flip,
     &  nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &  leftcut, leftcut1, rightcut, rightcut1,
     &  bottomcut, bottomcut1, topcut, topcut1,
     &  bundle_frac, naver, ntotdt,
     &  na_mean, te_mean, ti_mean, ua_mean, po_mean,
     &  e_na, e_te, e_ti, e_ua, e_po,
     &  sna_mean, sne_mean, smo_mean, smr_mean, smd_mean,
     &  she_mean, shi_mean, sch_mean,
     &  e_sna, e_sne, e_smo, e_smr, e_smd, e_she, e_shi, e_sch,
     &  na1_mean, te1_mean, ti1_mean, ua1_mean, po1_mean,
     &  e_na1, e_te1, e_ti1, e_ua1, e_po1,
     &  sna1_mean, sne1_mean, smo1_mean, smr1_mean, smd1_mean,
     &  she1_mean, shi1_mean, sch1_mean,
     &  e_sna1, e_sne1, e_smo1, e_smr1, e_smd1,
     &  e_she1, e_shi1, e_sch1)
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_version
      use b2mod_constants
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      character*256 filename
      integer nx, ny, ns, nx1, ny1, ns1, nncutmax, non_commensurate
      integer nnreg(0:2), nncut, periodic_bc, jsep, jsep1
      integer region(-1:nx,-1:ny,0:2), region1(-1:nx1,-1:ny1,0:2)
      integer leftcut(nncutmax), leftcut1(nncutmax)
      integer rightcut(nncutmax), rightcut1(nncutmax)
      integer bottomcut(nncutmax), bottomcut1(nncutmax)
      integer topcut(nncutmax), topcut1(nncutmax)
      real (kind=R8) ::
     *  bundle_frac(0:ns1-1,0:ns-1)
      real (kind=R8) ::
     . na_mean(-1:nx,-1:ny,0:ns-1),
     . te_mean(-1:nx,-1:ny), ti_mean(-1:nx,-1:ny),
     . ua_mean(-1:nx,-1:ny,0:ns-1), po_mean(-1:nx,-1:ny),
     . e_na(-1:nx,-1:ny,0:ns-1), e_te(-1:nx,-1:ny), e_ti(-1:nx,-1:ny),
     . e_ua(-1:nx,-1:ny,0:ns-1), e_po(-1:nx,-1:ny),
     . sna_mean(-1:nx,-1:ny,0:ns-1), sne_mean(-1:nx,-1:ny),
     . smo_mean(-1:nx,-1:ny,0:ns-1), she_mean(-1:nx,-1:ny),
     . smr_mean(-1:nx,-1:ny,0:ns-1), smd_mean(-1:nx,-1:ny,0:ns-1),
     . shi_mean(-1:nx,-1:ny), sch_mean(-1:nx,-1:ny),
     . e_sna(-1:nx,-1:ny,0:ns-1), e_sne(-1:nx,-1:ny),
     . e_smo(-1:nx,-1:ny,0:ns-1),
     . e_smr(-1:nx,-1:ny,0:ns-1), e_smd(-1:nx,-1:ny,0:ns-1),
     . e_she(-1:nx,-1:ny), e_shi(-1:nx,-1:ny),
     . e_sch(-1:nx,-1:ny)
      logical flip
*   ..output arguments (unspecified on entry)
      integer naver, ntotdt
      real (kind=R8) ::
     . na1_mean(-1:nx1,-1:ny1,0:ns1-1),
     . te1_mean(-1:nx1,-1:ny1), ti1_mean(-1:nx1,-1:ny1),
     . ua1_mean(-1:nx1,-1:ny1,0:ns1-1), po1_mean(-1:nx1,-1:ny1),
     . e_na1(-1:nx1,-1:ny1,0:ns1-1),
     . e_te1(-1:nx1,-1:ny1), e_ti1(-1:nx1,-1:ny1),
     . e_ua1(-1:nx1,-1:ny1,0:ns1-1), e_po1(-1:nx1,-1:ny1),
     . sna1_mean(-1:nx1,-1:ny1,0:ns1-1),
     . sne1_mean(-1:nx1,-1:ny1),
     . smo1_mean(-1:nx1,-1:ny1,0:ns1-1),
     . smr1_mean(-1:nx1,-1:ny1,0:ns1-1),
     . smd1_mean(-1:nx1,-1:ny1,0:ns1-1),
     . she1_mean(-1:nx1,-1:ny1),
     . shi1_mean(-1:nx1,-1:ny1),
     . sch1_mean(-1:nx1,-1:ny1),
     . e_sna1(-1:nx1,-1:ny1,0:ns1-1),
     . e_sne1(-1:nx1,-1:ny1),
     . e_smo1(-1:nx1,-1:ny1,0:ns1-1),
     . e_smr1(-1:nx1,-1:ny1,0:ns1-1),
     . e_smd1(-1:nx1,-1:ny1,0:ns1-1),
     . e_she1(-1:nx1,-1:ny1),
     . e_shi1(-1:nx1,-1:ny1),
     . e_sch1(-1:nx1,-1:ny1)

*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2YT_AVERAGE transfers the averages plasma state to a mesh of
*     different size.
*
*
*  3. description (see also routine b2cdca)
*
*     This routine transfers the averaged plasma state to a mesh of
*     different size.  Because of possible cuts in the geometry and
*     associated discontinuities, the lowest order (nearest-neighbour)
*     form of interpolation is employed.
*     In the case of a non-conformal mapping, then we do a 2-D
*     interpolation of values in the computational domain, per region.
*
*
*  5. parameters (see also routine b2cdcv)
*
*     nx, ny, ns, nx1, ny1 - integers, input.
*     na_mean, te_mean, ti_mean, ua_mean, po_mean,
*     sna_mean, sne_mean, smo_mean, smr_mean, smd_mean, she_mean,
*     shi_mean, sch_mean,
*     e_na, e_te, e_ti, e_ua, e_po,
*     e_sna, e_sne, e_smo, e_smr, e_smd,
*     e_she, e_shi, e_sch - real*8 arrays, input.
*     See the calling routine for description.
*
*     ne1, na1, ua1, uadia1, te1, ti1, po1 - real*8 arrays, output.
*     na1_mean, te1_mean, ti1_mean, ua1_mean, po1_mean,
*     sna1_mean, sne1_mean, smo1_mean, smr1_mean, smd1_mean,
*     she1_mean, shi1_mean, sch1_mean,
*     e_na1, e_te1, e_ti1, e_ua1, e_po1,
*     e_sna1, e_sne1, e_smo1, e_smr1, e_smd1,
*     e_she1, e_shi1, e_sch1 - real*8 arrays, output.
*     Transformed averaged plasma state.
*
*     (parameters in common:)
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer is, is1, n1, n2, idum(0:2)
      real (kind=R8) ::
     *   wrk0(-1:nx,-1:ny1), wrk1(-1:nx1,-1:ny1), na_new, na_min, sum,
     *   zpdum(0:ns-1,4), zero(-1:nx,-1:ny), rza(-1:nx,-1:ny,0:ns-1),
     *   swrk2(-1:nx,-1:ny,0:1), swrk2s(-1:nx,-1:ny,0:1,0:ns-1),
     *   swrk4(-1:nx,-1:ny,0:3), swrk4s(-1:nx,-1:ny,0:3,0:ns-1),
     *   ne_mean(-1:nx,-1:ny), ni_mean(-1:nx,-1:ny,0:1)
      character*120 :: label
      save na_new, na_min
      data na_new/1.0e14_R8/, na_min/1.0e4_R8/
*   ..procedures
      external xertst, ipgetr, sfill, b2saxpy
      external b2xtvc, b2xpne, b2xpni

*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2yt_average')
*   ..test nx, ny, ns, nx1, ny1, ns1
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      call xertst (1.le.ns, 'faulty argument ns')
      call xertst (0.le.nx1.and.0.le.ny1, 'faulty argument nx1, ny1')
      call xertst (1.le.ns1, 'faulty argument ns1')
*   ..read input parameters
      call ipgetr ('b2mndr_na_min', na_min)
      call ipgetr ('b2mndr_na_new', na_new)
      call xertst (na_min.gt.0.0_R8,
     & 'faulty internal parameter na_min')
      call xertst (na_new.gt.0.0_R8,
     & 'faulty internal parameter na_new')
      call xertst (na_new.gt.na_min,
     & 'na_new must be larger than na_min!')
* ..initialize original arrays for averaged state
      n2 = (nx+2)*(ny+2)
      call sfill (n2*ns, 0.0_R8, na_mean, 1)
      call sfill (n2*ns, 0.0_R8, ua_mean, 1)
      call sfill (n2, 0.0_R8, te_mean, 1)
      call sfill (n2, 0.0_R8, ti_mean, 1)
      call sfill (n2, 0.0_R8, po_mean, 1)
      call sfill (n2*ns, 0.0_R8, e_na, 1)
      call sfill (n2*ns, 0.0_R8, e_ua, 1)
      call sfill (n2, 0.0_R8, e_te, 1)
      call sfill (n2, 0.0_R8, e_ti, 1)
      call sfill (n2, 0.0_R8, e_po, 1)
      call sfill (n2*ns, 0.0_R8, sna_mean, 1)
      call sfill (n2*ns, 0.0_R8, sne_mean, 1)
      call sfill (n2*ns, 0.0_R8, smo_mean, 1)
      call sfill (n2*ns, 0.0_R8, smr_mean, 1)
      call sfill (n2*ns, 0.0_R8, smd_mean, 1)
      call sfill (n2, 0.0_R8, she_mean, 1)
      call sfill (n2, 0.0_R8, shi_mean, 1)
      call sfill (n2, 0.0_R8, sch_mean, 1)
      call sfill (n2*ns*2, 0.0_R8, e_sna, 1)
      call sfill (n2*ns*2, 0.0_R8, e_sne, 1)
      call sfill (n2*ns*4, 0.0_R8, e_smo, 1)
      call sfill (n2*ns*4, 0.0_R8, e_smr, 1)
      call sfill (n2*ns*4, 0.0_R8, e_smd, 1)
      call sfill (n2*4, 0.0_R8, e_she, 1)
      call sfill (n2*4, 0.0_R8, e_shi, 1)
      call sfill (n2*4, 0.0_R8, e_sch, 1)
* ..read 2d averaged plasma profiles from b2favere file
      open(unit=99,file=trim(filename),access='SEQUENTIAL',
     .        form='FORMATTED',status='OLD')
      call cfverr (99, b2faveri_version)
      call cfruin (99, 3, idum, 'nx,ny,ns')
      call cfruch (99, 120, label, 'label')
      call b2ruzd (99, b2faveri_version, ns,
     &  zpdum(0,1), zpdum(0,2), zpdum(0,3), zpdum(0,4), .false.)

      call cfruin (99, 1, idum, 'naver')
      naver = idum(0)
      call cfruin (99, 1, idum, 'ntotdt')
      ntotdt = idum(0)
      call cfrure (99, n2*ns, na_mean, 'na_mean')
      call cfrure (99, n2,    te_mean, 'te_mean')
      call cfrure (99, n2,    ti_mean, 'ti_mean')
      call cfrure (99, n2,    po_mean, 'po_mean')
      call cfrure (99, n2*ns, ua_mean, 'ua_mean')

      if (b2faveri_version.ge.'03.000.008') then
        call cfrure (99, n2*ns, sna_mean, 'sna_mean')
        call cfrure (99, n2,    sne_mean, 'sne_mean')
        call cfrure (99, n2,    she_mean, 'she_mean')
        call cfrure (99, n2,    shi_mean, 'shi_mean')
        call cfrure (99, n2,    sch_mean, 'sch_mean')
        call cfrure (99, n2*ns, smo_mean, 'smo_mean')
        if (b2faveri_version.ge.'03.000.009'.and.
     &      b2faveri_version.ne.'03.001.000') then
          call cfrure (99, n2*ns, smr_mean, 'smr_mean')
          call cfrure (99, n2*ns, smd_mean, 'smd_mean')
        end if
      else
        zero = 0.0_R8
        do is = 0, ns-1
          rza(:,:,is) = (zamin(is)+zamax(is))/2.0_R8
        end do
        call b2xpne ( nx, ny, ns, rza, na_mean, zero, ne_mean)
        call b2xpni ( nx, ny, ns, na_mean, ni_mean)
        call cfrure (99, n2*2*ns, swrk2s, 'sna_mean')
        sna_mean(:,:,:) = swrk2s(:,:,0,:) +
     .                    swrk2s(:,:,1,:) * na_mean(:,:,:)
        call cfrure (99, n2*2,    swrk2 , 'sne_mean')
        sne_mean(:,:) = swrk2(:,:,0) +
     .                  swrk2(:,:,1) * ne_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'she_mean')
        she_mean(:,:) = swrk4(:,:,0) +
     .                  swrk4(:,:,1) * ne_mean(:,:) +
     .                  swrk4(:,:,2) * te_mean(:,:) +
     .                  swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'shi_mean')
        shi_mean(:,:) = swrk4(:,:,0) +
     .                  swrk4(:,:,1) * ni_mean(:,:,0) +
     .                  swrk4(:,:,2) * ti_mean(:,:) +
     .                  swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'sch_mean')
        sch_mean(:,:) = swrk4(:,:,0) +
     .                  swrk4(:,:,1) * po_mean(:,:) +
     .                  swrk4(:,:,2) * ne_mean(:,:) +
     .                  swrk4(:,:,3) * ne_mean(:,:) * po_mean(:,:)
        call cfrure (99, n2*4*ns, swrk4s, 'smo_mean')
        do is = 0, ns-1
          smo_mean(:,:,is) = swrk4s(:,:,0,is) +
     .     swrk4s(:,:,1,is) * ua_mean(:,:,is) +
     .     swrk4s(:,:,2,is) * na_mean(:,:,is) * mp * am(is) +
     .     swrk4s(:,:,1,is) * ua_mean(:,:,is) *
     .                        na_mean(:,:,is) * mp * am(is)
        end do
      end if

      call cfrure (99, n2*ns, e_na, 'e_na')
      call cfrure (99, n2,    e_te, 'e_te')
      call cfrure (99, n2,    e_ti, 'e_ti')
      call cfrure (99, n2,    e_po, 'e_po')
      call cfrure (99, n2*ns, e_ua, 'e_ua')

      if (b2faveri_version.ge.'03.000.008') then
        call cfrure (99, n2*ns, e_sna, 'e_sna')
        call cfrure (99, n2,    e_sne, 'e_sne')
        call cfrure (99, n2,    e_she, 'e_she')
        call cfrure (99, n2,    e_shi, 'e_shi')
        call cfrure (99, n2,    e_sch, 'e_sch')
        call cfrure (99, n2*ns, e_smo, 'e_smo')
        if (b2faveri_version.ge.'03.000.009'.and.
     &      b2faveri_version.ne.'03.001.000') then
          call cfrure (99, n2*ns, e_smr, 'e_smr')
          call cfrure (99, n2*ns, e_smd, 'e_smd')
        end if
      else
        call cfrure (99, n2*2*ns, swrk2s, 'e_sna')
        e_sna(:,:,:) = swrk2s(:,:,0,:) +
     .                 swrk2s(:,:,1,:) * na_mean(:,:,:)
        call cfrure (99, n2*2,    swrk2 , 'e_sne')
        e_sne(:,:) = swrk2(:,:,0) +
     .               swrk2(:,:,1) * ne_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'e_she')
        e_she(:,:) = swrk4(:,:,0) +
     .               swrk4(:,:,1) * ne_mean(:,:) +
     .               swrk4(:,:,2) * te_mean(:,:) +
     .               swrk4(:,:,3) * ne_mean(:,:) * te_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'e_shi')
        e_shi(:,:) = swrk4(:,:,0) +
     .               swrk4(:,:,1) * ni_mean(:,:,0) +
     .               swrk4(:,:,2) * ti_mean(:,:) +
     .               swrk4(:,:,3) * ni_mean(:,:,0) * ti_mean(:,:)
        call cfrure (99, n2*4,    swrk4 , 'e_sch')
        e_sch(:,:) = swrk4(:,:,0) +
     .               swrk4(:,:,1) * po_mean(:,:) +
     .               swrk4(:,:,2) * ne_mean(:,:) +
     .               swrk4(:,:,3) * ne_mean(:,:) * po_mean(:,:)
        call cfrure (99, n2*4*ns, swrk4s, 'e_smo')
        do is = 0, ns-1
          e_smo(:,:,is) = swrk4s(:,:,0,is) +
     .     swrk4s(:,:,1,is) * ua_mean(:,:,is) +
     .     swrk4s(:,:,2,is) * na_mean(:,:,is) * mp * am(is) +
     .     swrk4s(:,:,1,is) * ua_mean(:,:,is) *
     .                        na_mean(:,:,is) * mp * am(is)
        end do
      end if

      close(99)

* ..initialize transformed arrays for averaged state
      n1 = (nx1+2)*(ny1+2)
      call sfill (n1*ns1, 0.0_R8, na1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, ua1_mean, 1)
      call sfill (n1, 0.0_R8, te1_mean, 1)
      call sfill (n1, 0.0_R8, ti1_mean, 1)
      call sfill (n1, 0.0_R8, po1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, e_na1, 1)
      call sfill (n1*ns1, 0.0_R8, e_ua1, 1)
      call sfill (n1, 0.0_R8, e_te1, 1)
      call sfill (n1, 0.0_R8, e_ti1, 1)
      call sfill (n1, 0.0_R8, e_po1, 1)
      call sfill (n1*ns1, 0.0_R8, sna1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, sne1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, smo1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, smr1_mean, 1)
      call sfill (n1*ns1, 0.0_R8, smd1_mean, 1)
      call sfill (n1, 0.0_R8, she1_mean, 1)
      call sfill (n1, 0.0_R8, shi1_mean, 1)
      call sfill (n1, 0.0_R8, sch1_mean, 1)
      call sfill (n1*ns1*2, 0.0_R8, e_sna1, 1)
      call sfill (n1*ns1*2, 0.0_R8, e_sne1, 1)
      call sfill (n1*ns1*4, 0.0_R8, e_smo1, 1)
      call sfill (n1*ns1*4, 0.0_R8, e_smr1, 1)
      call sfill (n1*ns1*4, 0.0_R8, e_smd1, 1)
      call sfill (n1*4, 0.0_R8, e_she1, 1)
      call sfill (n1*4, 0.0_R8, e_shi1, 1)
      call sfill (n1*4, 0.0_R8, e_sch1, 1)

*   ..transform averaged field variables
      do is1 = 0, ns1-1
        sum = 0.0_R8
        do is = 0, ns-1
          sum = sum + bundle_frac(is1,is)
          if (bundle_frac(is1,is).ge.1.0_R8) then  ! simple case
            call b2xtvc (nx, ny, nx1, ny1,
     &        na_mean(-1,-1,is), na1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_na(-1,-1,is), e_na1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        ua_mean(-1,-1,is), ua1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_ua(-1,-1,is), e_ua1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        sna_mean(-1,-1,is), sna1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_sna(-1,-1,is), e_sna1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        smo_mean(-1,-1,is), smo1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_smo(-1,-1,is), e_smo1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
           call b2xtvc (nx, ny, nx1, ny1,
     &        smr_mean(-1,-1,is), smr1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_smr(-1,-1,is), e_smr1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
           call b2xtvc (nx, ny, nx1, ny1,
     &        smd_mean(-1,-1,is), smd1_mean(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
            call b2xtvc (nx, ny, nx1, ny1,
     &        e_smd(-1,-1,is), e_smd1(-1,-1,is1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1,
     &        wrk0)
          else if (bundle_frac(is1,is).gt.0.0_R8) then
            call b2xtvc (nx, ny, nx1, ny1, na_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        na1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, ua_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        ua1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_na(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_na1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_ua(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_ua1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, sna_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        sna1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_sna(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_sna1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, smo_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        smo1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_smo(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_smo1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, smr_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        smr1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_smr(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_smr1(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, smd_mean(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        smd1_mean(-1,-1,is1), 1)
            call b2xtvc (nx, ny, nx1, ny1, e_smd(-1,-1,is), wrk1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
            call b2saxpy (n2, bundle_frac(is1,is), wrk1, 1,
     &        e_smd1(-1,-1,is1), 1)
          endif
        enddo
        if (sum.eq.0.0_R8) then  ! this is a new species
          call sfill (n2, na_new, na1_mean(-1,-1,is1), 1)
        endif
      enddo
      call b2xtvc (nx, ny, nx1, ny1, te_mean, te1_mean,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, ti_mean, ti1_mean,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, po_mean, po1_mean,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_te, e_te1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_ti, e_ti1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_po, e_po1,
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1,
     &        sne_mean(-1,-1), sne1_mean(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_sne(-1,-1), e_sne1(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1,
     &        she_mean(-1,-1), she1_mean(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_she(-1,-1), e_she1(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1,
     &        shi_mean(-1,-1), shi1_mean(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_shi(-1,-1), e_shi1(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1,
     &        sch_mean(-1,-1), sch1_mean(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)
      call b2xtvc (nx, ny, nx1, ny1, e_sch(-1,-1), e_sch1(-1,-1),
     &        nncutmax, non_commensurate, flip,
     &        nnreg, nncut, periodic_bc, jsep, jsep1, region, region1,
     &        leftcut, leftcut1, rightcut, rightcut1,
     &        bottomcut, bottomcut1, topcut, topcut1, wrk0)

* ..return
      call subend ()
      return

*-----------------------------------------------------------------------
*.end b2yt_average

      end subroutine b2yt_average
