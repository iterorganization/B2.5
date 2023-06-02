C 2002.11.01
      module b2mod_ppout
      use b2mod_types
#ifdef NAGF95
      use f90_unix           ! IGNORE
#endif

      integer, save ::
     *  nxpp, nypp, nspp,
     *  iyeard

      integer, allocatable, save ::
     *  isyear(:)

#ifdef JET

      logical, save ::
     *  labort

      integer, save ::
     *  np, iopen, nxw, nc, nrow, jprgt, jplft,
     *  nsepx, nbcore, nbwall, ism, ism0, nz, nzs, limrec,
     *  iform

      real (kind=R8), save ::
     *     rpx, zpx, r0, z0, epsang

      real (kind=R8), allocatable, save ::
     *  xnorm(:,:,:),
     *  crxe(:,:,:), crye(:,:,:),
     *  cxpp(:,:),   cypp(:,:),
     *  psi(:,:),    zeff(:,:),
     *  rho(:,:),    hrho(:,:),
     *  theta(:,:),  htheta(:,:),
     *  bfi(:,:),    sh(:,:),
     *  dm(:,:),
     *  rsepx(:),    zsepx(:),
     *  zmass(:),    sqezrad(:,:,:),
     *  rz(:,:),     pz(:,:), pmach(:,:,:),
     *  rbcore(:),   zbcore(:),
     *  rbwall(:,:), zbwall(:,:),
     *  rtemp(:)

      integer, allocatable, save ::
     *  izs2is(:),
     *  izch(:),
     *  itagdv(:,:), itemp(:)

      character*80 adasid

      character*132, save :: labgeo

c   ..functions
      integer ixm1, ixp1
      character*120 gmtry_file

      real (kind=R8) ::
     *  vdot, vlen, thetasepx
#elif defined(IMAS)
      character*132, save :: labgeo
#endif

      contains
c
c=======================================================================
      subroutine alloc_b2mod_ppout(nxd, nyd, nsd)

      implicit none

      integer nxd , nyd , nsd

c     ------------------------------------------------------------------
c     ALLOC_PPOUT allocates PPOUT module arrays, etc
c     ------------------------------------------------------------------

c   ..initialise
      nxpp   = nxd
      nypp   = nyd
      nspp   = nsd

      iyeard = 89

c   ..allocate JET-independent arrays
      allocate(isyear(0:nsd-1))

#ifdef JET
      epsang = 0.01_R8
      limrec = 12000

      adasid = '/u/adas/adas'

      iform  = 0

c   ..allocate
      allocate(xnorm(-3:nxd+2,-1:nyd,0:1))
      allocate(crxe(-3:nxd+2,-1:nyd,0:3))
      allocate(crye(-3:nxd+2,-1:nyd,0:3))

      allocate(cxpp(-1:nxd,-1:nyd))
      allocate(cypp(-1:nxd,-1:nyd))
      allocate(psi(-1:nxd,-1:nyd))

      allocate(rho(-1:nxd,-1:nyd))
      allocate(hrho(-1:nxd,-1:nyd))
      allocate(theta(-1:nxd,-1:nyd))
      allocate(htheta(-1:nxd,-1:nyd))

      allocate(bfi(-1:nxd,-1:nyd))
      allocate(sh(-1:nxd,-1:nyd))

      allocate(dm(-1:nxd,-1:nyd))
      allocate(zeff(-1:nxd,-1:nyd))

      allocate(rsepx(1:nxd+3))
      allocate(zsepx(1:nxd+3))

      allocate(rbcore(1:nxd))
      allocate(zbcore(1:nxd))
      allocate(rbwall(1:2*(nxd+nyd+4),1:2))
      allocate(zbwall(1:2*(nxd+nyd+4),1:2))

      allocate(zmass(1:nsd))
      allocate(izch(1:nsd))
      allocate(izs2is(1:nsd))
      allocate(sqezrad(-1:nxd,-1:nyd,1:nsd))
      allocate(rz(-1:nxd,-1:nyd))   ! was wrongly declared with nsd as well
      allocate(pz(-1:nxd,-1:nyd))   ! was wrongly declared with nsd as well
      allocate(pmach(-1:nxd,-1:nyd,0:nsd-1))  ! was 1:nsd

      allocate(itagdv(-1:nxd,-1:nyd))

      allocate(rtemp(1:1))
      allocate(itemp(1:1))
#endif

      return
      end subroutine alloc_b2mod_ppout
c
c=======================================================================
      subroutine dealloc_b2mod_ppout
      implicit none

      deallocate(isyear)

#ifdef JET
c   ..deallocate
      deallocate(xnorm)
      deallocate(crxe)
      deallocate(crye)

      deallocate(cxpp)
      deallocate(cypp)
      deallocate(psi)

      deallocate(rho)
      deallocate(hrho)
      deallocate(theta)
      deallocate(htheta)

      deallocate(bfi)
      deallocate(sh)

      deallocate(dm)
      deallocate(zeff)

      deallocate(rsepx)
      deallocate(zsepx)

      deallocate(rbcore)
      deallocate(zbcore)
      deallocate(rbwall)
      deallocate(zbwall)

      deallocate(zmass)
      deallocate(izch)
      deallocate(izs2is)
      deallocate(sqezrad)
      deallocate(rz)
      deallocate(pz)
      deallocate(pmach)

      deallocate(itagdv)

      deallocate(rtemp)
      deallocate(itemp)
#endif

      return
      end subroutine dealloc_b2mod_ppout
c
c=======================================================================
#ifdef JET
      subroutine set_b2mod_ppout

      use b2mod_constants
      use b2mod_geo
      use b2mod_rates
      use b2mod_indirect
      use b2mod_plasma
      use b2mod_sources
      use b2mod_b2cmpa

      implicit none

c     ------------------------------------------------------------------
c     SET_B2MOD_PPOUT derives quantities for the PPOUT module.
c     ------------------------------------------------------------------

c   ..locals
      integer ix , iy , is , i , j
cdpc+xpb
      integer imin, imax
      real (kind=R8) rmin, rmax
cdpc+xpb-end

      logical lok

c   ..conditions for not producing a postprocessor file
      if( bottomcut(1).ne.-1 .or. nncut.ne.1 ) then
          labort = .true.
          return
      else
          labort = .false.
      endif

c   ..first open ring
      iopen = topcut(1) + 2

c   ..sol wall ring
      nxw = nypp + 2

c   ..number of rings
      nc = nxw + iopen - 1

c   ..number of rows
      nrow = nxpp + 2

c   ..first row intersecting core (on outer side)
      jprgt = nxpp - (rightcut(1) - 1) + 1

c   ..first row intersecting core (on inner side)
      jplft = nxpp - leftcut(1) + 1

c   ..number of EDGE2D mesh points (periodic points in core)
      np = ( ( nxpp + 2 ) * ( nypp + 2 ) ) + ( iopen - 1 )

c   ..dominant ion in plasma and its neutral
      ism  = -1         ! index of main ion
      ism0 = -1         ! index of main ion neutral

      do is = 0, nspp - 1
         if( nint(zn(is)).eq.1 ) then
             if( zamax(is).eq.1.0_R8 ) then
                 ism  = is
             else if( is_neutral(is) ) then
                 ism0 = is
             endif
         endif
      enddo

      if( ism.eq.-1 ) then
          write(*,9030)
          labort = .true.
          return
      endif

c   ..collect non-dominant species
      nzs     = 0        ! no. of non-dominant species
      nz      = 0        ! no. of states of all non-dominant species
                         ! (excludes neutrals)
      sqezrad = 0.0_R8   ! Radiation density for each species

      do is = 0 , nspp - 1
         if( is_neutral(is) .and. is.ne.ism0 ) then
             nzs         = nzs + 1
             izs2is(nzs) = is
             izch(nzs)   = nint(zn(is))
             zmass(nzs)  = am(is)
         endif
         if( .not.is_neutral(is) .and.
    &         nint(zn(is)).ne.nint(zn(ism)) ) then
             do ix = -1 , nxpp
                do iy = -1 , nypp
                   sqezrad(ix,iy,nzs) = sqezrad(ix,iy,nzs)
     *                                - ( rqahe(ix,iy,is)/vol(ix,iy) )
                enddo
             enddo
         endif
      enddo
      nz = nspp - nzs - 2

c   ..mach number of each species
      rz    = 0.0_R8
      pz    = 0.0_R8
      pmach = 0.0_R8

      call b2xppz( nxpp, nypp, nspp, ne, na, te, ti, pz )
      call b2xprz( nxpp, nypp, nspp, mp, am, na, rz )
      call b2xzef( nx, ny, ns, rz2, na, ne, zeff)

      do is = 0, nspp - 1
         call b2xpmc( nxpp, nypp
     *              , ua(-1,-1,is), pz(-1,-1), rz(-1,-1)
     *              , pmach(-1,-1,is) )
      enddo

c   ..geometric centres of mesh polygons & flux
      do ix = -1 , nxpp
         do iy = -1 , nypp
            cxpp(ix,iy) = 0.25_R8 *
     *                  ( crx(ix,iy,0) + crx(ix,iy,1) + crx(ix,iy,2)
     *                  + crx(ix,iy,3) )
            cypp(ix,iy) = 0.25_R8 *
     *                  ( cry(ix,iy,0) + cry(ix,iy,1) + cry(ix,iy,2)
     *                  + cry(ix,iy,3) )
            psi(ix,iy) = 0.25_R8 *
     *                  ( fpsi(ix,iy,0,0) + fpsi(ix,iy,1,0)
     *                  + fpsi(ix,iy,2,0) + fpsi(ix,iy,3,0) )
         enddo
      enddo

c   ..separatrix
      nsepx = 0
      iy    = topcut(1)

      do ix = nxpp , -1 , -1
         nsepx = nsepx + 1
         rsepx(nsepx) = crx(ix,iy,1)
         zsepx(nsepx) = cry(ix,iy,1)
      enddo

      nsepx = nsepx + 1
      rsepx(nsepx) = crx(-1,iy,0)
      zsepx(nsepx) = cry(-1,iy,0)

c   ..plasma core boundary
      nbcore = 0
      iy     = -1

      do ix = rightcut(1) - 1, leftcut(1), -1
         nbcore = nbcore + 1
         rbcore(nbcore) = crx(ix,iy,1)
         zbcore(nbcore) = cry(ix,iy,1)
      enddo

      nbcore = nbcore + 1
      rbcore(nbcore) = crx(leftcut(1),iy,0)
      zbcore(nbcore) = cry(leftcut(1),iy,0)

c   ..plasma wall boundary
      nbwall = 0

      ix = nxpp
      do iy = -1, nypp
         nbwall = nbwall + 1
         rbwall(nbwall,1) = crx(ix,iy,1)
         zbwall(nbwall,1) = cry(ix,iy,1)
         rbwall(nbwall,2) = crx(ix,iy,3)
         zbwall(nbwall,2) = cry(ix,iy,3)
      enddo

      iy = nypp
      do ix = nxpp, -1, -1
         nbwall = nbwall + 1
         rbwall(nbwall,1) = crx(ix,iy,3)
         zbwall(nbwall,1) = cry(ix,iy,3)
         rbwall(nbwall,2) = crx(ix,iy,2)
         zbwall(nbwall,2) = cry(ix,iy,2)
      enddo

      ix = -1
      do iy = nypp, -1, -1
         nbwall = nbwall + 1
         rbwall(nbwall,1) = crx(ix,iy,2)
         zbwall(nbwall,1) = cry(ix,iy,2)
         rbwall(nbwall,2) = crx(ix,iy,0)
         zbwall(nbwall,2) = cry(ix,iy,0)
      enddo

      iy = -1
      do ix = -1, leftcut(1) - 1
         nbwall = nbwall + 1
         rbwall(nbwall,1) = crx(ix,iy,0)
         zbwall(nbwall,1) = cry(ix,iy,0)
         rbwall(nbwall,2) = crx(ix,iy,1)
         zbwall(nbwall,2) = cry(ix,iy,1)
      enddo

      iy = -1
      do ix = rightcut(1), nxpp
         nbwall = nbwall + 1
         rbwall(nbwall,1) = crx(ix,iy,0)
         zbwall(nbwall,1) = cry(ix,iy,0)
         rbwall(nbwall,2) = crx(ix,iy,1)
         zbwall(nbwall,2) = cry(ix,iy,1)
      enddo

c   ..x-point
      rpx = 0.0_R8
      zpx = 0.0_R8

      xpnt: do i = 1 , nsepx
         do j = nsepx , 1 , -1
            if( i.ne.j .and.
     *          rsepx(i).eq.rsepx(j) .and.
     *          zsepx(i).eq.zsepx(j)       ) then
                rpx = rsepx(i)
                zpx = zsepx(i)
                exit xpnt
            endif
         enddo
      enddo xpnt

c   ..plasma centre

cdpc+xpb
      rmin=1e20_R8
      rmax=-1e20_R8
      imin=-2
      imax=-2
      do ix=-1,nxpp
         if(on_closed_surface(ix,-1)) then
            if(crx(ix,-1,2).lt.rmin) then
               rmin=crx(ix,-1,2)
               imin=ix
            endif
            if(crx(ix,-1,2).gt.rmax) then
               rmax=crx(ix,-1,2)
               imax=ix
            endif
         endif
      enddo
      if(imin.eq.-2.or.imax.eq.-2) then
cdpc+xpb-end
         i    = 0
         r0   = 0.0_R8
         z0   = 0.0_R8
         iy = -1
         do ix = rightcut(1) - 1 , leftcut(1) , -1
            i  = i + 1
            r0 = r0 + crx(ix,iy,1)
            z0 = z0 + cry(ix,iy,1)
            if( ix.eq.leftcut(1) ) then
               i  = i + 1
               r0 = r0 + crx(ix,iy,0)
               z0 = z0 + cry(ix,iy,0)
            endif
         enddo
         if (i.gt.0) then
           r0 = r0 / float(i)
           z0 = z0 / float(i)
           write(*,*)
     *        'Used initial approach for calculating magnetic axis'
         else
           r0 = crx(-1,-1,0)
           z0 = cry(-1,-1,0)
         endif
cdpc+xpb
      else
        r0=(rmin+rmax)/2.0_R8
        z0=(cry(imin,-1,1)+cry(imax,-1,1))/2.0_R8
         write(*,*)
     *        'Used new approach for calculating magnetic axis'
      endif
      write(*,*) ' Rmag, Zmag  = ',r0, z0

c   ..non-orthogonal points (1=Outer side, 2=Inner side)
      do ix = -1 , nxpp
         do iy = -1 , nypp
            if( abs( qz(ix,iy,0) ).lt.epsang ) then
                itagdv(ix,iy) = 0           ! Orthogonal
            else if( ix.gt.(nxpp+2)/2 ) then
                itagdv(ix,iy) = 1           ! Outer side
            else
                itagdv(ix,iy) = 2           ! Inner side
            endif
         enddo
      enddo

c   ..rho
      do iy = -1 , nypp
         do ix = -1 , nxpp
            rho(ix,iy)  = 100 * ( iy + 1 + 10 )
         enddo
      enddo

c   ..hrho
      do iy = -1 , nypp
         do ix = -1 , nxpp
            if( bottomiy(ix,iy).eq.-2 .and. topiy(ix,iy).le.nypp ) then
                hrho(ix,iy) = hy(ix,iy) * qz(ix,iy,1) /
     *           ( rho(topix(ix,iy),topiy(ix,iy)) - rho(ix,iy) )
            else if( topiy(ix,iy).eq.nypp+1 .and.
     *               bottomiy(ix,iy).ge.-1) then
                hrho(ix,iy) = hy(ix,iy) * qz(ix,iy,1) /
     *           ( rho(ix,iy) - rho(bottomix(ix,iy),bottomiy(ix,iy)) )
            else if( topiy(ix,iy).le.nypp .and.
     *               bottomiy(ix,iy).ge.-1) then
                hrho(ix,iy) = 2.0_R8 * hy(ix,iy) * qz(ix,iy,1) /
     *           ( rho(topix(ix,iy),topiy(ix,iy)) -
     *             rho(bottomix(ix,iy),bottomiy(ix,iy)) )
            else
                hrho(ix,iy) = 0.0_R8
            endif
         enddo
      enddo

c   ..theta
      call calc_theta()

c   ..htheta
      do iy = -1 , nypp
         do ix = -1 , nxpp
            if( leftix(ix,iy).eq.-2 .and. rightix(ix,iy).le.nxpp ) then
                htheta(ix,iy) = hx(ix,iy) /
     *           ( theta(ix,iy) - theta(rightix(ix,iy),rightiy(ix,iy)) )
            else if( rightix(ix,iy).eq.nxpp+1 .and.
     *               topix(ix,iy).ge.-1 ) then
                htheta(ix,iy) = hx(ix,iy) /
     *           ( theta(leftix(ix,iy),leftiy(ix,iy)) - theta(ix,iy) )
            else if( rightix(ix,iy).le.nxpp .and.
     *               topix(ix,iy).ge.-1 ) then
                htheta(ix,iy) = 2.0_R8 * hx(ix,iy) /
     *           ( theta(leftix(ix,iy) ,leftiy(ix,iy)) -
     *             theta(rightix(ix,iy),rightiy(ix,iy)) )
            else
                htheta(ix,iy) = 0.0_R8
            endif
         enddo
      enddo

c   ..B-tor
      do ix = -1 , nxpp
         do iy = -1 , nypp
            bfi(ix,iy) = bb(ix,iy,2)
         enddo
      enddo

c   ..Pitch angle
      do ix = -1 , nxpp
         do iy = -1 , nypp
            sh(ix,iy) = bb(ix,iy,0) / bb(ix,iy,3)
         enddo
      enddo

c   ..neutral molecular density (currently awaiting Eirene module)
      do ix = -1 , nxpp
         do iy = -1 , nypp
            dm(ix,iy) = 0.0_R8
         enddo
      enddo

c   ..formats
      return
 9000 format( ' *** error : nuclear charges not in ascending order !!'
     *      / '             is     =' , i6 , 'zamin(is) =' , i6
     *      / '             zn(is) =' , i6 , '   zn(is) =' , i6 )
 9010 format( ' *** error : atomic charges not in ascending order !!'
     *      / '             is     =' , i6 , '   zn(is) =' , i6
     *      / '          zamin(is) =' , i6 , 'zamax(is) =' , i6 )
 9020 format( ' *** error : no neutral for species !!'
     *      / '             is     =' , i6 , '   zn(is) =' , i6 )
 9030 format( ' *** error : no hydrogen species.' )

      end subroutine set_b2mod_ppout
c
c=======================================================================
      subroutine calc_theta

      use b2mod_geo
      use b2mod_indirect

      implicit none

c     ------------------------------------------------------------------
c     CALC_THETA sets up coordinate theta which is orthogonal to
c     the mesh rings
c     ------------------------------------------------------------------

c   ..locals
      integer ix , iy , iv , it , ixs , ixn , iys ,
     *        ivbott , ivbotn , ivtopt , ivtopn

      integer, allocatable ::
     *  iermat(:,:)

      real (kind=R8) ::
     *  xt1(3),
     *  xt2(3),
     *  xt3(3),
     *  efac,
     *  zmin,
     *  frac,
     *  t

      allocate(iermat(-1:nxpp,-1:nypp))

c   ..data
      data efac/1.0_R8/, zmin/-2.0_R8/
      data xt1/3*0._R8/, xt2/3*0._R8/, xt3/3*0._R8/

c   ..set up unit normals to field lines
      do ix = -1 , nxpp
         do iy = -1 , nypp
            xnorm(ix,iy,0) = ( cry(ix,iy,2) + cry(ix,iy,0) -
     *                         cry(ix,iy,3) - cry(ix,iy,1) )
     *                       / hx(ix,iy) / 2.0_R8
            xnorm(ix,iy,1) = ( crx(ix,iy,3) + crx(ix,iy,1) -
     *                         crx(ix,iy,2) - crx(ix,iy,0) )
     *                       / hx(ix,iy) / 2.0_R8
         enddo
      enddo

c   ..copy the grid to start extending the grid
      do iv = 0 , 3
         do ix = -1 , nxpp
            do iy = -1 , nypp
               crxe(ix,iy,iv) = crx(ix,iy,iv)
               crye(ix,iy,iv) = cry(ix,iy,iv)
            enddo
         enddo
      enddo

c   ..add a row of extra cells at each target (West then East)
c     by linear extrapolation
      do it = 1 , 2

         if( it.eq.1 ) then
             ix     = -1         ! target row
             ixn    = -2         ! new row
             ivbott =  1         ! for ixn : bottom vertex at target
             ivbotn =  0         !           bottom vertex opposite
             ivtopt =  3         !           top    vertex at target
             ivtopn =  2         !           top    vertex opposite
         else
             ix     =  nxpp      ! target row
             ixn    =  nxpp+1    ! new row
             ivbott =  0         ! for ixn : bottom vertex at target
             ivbotn =  1         !           bottom vertex opposite
             ivtopt =  2         !           top    vertex at target
             ivtopn =  3         !           top    vertex opposite
         endif

         do iy = -1 , nypp
            crxe(ixn,iy,ivbott) = crx(ix,iy,ivbotn)
            crye(ixn,iy,ivbott) = cry(ix,iy,ivbotn)
            crxe(ixn,iy,ivtopt) = crx(ix,iy,ivtopn)
            crye(ixn,iy,ivtopt) = cry(ix,iy,ivtopn)
            crxe(ixn,iy,ivbotn) = -efac*crx(ix,iy,ivbott) +
     *                              (1.0_R8 + efac)*crx(ix,iy,ivbotn)
            crye(ixn,iy,ivbotn) = -efac*cry(ix,iy,ivbott) +
     *                              (1.0_R8 + efac)*cry(ix,iy,ivbotn)
            crxe(ixn,iy,ivtopn) = -efac*crx(ix,iy,ivtopt) +
     *                              (1.0_R8 + efac)*crx(ix,iy,ivtopn)
            crye(ixn,iy,ivtopn) = -efac*cry(ix,iy,ivtopt) +
     *                              (1.0_R8 + efac)*cry(ix,iy,ivtopn)
            xnorm(ixn,iy,0)     = xnorm(ix,iy,0)
            xnorm(ixn,iy,1)     = xnorm(ix,iy,1)
         enddo

      enddo

c   ..add a further row of extra cells at each target (West then East)
c     by extrapolating, parallel to sepx outer boundary, to zmin
      iys = topcut(1)

      do it = 1 , 2

         if( it.eq.1 ) then
             ix     = -2         ! target row
             ixn    = -3         ! new row
             ivbott =  1         ! for ixn : bottom vertex at target
             ivbotn =  0         !           bottom vertex opposite
             ivtopt =  3         !           top    vertex at target
             ivtopn =  2         !           top    vertex opposite
         else
             ix     =  nxpp+1    ! target row
             ixn    =  nxpp+2    ! new row
             ivbott =  0         ! for ixn : bottom vertex at target
             ivbotn =  1         !           bottom vertex opposite
             ivtopt =  2         !           top    vertex at target
             ivtopn =  3         !           top    vertex opposite
         endif

         t = ( crye(ix,iys,ivbott)+crye(ix,iys,ivtopt) -
     *           crye(ix,iys,ivbotn)-crye(ix,iys,ivtopn) ) /
     *             ( crxe(ix,iys,ivbott)+crxe(ix,iys,ivtopt) -
     *                 crxe(ix,iys,ivbotn)-crxe(ix,iys,ivtopn) )

         do iy = -1 , nypp
            crxe(ixn,iy,ivbott) = crxe(ix,iy,ivbotn)
            crye(ixn,iy,ivbott) = crye(ix,iy,ivbotn)
            crxe(ixn,iy,ivtopt) = crxe(ix,iy,ivtopn)
            crye(ixn,iy,ivtopt) = crye(ix,iy,ivtopn)
            crxe(ixn,iy,ivbotn) = crxe(ix,iy,ivbotn) +
     *                              (zmin - crye(ix,iy,ivbotn))/t
            crye(ixn,iy,ivbotn) = zmin
            crxe(ixn,iy,ivtopn) = crxe(ix,iy,ivtopn) +
     *                              (zmin - crye(ix,iy,ivtopn))/t
            crye(ixn,iy,ivtopn) = zmin
            xnorm(ixn,iy,0)     = xnorm(ix,iys,0)
            xnorm(ixn,iy,1)     = xnorm(ix,iys,1)
         enddo

      enddo

c   ..for all original cells find intercept of ring normal which
c     passes through cell centre with sepx top boundary and
c     consequently assign theta
      do iy = -1 , nypp
         do ix = -1 , nxpp
            call follow_normal(ix,iy,ixs,frac,iermat(ix,iy))
            if( iermat(ix,iy).eq.0 )then
                theta(ix,iy) = (1.0_R8-frac)*thetasepx(ixs) +
     &                          frac*thetasepx(ixs+1)
            endif
         enddo
      enddo

c   ..cater for failed cells by interpolation
      do iy = 0 , nypp - 1
         do ix = -1 , nxpp
            if (bottomiy(ix,iy).ne.-2 .and.
     &          bottomiy(ix,iy).ne.nypp+1) then
              if( iermat(ix,iy).gt.0 .and.
     &            iermat(bottomix(ix,iy),bottomiy(ix,iy)).eq.0 .and.
     &            iermat(topix(ix,iy),topiy(ix,iy)).eq.0 )then
                xt1(1) = cxpp(ix,iy)
                xt1(2) = cypp(ix,iy)
                xt2(1) = cxpp(bottomix(ix,iy),bottomiy(ix,iy))
                xt2(2) = cypp(bottomix(ix,iy),bottomiy(ix,iy))
                call vsub(xt1,xt2,xt3)
                frac   = vlen(xt3)
                xt1(1) = cxpp(topix(ix,iy),topiy(ix,iy))
                xt1(2) = cypp(topix(ix,iy),topiy(ix,iy))
                call vsub(xt1,xt2,xt3)
                frac   = frac/vlen(xt3)
                theta(ix,iy) =
     &           (1.0_R8-frac)*theta(bottomix(ix,iy),bottomiy(ix,iy)) +
     &                   frac*theta(topix(ix,iy),topiy(ix,iy))
                iermat(ix,iy) = 0
              endif
            endif
            if (leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nxpp+1) then
              if( iermat(ix,iy).gt.0 .and.
     &            iermat(leftix(ix,iy) ,leftiy(ix,iy)).eq.0 .and.
     &            iermat(rightix(ix,iy),rightiy(ix,iy)).eq.0)then
                xt1(1) = cxpp(ix,iy)
                xt1(2) = cypp(ix,iy)
                xt2(1) = cxpp(leftix(ix,iy),leftiy(ix,iy))
                xt2(2) = cypp(leftix(ix,iy),leftiy(ix,iy))
                call vsub(xt1,xt2,xt3)
                frac   = vlen(xt3)
                xt1(1) = cxpp(rightix(ix,iy),rightiy(ix,iy))
                xt1(2) = cypp(rightix(ix,iy),rightiy(ix,iy))
                call vsub(xt1,xt2,xt3)
                frac   = frac/vlen(xt3)
                theta(ix,iy) =
     .           (1.0_R8-frac)*theta(leftix(ix,iy),leftiy(ix,iy)) +
     &                   frac*theta(rightix(ix,iy),rightiy(ix,iy))
                iermat(ix,iy) = 0
              endif
            endif
         enddo
      enddo

c   ..final check
      do iy = -1, nypp
         do ix = -1 , nxpp
            if( iermat(ix,iy).gt.0 )then
               write(*,9000) ix,iy
c               labort = .true.
            endif
         enddo
      enddo

c   ..formats
      return
 9000 format( ' *** error : failed to create theta for ix,iy =',2i4)

      end subroutine calc_theta
c
c=======================================================================
      subroutine follow_normal(ix,iy,ixs,frac,ier)

      use b2mod_geo
      use b2mod_indirect

      implicit none

      integer ix , iy , ixs, ier

      real (kind=R8) ::
     *  frac

c     ------------------------------------------------------------------
c     FOLLOW_NORMAL - Starting from cell centre (ix,iy) follow the ring
c     normal to find the intercept with the sepx ring top boundary. The
c     intercept is labelled by the cell row number ixs and the fraction
c     of distance along the cell top boundary
c     ------------------------------------------------------------------

c   ..locals
      integer iva(0:2) , ivb(0:2) , idir
      integer ixn , iyn , ixdbg , iydbg , nbnd, n

      real (kind=R8) ::
     *  ro(3),
     *  roalt(3),
     *  ra(3),
     *  ral(3),
     *  rb(3),
     *  rbl(3),
     *  rbsubra(3),
     *  rosubra(3),
     *  xdelta(3),
     *  xlambda(3),
     *  xt1(3),
     *  xt2(3),
     *  xt3(3),
     *  ang,
     *  t(0:2),
     *  fracn(0:2),
     *  s(0:2),
     *  smin

      data ixdbg/-10/, iydbg/-1/
      data      ro/3*0._R8/,   roalt/3*0._R8/
      data      ra/3*0._R8/,     ral/3*0._R8/
      data      rb/3*0._R8/,     rbl/3*0._R8/, rbsubra/3*0._R8/
      data rosubra/3*0._R8/,  xdelta/3*0._R8/, xlambda/3*0._R8/
      data     xt1/3*0._R8/,     xt2/3*0._R8/,     xt3/3*0._R8/

c   ..direction to sepx (-1 towards bottom, +1 towards top)
c   ..iva(n),ivb(n) specify the vertices corresponding to
c   ..the cell boundary n that can be intersected
      if( iy.le.topcut(1) ) then
          idir   =  1
          iva(0) =  2
          ivb(0) =  3
      else
          idir   = -1
          iva(0) =  0
          ivb(0) =  1
      endif
      iva(1) = 0
      ivb(1) = 2
      iva(2) = 1
      ivb(2) = 3

c   ..starting point
      ier   = 0
      ixn   = ix
      iyn   = iy
      ro(1) = 0.25_R8 *
     *       ( crxe(ixn,iyn,0) + crxe(ixn,iyn,1) + crxe(ixn,iyn,2) +
     *           crxe(ixn,iyn,3) )
      ro(2) = 0.25_R8 *
     *       ( crye(ixn,iyn,0) + crye(ixn,iyn,1) + crye(ixn,iyn,2) +
     *           crye(ixn,iyn,3) )

      ral(1) = 0._R8
      ral(2) = 0._R8
      rbl(1) = 0._R8
      rbl(2) = 0._R8
c
      if( ix.eq.ixdbg .and. iy.eq.iydbg ) then
          write(*,*) ' '
          write(*,*) 'b2mod_ppout -> follow_normal'
          write(*,*) ' Tracking from ix, iy :',ix,iy
      endif

c   ..track until sepx outer boundary is intersected
      do while ( iyn.ne.(topcut(1)+max(0,idir)) )

         if (ixn.lt.lbound(xnorm,1).or.
     &       ixn.gt.ubound(xnorm,1).or.
     &       iyn.lt.lbound(xnorm,2).or.
     &       iyn.gt.ubound(xnorm,2)) then
             write(*,*) 'Attempt to access xnorm with ',ixn,iyn
             ixn=max(lbound(xnorm,1),min(ubound(xnorm,1),ixn))
             iyn=max(lbound(xnorm,2),min(ubound(xnorm,2),iyn))
             write(*,*) 'ixn, iyn reset to ',ixn,iyn
         endif

         xlambda(1) = idir*xnorm(ixn,iyn,0)
         xlambda(2) = idir*xnorm(ixn,iyn,1)
         if( ix.eq.ixdbg .and. iy.eq.iydbg ) then
             write(*,*) ' '
             write(*,*) 'cell    :',ixn,iyn
             write(*,'(a,1p,2e20.12)') ' ro      :',ro(1),ro(2)
             write(*,'(a,1p,2e20.12)') ' roalt   :',roalt(1),roalt(2)
             write(*,'(a,1p,2e14.6)') ' xlambda :',xlambda(1),xlambda(2)
         endif

c   ..find intersections with cell boundary segments
         do n=0,2
            ra(1)    = crxe(ixn,iyn,iva(n))
            ra(2)    = crye(ixn,iyn,iva(n))
            rb(1)    = crxe(ixn,iyn,ivb(n))
            rb(2)    = crye(ixn,iyn,ivb(n))
            t(n)     = -1.e10_R8
            s(n)     = -1.e10_R8
            fracn(n) = -1.e10_R8

c         ..omit boundary containing previous intersect
            if( ra(1).ne.ral(1) .or. ra(2).ne.ral(2) .or.
     *          rb(1).ne.rbl(1) .or. rb(2).ne.rbl(2) ) then

c             ..intersection between r()=ro() + xlambda()*s
c             ..             and     r()=ra() +  xdelta()*t
                call vsub(ro,ra,rosubra)
                call vsub(rb,ra,rbsubra)
                call vunit(rbsubra,xdelta)
                call vcross(rosubra,xlambda,xt1)
                call vcross(rosubra,xdelta,xt3)
                call vcross(xdelta,xlambda,xt2)
                ang      = vdot(xt2,xt2)
                t(n)     = vdot(xt1,xt2) / ang
                s(n)     = vdot(xt3,xt2) / ang
                fracn(n) = t(n)/vlen(rbsubra)
                if( ix.eq.ixdbg .and. iy.eq.iydbg ) then
                    write(*,*) ' ..boundary n   :',n
                    write(*,*) '   iva,ivb      :',iva(n),ivb(n)
                    write(*,'(a,1p,3e14.6)') '    ra           :',
     *                              ra(1),ra(2)
                    write(*,'(a,1p,3e14.6)') '    rb           :',
     *                              rb(1),rb(2)
                    write(*,'(a,1p,3e14.6)') '    s,t,frac     :',
     *                              s(n),t(n),fracn(n)
                endif

            endif

         enddo

c      ..find which boundary is intersected within segment
         nbnd = -1
         n    =  0
         do while ( nbnd.eq.-1 .and. n.le.2 )
           if( fracn(n).gt.0._R8 .and. fracn(n).lt.1._R8 .and.
     *         s(n).gt.0._R8 ) then
               nbnd = n
           endif
           n = n + 1
         enddo

c      ..no intersect found - choose best candidate in
c      ..right direction
         if( nbnd.eq.-1 ) then
             smin = 1.e10_R8
             do n=0,2
               if( s(n).gt.0._R8 .and. s(n).lt.smin ) then
                  nbnd = n
                  smin = s(n)
               endif
             enddo
         endif

         if( ix.eq.ixdbg .and. iy.eq.iydbg ) then
             write(*,*) '==> selecting n :',nbnd
         endif

c      ..update position and cell
         if (nbnd.ge.0 .and. nbnd.le.2) then
           ral(1) = crxe(ixn,iyn,iva(nbnd))
           ral(2) = crye(ixn,iyn,iva(nbnd))
           rbl(1) = crxe(ixn,iyn,ivb(nbnd))
           rbl(2) = crye(ixn,iyn,ivb(nbnd))
           ro(1)  = ro(1) + xlambda(1)*s(nbnd)
           ro(2)  = ro(2) + xlambda(2)*s(nbnd)
           if( ix.eq.ixdbg .and. iy.eq.iydbg ) then
             call vsub(rbl,ral,rbsubra)
             call vunit(rbsubra,xdelta)
             roalt(1) = ral(1) + xdelta(1)*t(nbnd)
             roalt(2) = ral(2) + xdelta(2)*t(nbnd)
           endif
         endif
         if( nbnd.eq.0 ) then
             iyn = iyn + idir
         else if( nbnd.eq.1 ) then
             ixn = ixm1(ixn,iyn)
         else if( nbnd.eq.2 ) then
             ixn = ixp1(ixn,iyn)
         endif

c      ..next cell
         if( nbnd.eq.-1 ) then
             ier = 1
             return
         endif
      enddo

c   ..check the final intersection
      if( nbnd.eq.0 ) then
          ixs  = ixn
      else
          ixs  = 0
      endif
      frac = fracn(nbnd)

      return
      end subroutine follow_normal
c
c=======================================================================
      real (kind=R8) function thetasepx(ix)

      use b2mod_indirect

      implicit none

      integer ix

c     ------------------------------------------------------------------
c     THETASEPX - Calculates angles of each point of separatrix.
c     ------------------------------------------------------------------

      thetasepx = 100.0_R8 * (-ix + nxpp + 2 + 10.0_R8)

      return
      end function thetasepx
c
c=======================================================================
      integer function ixm1(ix,iy)

      use b2mod_indirect

      implicit none

      integer ix , iy

c     ------------------------------------------------------------------
c     IXM1 - For a given cell (ix,iy) return the ix value of the
c     adjacent cell in the left direction
c     ------------------------------------------------------------------

      if (ix.lt.lbound(leftix,1).or.
     &    ix.gt.ubound(leftix,1).or.
     &    iy.lt.lbound(leftiy,2).or.
     &    iy.gt.ubound(leftiy,2)) then
          write(*,*) 'Trying to access the left neighbour of ',ix,iy
        ixm1 = ix - 1
      else
        ixm1=leftix(ix,iy)
      endif

      return
      end function ixm1
c
c=======================================================================
      integer function ixp1(ix,iy)

      use b2mod_indirect

      implicit none

      integer ix , iy

c     ------------------------------------------------------------------
c     IXP1 - For a given cell (ix,iy) return the ix value of the
c     adjacent cell in the right direction
c     ------------------------------------------------------------------

      if (ix.lt.lbound(rightix,1).or.
     &    ix.gt.ubound(rightix,1).or.
     &    iy.lt.lbound(rightiy,2).or.
     &    iy.gt.ubound(rightiy,2)) then
          write(*,*) 'Trying to access the right neighbour of ',ix,iy
        ixp1 = ix + 1
      else
        ixp1=rightix(ix,iy)
      endif

      return
      end function ixp1
c
c=======================================================================
      subroutine vsub(a,b,c)

      implicit none

      integer, parameter ::
     *  ndim = 3

      real (kind=R8) ::
     *  a(ndim), b(ndim), c(ndim)

c     ------------------------------------------------------------------
c     VSUB - Calculate the vector : c = a - b
c     ------------------------------------------------------------------

      c = a - b

      return
      end subroutine vsub
c
c=======================================================================
      subroutine vunit(a,b)

      implicit none

      integer, parameter ::
     *  ndim = 3

      real (kind=R8) ::
     *  a(ndim), b(ndim)

c     ------------------------------------------------------------------
c     VUNIT - Calculate the vector : b = a / |a|
c     ------------------------------------------------------------------

c   ..locals
      integer i

      real (kind=R8) ::
     *  xmoda

c   ..calculation
      xmoda = sqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )

      do i = 1 , ndim
         b(i) = a(i) / xmoda
      enddo

      return
      end subroutine vunit
c
c=======================================================================
      subroutine vcross(a,b,c)

      implicit none

      integer, parameter ::
     *  ndim = 3

      real (kind=R8) ::
     *  a(ndim), b(ndim), c(ndim)

c     ------------------------------------------------------------------
c     VCROSS - Calculate the cross product : c = a X b
c     ------------------------------------------------------------------

      c(1) = +a(2)*b(3) - a(3)*b(2)
      c(2) = -a(1)*b(3) + a(3)*b(1)
      c(3) = +a(1)*b(2) - a(2)*b(1)

      return
      end subroutine vcross
c
c=======================================================================
      function vdot(a,b)

      implicit none

      integer, parameter ::
     *  ndim = 3

      real (kind=R8) ::
     *  a(ndim), b(ndim)

      real (kind=R8) ::
     *  vdot

c     ------------------------------------------------------------------
c     VDOT - Calculate the dot product : a . b
c     ------------------------------------------------------------------

      vdot = a(1)*b(1) + a(2)*b(2) + a(3)*b(3)

      return
      end function vdot
c
c=======================================================================
      function vlen(a)

      implicit none

      integer, parameter ::
     *  ndim = 3

      real (kind=R8) ::
     *  a(ndim)

      real (kind=R8) ::
     *  vlen

c     ------------------------------------------------------------------
c     VLEN - Calculate the length
c     ------------------------------------------------------------------

      vlen = sqrt( a(1)*a(1) + a(2)*a(2) + a(3)*a(3) )

      return
      end function vlen
c
c=======================================================================
      subroutine write_b2mod_ppout(nwrite)

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_B2MOD_PPOUT is a front-end to a series of routines for
c     writing data to the postprocessor transfer file.
c     ------------------------------------------------------------------
      call write_pproc_header(nwrite)

      call write_pproc_version(nwrite)

      call write_pproc_control(nwrite)

      call write_pproc_control_form(nwrite)

      if( labort ) return

      call write_pproc_geometry(nwrite)

      call write_pproc_timestep(nwrite)

      call write_pproc_timetrace(nwrite)

      call write_pproc_boundary(nwrite)

      call write_pproc_eof(nwrite)

      return
      end subroutine write_b2mod_ppout
c
c=======================================================================
      subroutine write_pproc_header(nwrite)

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_PPROC_HEADER writes header section of postprocessor
c     transfer file.
c     ------------------------------------------------------------------

c   ..locals
      character*8 tcode , tmachid, tishot , tdate

c   ..process
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*CATALOGUE IDENTIFIER')

      call ppname( tcode , tmachid, tishot , tdate )
      call ppcstr(nwrite,' '//tcode//' '//tmachid//tishot//' '//tdate )

      return
      end subroutine write_pproc_header
c
c=======================================================================
      subroutine write_pproc_version(nwrite)

      implicit none

      integer nwrite
#ifndef NO_GETENV
      integer lenval, ierror
#ifndef USE_PXFGETENV
      intrinsic get_environment_variable
#endif
#endif
*!!   (We avoid to declare, as external or intrinsic, the subroutine
*     getenv; it is external on some systems (SunOS, Unicos, ...) and
*     intrinsic on others (AIX, HP-UX, ...).  If getenv is declared
*     external on AIX or HP-UX it returns some non-printing garbage,
*     without error message (AIX/XLF Version 2.2, HP-UX Release 9.0))
c     ------------------------------------------------------------------
c     WRITE_PPROC_VERSION writes version section of postprocessor
c     transfer file.
c     ------------------------------------------------------------------

      character*7  grpsolps,locsolps
      character*16 version
      logical streql
      external streql

      grpsolps=' '
      locsolps=' '
#ifndef NO_GETENV
#ifdef USE_PXFGETENV
      CALL PXFGETENV ('GRPSOLPS', 0, grpsolps, lenval, ierror)
      CALL PXFGETENV ('LOCSOLPS', 0, locsolps, lenval, ierror)
#else
      call get_environment_variable('GRPSOLPS',
     .   status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('GRPSOLPS',
     .   value=grpsolps)
      call get_environment_variable('LOCSOLPS',
     .   status=ierror, length=lenval)
      if (ierror.eq.0) call get_environment_variable('LOCSOLPS',
     .   value=locsolps)
#endif
#endif
      if (.not.streql(grpsolps,' ').and.
     &    .not.streql(locsolps,' ')) then
        version = grpsolps//'_'//locsolps
      else
        version = ' '
      endif
c
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*CODE VERSION')

      call ppcstr(nwrite,' VERSION :'//version )

      return
      end subroutine write_pproc_version
c
c=======================================================================
      subroutine write_pproc_control(nwrite)

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_B2MOD_CONTROL writes data which controls the run to the
c     postprocessor transfer file.
c     ------------------------------------------------------------------

      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*CONTROL DATA')

      return
      end subroutine write_pproc_control
c
c=======================================================================
      subroutine write_pproc_control_form(nwrite)

      use b2mod_b2cmpa

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_B2MOD_CONTROL_FORM writes formatted control data for the
c     run to the postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer izs , is

      character*80 adasidh
      character*80, allocatable ::
     *   adasidz(:)

      character*8 cbuf8
      character*2 yearh, yearz

      character*80 cbuff

C   ..allocate
      allocate(adasidz(1:nspp))

c   ..initialise
      adasidh = adasid
      adasidz = adasid

c   ..header
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*CONTROL DATA (FORMATTED)')

c   ..ADAS main ion data
      write( yearh,'(i2)' ) isyear(ism)

      call ppcstr( nwrite, ' ADASIDH:' // adasidh )
      call ppcstr( nwrite, ' YEARH  :' // yearh )

c   ..ADAS impurity data
      do izs = 1 , max( nzs, 1 )
         is = izs2is(izs)

         cbuf8 = ' ' // ppaddsp( 'ADASIDZ', izs )
         call ppcstr( nwrite, cbuf8 // ':' // adasidz(izs) )

         write(yearz,'(i2)') isyear(is)
         cbuf8 = ppaddsp( ' YEARZ', izs )
         call ppcstr( nwrite, cbuf8 // ':' // yearz )
      enddo

      return
      end subroutine write_pproc_control_form
c
c=======================================================================
      subroutine write_pproc_geometry(nwrite)

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_B2MOD_GEOMETRY is a front-end to writing geometry data to
c     the postprocessor transfer file.
c     ------------------------------------------------------------------

      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*GEOMETRY')

      call write_gmtry(nwrite)

      return
      end subroutine write_pproc_geometry
c
c=======================================================================
      subroutine write_pproc_timestep(nwrite)

      use b2mod_types
      use b2mod_plasma
      use b2mod_constants
      use b2mod_elements
      use b2mod_b2cmpa

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_PPROC_TIMESTEP writes time step section of the
c     postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer izs, is , i

      character*2 cbuf2
      character*6 cbuf6
      character*8 cbuf8

      character*80 cbuff

c   ..header
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*TIME STEP DATA')
      write(cbuff,1000) 0.0_R8 , 1
      call ppcstr(nwrite,cbuff)

c   ..plasma state
      if( ism.ge.0 ) then
          call ppwure( nwrite
     *               , 'HMASS' , 'ION MASS' , 'AMU'
     *               , am(ism) , 1.0_R8 , 1  , 0 )
          itemp(1) = int((zamin(ism)+zamax(ism))/2.0_R8)
          call ppwuin( nwrite
     *               , 'HCH'   , 'ION ATOMIC NUMBER' , ''
     *               , itemp   , 1 , 1  , 0 )
      endif

      call ppwuin( nwrite
     *           , 'NZ' , 'NO. OF IMP. CHARGE STATES' , ''
     *           , izch , 1 , max(nzs,1), 0 )

      if( nz.gt.0 ) then
          call ppwure( nwrite
     *               , 'ZMASS' , 'IMP.MASS' , 'AMU'
     *               , zmass   , 1.0_R8 , nzs , 0 )
          call ppwuin( nwrite
     *               , 'ZCH'   , 'IMP.ATOMIC NUMBER' , ''
     *               , izch    , 1      , nzs , 0 )
      endif

      if( ism.ge.0 ) then
          call ppwure_conv_grid( nwrite
     *                         , 'DEN' , 'ION DENSITY' , 'M-3'
     *                         , na(-1,-1,ism) , 1.0_R8 , 1 , 0 )
      endif

      call ppwure_conv_grid( nwrite
     *                     , 'DENEL' , 'EL. DENSITY' , 'M-3'
     *                     , ne     , 1.0_R8 , 1 , 0 )
      call ppwure_conv_grid( nwrite
     *                     , 'TEV' , 'ION TEMPERATURE' , 'EV'
     *                     , ti     , 1.0_R8/ev , 1 , 0 )
      call ppwure_conv_grid( nwrite
     *                     , 'TEVE' , 'EL. TEMPERATURE' , 'EV'
     *                     , te     , 1.0_R8/ev , 1 , 0 )
      call ppwure_conv_grid( nwrite
     *                     , 'ZEFF' , 'Z EFFECTIVE' , ' '
     *                     , zeff   , 1.0_R8    , 1 , 0 )

      if( ism.ge.0 ) then
          call ppwure_conv_grid( nwrite
     *                         , 'VPI' , 'ION PAR.VEL' , 'M S-1'
     *                         , ua(-1,-1,ism) , -1.0_R8 , 1 , 1 )
          call ppwure_conv_grid( nwrite
     *                         , 'PMACH' , 'MACH NUMBER (IONS)' , ''
     *                         , pmach(-1,-1,ism) , -1.0_R8 , 1 , 1 )
      endif

c   ..neutral atoms
      if( ism0.ge.0 ) then
          call ppwure_conv_grid( nwrite
     *                         , 'DA' , 'NEUT. ATOM DENSITY' , 'M-3'
     *                         , na(-1,-1,ism0) , 1.0_R8 , 1 , 0 )
      endif

c   ..neutral molecules
      call ppwure_conv_grid( nwrite
     *                     , 'DM' , 'NEUT. MOL. DENSITY' , 'M-3'
     *                     , dm   , 1.0_R8 , 1 , 0 )

c   ..neutral impurities
      if( nz.gt.0 ) then
          i = 0
          do is = 0 , nspp-1
             if( is_neutral(is) .and.
     &           nint(zn(is)).ne.nint(zn(ism)) ) then
                 i = i + 1
                 cbuf2 = is_codes(is)
                 call chcase(1,cbuf2)
                 cbuf8 = ppaddsp( 'DZ', i )
                 call ppwure_conv_grid( nwrite , cbuf8
     *                                , cbuf2 // '- NEUT. IMP. DENSITY'
     *                                , 'M-3'
     *                                , na(-1,-1,is) , 1.0_R8 , 1 , 0 )
             endif
          enddo
      endif

c   ..impurity species
      do izs = 1 , nzs
          is = izs2is(izs)
          cbuf2 = is_codes(is)
          call chcase(1,cbuf2)
          cbuf8 = ppaddsp( 'SQEZRAD', izs )
          call ppwure_conv_grid( nwrite , cbuf8
     *                         , cbuf2 // '- IMPURITY RADIATION'
     *                         , 'W M-3'
     *                         , sqezrad(-1,-1,izs) , 1.0_R8 , 1 , 0 )
      enddo

c   ..impurity charged states
      if( nz.gt.0 ) then
          i = 0
          do is = 0 , nspp-1
             if( .not.is_neutral(is) .and.
     *            nint(zn(is)).ne.nint(zn(ism)) ) then
                 i = i + 1
                 write(cbuf6,1010) 'IMP ' , i

                 write(cbuf8,1010) 'DENZ' , i
                 call ppwure_conv_grid( nwrite , cbuf8
     *                                , cbuf6//' DENSITY' , 'M-3'
     *                                , na(-1,-1,is) , 1.0_R8 , 1 , 0 )

                 write(cbuf8,1010) 'VPZ' , i
                 call ppwure_conv_grid( nwrite , cbuf8
     *                                , cbuf6//' PAR.VEL' , 'M S-1'
     *                                , ua(-1,-1,is) , -1.0_R8 , 1 , 1 )
             endif
          enddo
      endif

c   ..formats
      return
 1000 format( '-TIME  :' , E11.5 , ' STEP  :' , I8 )
 1010 format( a , i3.3 )

      end subroutine write_pproc_timestep
c
c=======================================================================
      subroutine write_pproc_timetrace(nwrite)

      use b2mod_types
      use b2mod_constants
      use b2mod_b2cmpa

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_PPROC_TIMETRACE writes time trace section of the
c     postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      real (kind=R8) ::
     *  tvec(1)

c   ..header
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*TIME TRACES')

c   ..traces
      tvec(1) = 0.0_R8
      call ppwure( nwrite
     *           , 'TVEC', 'TRACE TIMES' , 'MS'
     *           , tvec  , 1.0_R8 , 1    , 0 )

c   ..formats
      return
      end subroutine write_pproc_timetrace
c
c=======================================================================
      subroutine write_pproc_boundary(nwrite)

      use b2mod_types
      use b2mod_plasma
      use b2mod_constants
      use b2mod_b2cmpa

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_PPROC_BOUNDARY writes boundary section of the
c     postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..header
      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*BOUNDARY')

c   ..plasma wall boundary
      call ppwure( nwrite
     *           , 'RPLAM1' , 'R-PLASMA BOUNDARY 1' , 'M'
     *           , rbwall(1,1), 1.0_R8 , nbwall  , 0 )
      call ppwure( nwrite
     *           , 'ZPLAM1' , 'Z-PLASMA BOUNDARY 1' , 'M'
     *           , zbwall(1,1), 1.0_R8 , nbwall  , 0 )
      call ppwure( nwrite
     *           , 'RPLAM2' , 'R-PLASMA BOUNDARY 2' , 'M'
     *           , rbwall(1,2), 1.0_R8 , nbwall  , 0 )
      call ppwure( nwrite
     *           , 'ZPLAM2' , 'Z-PLASMA BOUNDARY 2' , 'M'
     *           , zbwall(1,2), 1.0_R8 , nbwall  , 0 )

c   ..plasma core boundary
      call ppwure( nwrite
     *           , 'RBCORE' , 'R-CORE BOUNDARY' , 'M'
     *           , rbcore(1), 1.0_R8 , nbcore  , 0 )
      call ppwure( nwrite
     *           , 'ZBCORE' , 'Z-CORE BOUNDARY' , 'M'
     *           , zbcore(1), 1.0_R8 , nbcore  , 0 )

      return
      end subroutine write_pproc_boundary
c
c=======================================================================
      subroutine write_pproc_eof(nwrite)

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_PPROC_EOF writes terminal section of postprocessor
c     transfer file.
c     ------------------------------------------------------------------

      call ppcstr(nwrite,'*')
      call ppcstr(nwrite,'*EOF')

      return
      end subroutine write_pproc_eof
c
c=======================================================================
      subroutine write_gmtry(nwrite)

      use b2mod_geo

      implicit none

      integer nwrite

      character*80 cbuff

c     ------------------------------------------------------------------
c     WRITE_GMTRY writes geometry data of postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..name of geometry file
      call ppcstr( nwrite, ' '//gmtry_file() )

c   ..basic values
      call ppcstr( nwrite, 'NP,IOPEN,NXW,NC,NROW,JPRGT,JPLFT,REFCT' )
      if( iform.eq.0 ) then
          write(nwrite) np, iopen, nxw, nc, nrow, jprgt, jplft, 0
      else
          cbuff = ' '
          write(cbuff,1000) np, iopen, nxw, nc, nrow, jprgt, jplft, 0
          call ppcstr( nwrite, cbuff )
      endif

c   ..x-point
      rtemp(1) = rpx
      call ppwure( nwrite
     *           , 'RPX'   , 'X-POINT R-COORD' , 'M'
     *           , rtemp   , 1.0_R8 , 1 , 0 )
      rtemp(1) = zpx
      call ppwure( nwrite
     *           , 'ZPX'   , 'X-POINT Z-COORD' , 'M'
     *           , rtemp   , 1.0_R8 , 1 , 0 )

c   ..plasma centre
      rtemp(1) = r0
      call ppwure( nwrite
     *           , 'R0'   , 'AXIAL R-COORD' , 'M'
     *           , rtemp  , 1.0_R8 , 1 , 0 )
      rtemp(1) = z0
      call ppwure( nwrite
     *           , 'Z0'   , 'AXIAL Z-COORD' , 'M'
     *           , rtemp  , 1.0_R8 , 1 , 0 )

c   ..separatrix
      call ppwure( nwrite
     *           , 'RSEPX' , 'R FOR SEPARATRIX' , 'M'
     *           , rsepx   , 1.0_R8 , nsepx , 0 )
      call ppwure( nwrite
     *           , 'ZSEPX' , 'Z FOR SEPARATRIX' , 'M'
     *           , zsepx   , 1.0_R8 , nsepx , 0 )

c   ..(r,z) coordinates (convert to postprocessor grid)
      call ppwure_conv_grid( nwrite
     *                     , 'RMESH' , 'R-COORD' , 'M'
     *                     , cxpp    , 1.0_R8 , 0 , 0 )
      call ppwure_conv_grid( nwrite
     *                     , 'ZMESH' , 'Z-COORD' , 'M'
     *                     , cypp    , 1.0_R8 , 0 , 0 )

c   ..non-orthoginal mapping
      call ppwuin_conv_grid( nwrite
     *                     , 'ITAGDV' , '' , ''
     *                     , itagdv   , 1  , 0 )

c   ..flux
      call ppwure_conv_grid( nwrite
     *                     , 'PSI' , 'FLUX' , 'WEBERS'
     *                     , psi   , 1.0_R8 , 1 , 0 )

c   ..rho
      call ppwure_conv_grid( nwrite
     *                     , 'RHO' , 'RHO' , ''
     *                     , rho   , 1.0_R8 , 0 , 0 )

c   ..theta
      call ppwure_conv_grid( nwrite
     *                     , 'THETA' , 'THETA' , ''
     *                     , theta   , 1.0_R8 , 0 , 0 )

c   ..hrho (1.0E+02 factor due to postprocessor patch)
      call ppwure_conv_grid( nwrite
     *                     , 'HRHO' , 'HRHO' , ''
     *                     , hrho   , 1.0E+02_R8 , 0 , 0 )

c   ..htheta (1.0E+02 factor due to postprocessor patch)
      call ppwure_conv_grid( nwrite
     *                     , 'HTETA' , 'HTETA' , ''
     *                     , htheta  , 1.0E+02_R8  , 0 , 0 )

c   ..b-toroidal
      call ppwure_conv_grid( nwrite
     *                     , 'BFI' , 'B-TOROIDAL' , 'T'
     *                     , bfi   , 1.0_R8  , 0 , 0 )

c   ..pitch angle
      call ppwure_conv_grid( nwrite
     *                     , 'SH' , 'BPOL/BTOT' , ''
     *                     , sh   , 1.0_R8  , 0 , 0 )

c   ..polygons
      call write_poly( nwrite )
      return

c   ..formats
 1000 format( 7i4 , '   0' )

      end subroutine write_gmtry
c
c=======================================================================
      subroutine ppcstr(nwrite,cstr)

      integer nwrite
      character*(*) cstr

c     ------------------------------------------------------------------
c     PPCSTR is a general routine to write a line to the postprocessor
c     transfer file.
c     ------------------------------------------------------------------

c   ..locals
      character*80 cbuff

c   ..initialisation
      cbuff = ' '
      cbuff = cstr

c   ..process
      if( iform.eq.0 ) then
          write(nwrite) cbuff
      else
          write(nwrite,'(a)') cbuff
      endif

      return
      end subroutine ppcstr
c
c=======================================================================
      subroutine ppname(tcode,tmachid,tishot,tdate)

      use b2mod_geo

      implicit none

c     ------------------------------------------------------------------
c     PPNAME is a general routine to extract B2 code name, shot and
c     run date.
c     ------------------------------------------------------------------

      character*(*) tcode , tmachid, tishot , tdate

c   ..locals
      integer        i , isep1 , isep2
      character*120  buffer
      character*1    csep

      data csep/ '.' /

c   ..initialisation
      tcode  = 'solps'

      isep1  = 0
      isep2  = 0

c   ..find first two separators in geometry file name
      buffer = gmtry_file()

      seps: do i = 1 , len(buffer)
          if( buffer(i:i).eq.csep ) then
              if( isep1.eq.0 ) then
                  isep1 = i
              else
                  isep2 = i
                  exit seps
              endif
          endif
      enddo seps

c   ..shot (JET geometry file is of the form "jet.shot.etc.etc")
      tmachid = buffer(1:3)
      call chcase(0,tmachid)
      if( tmachid.ne.'jet' .or. isep1.eq.0 .or. isep2.eq.0 ) then
          write(*,9000) gmtry_file()
          labort = .true.
      else
          tishot  = buffer(isep1+1:isep2-1)
      endif

c   ..process
      call date(buffer)

      if( buffer(1:1).eq.' ' ) buffer(1:1) = '0'
      if( buffer(8:8).eq.' ' ) buffer(8:8) = '0'

      tdate = buffer(4:6) // buffer(1:2) // buffer(8:9)
      call chcase(0,tdate)
      return

c   ..formats
 9000 format( ' *** error : non-JET geometry file!!'
     *      / '             ' , a )

      end subroutine ppname
c
c=======================================================================
      subroutine ppwure_conv_grid( nwrite
     *                           , name   , desc , units
     *                           , rdat   , rsf  , iextrp , istag )

      implicit none

      integer nwrite , iextrp , istag

      character*(*) name , desc , units

      real (kind=R8) ::
     *  rdat(-1:nxpp,-1:nypp),
     *  rsf

c     ------------------------------------------------------------------
c     PPWURE_CONV_GRID converts real*8 data (on B2 grid) to a data
c     (on EDGE2D grid) then writes to the postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      real (kind=R8), allocatable ::
     *  rconv(:)

      integer, allocatable ::
     *  nj(:),
     *  kory(:,:)

c   ..initialise
      allocate( rconv(1:np) )
      allocate( nj(1:nc) )
      allocate( kory(1:nc,1:nrow) )

c   ..convert B2 mesh to EDGE2D mesh (necessary for postprocessor)
      call ppre_conv_grid( rdat , rconv , iextrp , nj , kory )

c   ..stagger EDGE2D mesh (if required)
      if( istag.ne.0 ) then
          call ppre_stag_grid( rconv, rconv , nj , kory )
      endif

c   ..process
      call ppwure( nwrite
     *           , name   , desc , units
     *           , rconv  , rsf  , np    , istag )

      return
      end subroutine ppwure_conv_grid
c
c=======================================================================
      subroutine ppwuin_conv_grid( nwrite
     *                           , name   , desc , units
     *                           , idat   , isf  , istag )

      implicit none

      integer nwrite , istag

      character*(*) name , desc , units

      integer ::
     *  idat(-1:nxpp,-1:nypp),
     *  isf

c     ------------------------------------------------------------------
c     PPWUIN_CONV_GRID converts integer data (on B2 grid) to a data
c     (on EDGE2D grid) then writes to the postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer, allocatable ::
     *  iconv(:)

c   ..initialise
      allocate( iconv(1:np) )

c   ..convert B2 mesh to EDGE2D mesh (necessary for postprocessor)
      call ppin_conv_grid( idat , iconv )

c   ..process
      call ppwuin( nwrite
     *           , name   , desc , units
     *           , iconv  , isf  , np    , istag )

      return
      end subroutine ppwuin_conv_grid
c
c=======================================================================
      subroutine ppwure( nwrite
     *                 , name   , desc , units
     *                 , rdat   , rsf  , ndata , istag )

      implicit none

      integer nwrite , ndata , istag

      character*(*) name , desc , units

      real (kind=R8) ::
     *  rdat(1:ndata),
     *  rsf

c     ------------------------------------------------------------------
c     PPWURE writes real*8 data (on EDGE2D grid) to postprocessor
c     trabsfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer i
      character*1   ftype , fstag
      character*8   fname , funits
      character*32  fdesc
      character*80  cbuff

c   ..initialisation
      fname  = name
      fdesc  = desc
      funits = units
      ftype  = 'R'

      if( istag.eq.0 ) then
          fstag = ' '
      else
          fstag = 'S'
      endif

c   ..write data
      if( iform.eq.0 ) then
          cbuff = ' '
          write(cbuff,1000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite) cbuff
          write(nwrite) ( sngl( rdat(i)*rsf ), i=1,ndata)
      else if( iform.eq.1 ) then
          write(nwrite,2000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite,2010) ( sngl( rdat(i)*rsf ), i=1,ndata)
      else
          write(nwrite,2000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite,3010) ( sngl( rdat(i)*rsf ), i=1,ndata)
      endif
      return

c   ..formats
 1000 format( 3a,i8,2a,1x,a )
 2000 format( 2(a,1x),a,i8,1x,3(a,1x) )
 2010 format( 10z8 )
 3010 format( 1p,7e11.3 )

      end subroutine ppwure
c
c=======================================================================
      subroutine ppwuin( nwrite
     *                 , name   , desc , units
     *                 , idat   , isf  , ndata , istag )

      implicit none

      integer nwrite , ndata , istag

      character*(*) name , desc , units

      integer ::
     *  idat(1:ndata),
     *  isf

c     ------------------------------------------------------------------
c     PPWUIN writes integer data (on EDGE2D grid) to postprocessor
c     trabsfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer i
      character*1   ftype , fstag
      character*8   fname , funits
      character*32  fdesc
      character*80  cbuff

c   ..initialisation
      fname  = name
      fdesc  = desc
      funits = units
      ftype  = 'I'

      if( istag.eq.0 ) then
          fstag = ' '
      else
          fstag = 'S'
      endif

c   ..write data
      if( iform.eq.0 ) then
          cbuff = ' '
          write(cbuff,1000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite) cbuff
          write(nwrite) ( ( idat(i)*isf ), i=1,ndata)
      else if( iform.eq.1 ) then
          write(nwrite,2000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite,2010) ( ( idat(i)*isf ), i=1,ndata)
      else
          write(nwrite,2000) '#', fname, ftype, ndata, fdesc, funits,
     *                       fstag
          write(nwrite,3010) ( ( idat(i)*isf ), i=1,ndata)
      endif
      return

c   ..formats
 1000 format( 3a,i8,2a,1x,a )
 2000 format( 2(a,1x),a,i8,1x,3(a,1x) )
 2010 format( 10z8 )
 3010 format( 7i10 )

      end subroutine ppwuin
c
c=======================================================================
      subroutine ppre_conv_grid( rdat , rconv , iextrp , nj , kory )

      use b2mod_indirect

      implicit none

      real (kind=R8) ::
     *  rdat(-1:nxpp,-1:nypp),
     *  rconv(1:np)

      integer ::
     *  iextrp,
     *  nj(1:nc),
     *  kory(1:nc,1:nrow)

c     ------------------------------------------------------------------
c     PPRE_CONV_GRID converts real*8 data (on B2 grid) to a data
c     (on EDGE2D grid).  Linearly extrapolate end-points of SOL and
c     private rings.
c     ------------------------------------------------------------------

c   ..locals
      integer ix , iy , k, kfirst, klast, iring

c   ..initialise
      k     = 0
      iring = 0

c   ..core rings
      do iy = bottomcut(1) , topcut(1) - 1
         iring = iring + 1
         nj(iring) = 0
         do ix = rightcut(1) - 1 , leftcut(1) , -1
            k  = k + 1
            rconv(k) = rdat(ix,iy)
            nj(iring) = nj(iring) + 1
            kory(iring,nj(iring)) = k
         enddo
         k = k + 1
         rconv(k) = rdat(rightcut(1)-1,iy)  ! periodic point
         nj(iring) = nj(iring) + 1
         kory(iring,nj(iring)) = k
      enddo

c   ..sol rings
      do iy = topcut(1) , nypp
         iring = iring + 1
         nj(iring) = 0
         do ix = nxpp , -1 , -1
            k  = k + 1
            rconv(k) = rdat(ix,iy)
            nj(iring) = nj(iring) + 1
            kory(iring,nj(iring)) = k
            if( ix.eq.nxpp ) kfirst = k
            if( ix.eq.-1   ) klast  = k
         enddo
         if( iextrp.ne.0 ) then
             rconv(kfirst) = rconv(kfirst)
     *                     + ( rconv(kfirst) - rconv(kfirst+1) )
             rconv(klast)  = rconv(klast)
     *                     + ( rconv(klast)  - rconv(klast-1)  )
         endif
      enddo

c   ..private rings
      do iy = bottomcut(1) , topcut(1) - 1
         iring = iring + 1
         nj(iring) = 0
         do ix = nxpp , -1 , -1
            if( ix.ge.rightcut(1) .or. ix.lt.leftcut(1) ) then
                k  = k + 1
                rconv(k) = rdat(ix,iy)
                nj(iring) = nj(iring) + 1
                kory(iring,nj(iring)) = k
                if( ix.eq.nxpp ) kfirst = k
                if( ix.eq.-1   ) klast  = k
            endif
         enddo
         if( iextrp.ne.0 ) then
             rconv(kfirst) = rconv(kfirst)
     *                     + ( rconv(kfirst) - rconv(kfirst+1) )
             rconv(klast)  = rconv(klast)
     *                     + ( rconv(klast)  - rconv(klast-1)  )
         endif
      enddo

      return
      end subroutine ppre_conv_grid
c
c=======================================================================
      subroutine ppin_conv_grid( idat , iconv )

      use b2mod_indirect

      implicit none

      integer ::
     *  idat(-1:nxpp,-1:nypp),
     *  iconv(1:np)

c     ------------------------------------------------------------------
c     PPIN_CONV_GRID converts integer data (on B2 grid) to a data
c     (on EDGE2D grid).
c     ------------------------------------------------------------------

c   ..locals
      integer ix , iy , k

c   ..initialise
      k  = 0

c   ..core rings
      do iy = bottomcut(1) , topcut(1) - 1
         do ix = rightcut(1) - 1 , leftcut(1) , -1
            k  = k + 1
            iconv(k) = idat(ix,iy)
         enddo
         k = k + 1
         iconv(k) = idat(rightcut(1)-1,iy)  ! periodic point
      enddo

c   ..sol rings
      do iy = topcut(1) , nypp
         do ix = nxpp , -1 , -1
            k  = k + 1
            iconv(k) = idat(ix,iy)
         enddo
      enddo

c   ..private rings
      do iy = bottomcut(1) , topcut(1) - 1
         do ix =  nxpp , -1 , -1
            if( ix.ge.rightcut(1) .or. ix.lt.leftcut(1) ) then
                k  = k + 1
                iconv(k) = idat(ix,iy)
            endif
         enddo
      enddo

      return
      end subroutine ppin_conv_grid
c
c=======================================================================
      subroutine ppre_stag_grid( runstag , rstag , nj , kory )

      use b2mod_indirect

      implicit none

      real (kind=R8) ::
     *  runstag(1:np),
     *  rstag(1:np)

      integer ::
     *  nj(nc),
     *  kory(1:nc,1:nrow)

c     ------------------------------------------------------------------
c     PPRE_STAG_GRID converts real*8 unstaggered data to REAL*8
c     staggered data.
c     ------------------------------------------------------------------

c   ..locals
      integer i , j , k , km1

c   ..calculate staggered values for each ring
      do i = 1 , nc
         do j = nj(i) , 2 , -1
            k   = kory(i,j)
            km1 = kory(i,j-1)
            rstag(k) = 0.5_R8 * ( runstag(k) + runstag(km1) )
         enddo
         k = kory(i,1)
         if(i.le.topcut(1)-bottomcut(1)) then
            km1 = kory(i,nj(i))
            rstag(k) = 0.5_R8 * ( runstag(k) + runstag(km1) )
         else
            rstag(k) = 0.0_R8
         endif
      enddo

      return
      end subroutine ppre_stag_grid
c
c=======================================================================
      character*8 function ppaddsp( name , izs )

      implicit none

      character*(*) name
      integer       izs

c     ------------------------------------------------------------------
c     PPADDSP turn label into label(izs).
c     ------------------------------------------------------------------

c   ..locals
      integer     i , lname
      character*2 cbuf2

c   ..initialise
      loop: do i = len(name) , 1 , -1
        if( name(i:i).ne.' ' ) exit loop
      enddo loop

      lname = i

c   ..process
      if( nzs.gt.1 ) then
          if( lname.ge.7 ) lname = 5
          write(cbuf2,'(a,i1)') '_',izs
          ppaddsp = name(1:lname)//cbuf2
      else
          ppaddsp = name
      endif

      return
      end function ppaddsp
c
c=======================================================================
      subroutine write_poly( nwrite )

      use b2mod_indirect

      implicit none

      integer nwrite

c     ------------------------------------------------------------------
c     WRITE_POLY writes polygons data to postprocessor transfer file.
c     ------------------------------------------------------------------

c   ..locals
      integer, allocatable ::
     *   korpg(:),
     *   nvert(:)

      real (kind=R8), allocatable ::
     *   rvert(:,:),
     *   zvert(:,:)

      integer i, ix , iy , k , npoly , npoly0

c   ..allocate
      allocate(korpg(1:np))
      allocate(nvert(1:np))

      allocate(rvert(1:5,1:np))
      allocate(zvert(1:5,1:np))

c   ..initialise
      do k = 1 , np
         korpg(k) = 0
         nvert(k) = 0

         do i = 1 , 5
            rvert(i,k) = 0.0_R8
            zvert(i,k) = 0.0_R8
         enddo
      enddo

      k     = 0
      npoly = 0

c   ..core rings
      do iy = bottomcut(1) , topcut(1) - 1
         do ix = rightcut(1) - 1 , leftcut(1) , -1
            k = k + 1
            if( iy.ne.bottomcut(1) ) then
                npoly = npoly + 1
                if( ix.eq.rightcut(1) - 1 ) npoly0 = npoly

                korpg(k) = npoly
                call mkvert( ix, iy, npoly, nvert, rvert, zvert )
            endif
         enddo
         k = k + 1
         if( iy.ne.bottomcut(1) ) korpg(k) = npoly0
      enddo

c   ..sol rings
      do iy = topcut(1) , nypp
         do ix = nxpp , -1 , -1
            k = k + 1
            if( iy.ne.nypp .and. ix.ne.nxpp .and. ix.ne.-1 ) then
                npoly = npoly + 1
                korpg(k) = npoly
                call mkvert( ix, iy, npoly, nvert, rvert, zvert )
            endif
         enddo
      enddo

c   ..private rings
      do iy = bottomcut(1) , topcut(1) - 1
         do ix =  nxpp , -1 , -1
            if( ix.ge.rightcut(1) .or. ix.lt.leftcut(1) ) then
                k = k + 1
                if( iy.ne.bottomcut(1) .and.
     *              ix.ne.nxpp .and. ix.ne.-1 ) then
                    npoly = npoly + 1
                    korpg(k) = npoly
                    call mkvert( ix, iy, npoly, nvert, rvert, zvert )
                endif
            endif
         enddo
      enddo

c   ..write to file
      call ppwuin( nwrite
     *           , 'KORPG' , '' , ''
     *           , korpg   , 1 , np , 0 )
      call ppwuin( nwrite
     *           , 'NVERTP' , '' , ''
     *           , nvert    , 1 , npoly , 0 )

      call ppwure( nwrite
     *           , 'RVERTP' , '' , ''
     *           , rvert    , 1.0_R8 , 5*npoly , 0 )
      call ppwure( nwrite
     *           , 'ZVERTP' , '' , ''
     *           , zvert    , 1.0_R8 , 5*npoly , 0 )

      return
      end subroutine write_poly
c
c=======================================================================
      subroutine mkvert( ix , iy , npoly , nvert , rvert , zvert )

      use b2mod_geo

      implicit none

      integer ix , iy , npoly

      integer ::
     *  nvert(np)
      real (kind=R8) ::
     *  rvert(5,np),
     *  zvert(5,np)

c     ------------------------------------------------------------------
c     MKPOLY make polygon for an EDGE2D grid
c     ------------------------------------------------------------------

c   ..locals
      integer iv
      integer ivmap(4)

c   ..mapping EDGE2D vertices to B2 (ivmap(edge2d) = b2)
      data ivmap / 1 , 3 , 2 , 0/

c   ..build polygon
      nvert(npoly) = 4     ! number of vertices

      do iv = 1 , nvert(npoly)
         rvert(iv,npoly) = crx(ix,iy,ivmap(iv))
         zvert(iv,npoly) = cry(ix,iy,ivmap(iv))
      enddo

c   ..centre of polygon
      rvert(5,npoly) = cxpp(ix,iy)
      zvert(5,npoly) = cypp(ix,iy)

      return
      end subroutine mkvert
c
c=======================================================================
      character*120 function gmtry_file

c     ------------------------------------------------------------------
c     GMTRY_FILE returns name of the geometry file.
c     ------------------------------------------------------------------

c   ..locals
      integer i, l1, l2

c   ..initialise
      i1 = 0
      i2 = 0

c   ..form geometry file name
      do i = len(labgeo), 1, -1
         if( labgeo(i:i).ne.' ' .and. i2.eq.0 ) i2 = i
         if( labgeo(i:i).eq.' ' .and. i2.gt.0 ) then
             i1 = i
             exit
         endif
      enddo

      if( i2.gt.0 ) then
          gmtry_file = labgeo(i1+1:i2)
      else
          gmtry_file = ' '
      endif
c conn
      return
      end function gmtry_file

c
c=======================================================================
#endif
      subroutine read_b2mod_ppout(type, ninp, label)

      use b2mod_version
      implicit none

      integer ninp
      character*(*) type
      character*(*) label

      external cfruin

c     ------------------------------------------------------------------
c     READ_B2MOD_PPOUT_DATA reads data for use by ppout.
c     (in) type  = rates --- read atomic rates data for use by ppout.
c     (in) ninp  = stream number
c     (in) label = special label data
c     ------------------------------------------------------------------

c   ..initialise
      isyear = iyeard

c   ..atomic rates method
      if( type.eq.'rates' ) then
          if( label(61:64).eq.'ADAS' ) then
#ifdef JET
              call cfruin(ninp, nspp, isyear, 'ppout')
#else
              if (b2frates_version.ge.'01.001.024') then
                call cfruin(ninp, nspp, isyear, 'ppout')
              endif
#endif
          endif
      endif

      return
      end subroutine read_b2mod_ppout
c
c=======================================================================
      end module b2mod_ppout
