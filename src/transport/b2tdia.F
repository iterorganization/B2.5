      subroutine b2tdia(nCv, nFc, nVx, ns, geo, mpg,
     &   pz, fna, uadia, facdrift, fchdia)       !srv 16.19.17 {
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_constants
      use b2us_geo
      use b2us_map
      implicit none
*   ..input arguments (unchanged on exit)
      integer nCv, nFc, nVx, ns
      type (geometry), intent (in) :: geo
      type (mapping), intent (in) :: mpg
      real (kind=R8) ::
     &  pz(nCv), facdrift(nFc), fna(nFc,0:1,0:ns-1), 
     &  uadia(nFc,0:1,0:ns-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchdia(nFc,0:1)
*-----------------------------------------------------------------------
*.documentation
*
*  1. purpose
*
*     B2TDIA computes the X and Y components of the divergent part of
*     the diamagnetic current
*
*  8. further comments
*
*-----------------------------------------------------------------------
*.declarations

*   ..local variables
      integer iFc
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm, wght(nFc,2), gOnedBsq(nFc,0:1),
     &  pzbbf(nFc)
      real (kind=R8), save :: dia_cur = 1.0e0_R8                         !srv 17.05.99
#ifdef CLASSICAL_CURRENT
cwdk  option CLASSICAL_CURRENT not yet converted!
      integer is
      real (kind=R8) ::
     *   gampvd(-1:nx,-1:ny), gamrvd(-1:nx,-1:ny),
     *   gampvdc(-1:nx,-1:ny)
      external intvertex
#endif
*   ..procedures
      intrinsic max, abs
      real (kind=R8) ::
     *   smin, smax
      external subini, subend, xertst
      external b2xvsg, b2xvff, b2xvfx
*-----------------------------------------------------------------------
*.computation

* ..preliminaries
*   ..subprogram start-up calls
      call subini ('b2tdia')
*   ..set internal parameters on first call
      if (ncall.eq.0) then                                               !srv
       call ipgetr ('b2tfhe_dia_cur', dia_cur)                           !srv
      end if
*   ..extensive tests on first few calls
      if (ncall.lt.3) then
*   ..test nx, ny, ns
       call xertst (0.le.nCv.and.0.le.nFc, 'faulty argument nCv, nFc')
       call xertst (1.le.ns, 'faulty argument ns')
      endif
      
      facdriftm = maxval(facdrift)
      if (facdriftm.ne.0.0_R8 .and. dia_cur.ne.0.0e0_R8) then            !srv 16.10.17
*     ..interpolate to faces
        wght = 1.0_R8
        call intface(nCv,nFc,mpg%fcCv,wght,pz*geo%cvBb(:,2),pzbbf)

*     ..compute gradients of OnedBsq on cell faces
        call grad(nCv, nFc, nVx, 1, geo, mpg, geo%cvOnedBsq,
     &    geo%vxOnedBsq, gOnedBsq)

*     ..DIVERGENT PART OF DIAMAGNETIC CURRENT
        fchdia(:,0) =  
     &     dia_cur*facdrift*pzbbf*gOnedBsq(:,1)*geo%fcS*geo%fcQalf(:,0)
        fchdia(:,1) = 
     &    -dia_cur*facdrift*pzbbf*gOnedBsq(:,0)*geo%fcS*geo%fcQalf(:,1)

#ifdef CLASSICAL_CURRENT
crfs inertia terms
        call xerrab ('Not yet adapted to extended grids format')

        gampvd=0.0_R8
        gamrvd=0.0_R8
        do is = 0, ns-1
          if (.not.is_neutral(is)) then
            do iy = -1, ny
              do ix = -1, nx
                if(leftix(ix,iy).ne.-2.and.rightix(ix,iy).ne.nx+1.and.
     &             bottomiy(ix,iy).ne.-2.and.topiy(ix,iy).ne.ny+1.and.
     &             region(ix,iy,0).ne.0)
     >            gampvd(ix,iy)=gampvd(ix,iy)+am(is)*mp*
     1            0.5_R8*(fna(ix,iy,0,0,is)*uadia(ix,iy,0,0,is)+
     2            fna(rightix(ix,iy),rightiy(ix,iy),0,0,is)
     2            *uadia(rightix(ix,iy),rightiy(ix,iy),0,0,is))
              enddo
            enddo
            do ix = -1, nx
              do iy = -1, ny
                if(rightix(ix,iy).ne.nx+1 .and. topiy(ix,iy).ne.ny+1
     &       .and. leftix(ix,iy).ne.-2 .and. bottomiy(ix,iy).ne.-2
     &       .and. region(ix,iy,0).ne.0) then
cxpb  Must check for case of the island X-point
                  if(topiy(topix(ix,iy),topiy(ix,iy)).ne.iy) then     	 !xpb
cxpb  Standard case
                    gamrvd(ix,iy)=gamrvd(ix,iy)+am(is)*mp*
     1               0.5_R8*(fna(ix,iy,1,1,is)+
     2               fna(topix(ix,iy),topiy(ix,iy),1,1,is))*
     1               0.5_R8*(uadia(ix,iy,0,0,is)+
     2               uadia(rightix(ix,iy),rightiy(ix,iy),0,0,is))
                  else                                                   !xpb
cxpb  Island X-point cells
                    gamrvd(ix,iy)=gamrvd(ix,iy)+am(is)*mp*
     1               fna(ix,iy,1,1,is)*                                    !xpb
     2               0.5_R8*(uadia(ix,iy,0,0,is)+
     2               uadia(rightix(ix,iy),rightiy(ix,iy),0,0,is))
                  endif
                endif                                           	 !xpb
              enddo
            enddo
          endif
        enddo
        call intcor(nx,ny,vol,gampvd,gampvdc)
        do iy = -1, ny
          do ix = -1, nx
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1 .and.
     &         bottomiy(ix,iy).ne.-2 .and. topiy(ix,iy).ne.ny+1 .and.
     &         region(ix,iy,0).ne.0) then
              fchdia(ix,iy,1,1) = fchdia(ix,iy,1,1)+facdrift(ix,iy)*
     1         0.5_R8/wbbv(ix,iy,3)*
     2         ((gampvdc(rightix(ix,iy),rightiy(ix,iy))
     3          -gampvdc(ix,iy))/
     6         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy)))+
     2         (gamrvd(ix,iy)
     3          -gamrvd(bottomix(ix,iy),bottomiy(ix,iy)))/
     8         (hy1(ix,iy)+hy1(bottomix(ix,iy),bottomiy(ix,iy))))
            endif
          enddo
        enddo
! /* CLASSICAL_CURRENT */
#endif
      else
        fchdia = 0.0_R8
      endif

* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tdia

      end subroutine b2tdia                                              !srv 16.10.17 }
