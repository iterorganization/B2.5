      subroutine b2tdia(nx, ny, vol, gs, bb, pz, facdrift, fchdia)       !srv 16.19.17 {
      use b2mod_types
      use b2mod_indirect
      use b2mod_geo_corner
      implicit none
*   ..input arguments (unchanged on exit)
      integer nx, ny
      real (kind=R8) ::
     &  vol(-1:nx,-1:ny), gs(-1:nx,-1:ny,0:2), bb(-1:nx,-1:ny,0:3),
     &  pz(-1:nx,-1:ny), facdrift(-1:nx,-1:ny)
*   ..output arguments (unspecified on entry)
      real (kind=R8) ::
     &  fchdia(-1:nx,-1:ny,0:1)
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
      integer n2, ix, iy
      integer, save :: ncall = 0
      real (kind=R8) ::
     &  facdriftm
      real (kind=R8), save :: dia_cur = 1.0e0_R8                         !srv 17.05.99
*   ..procedures
      intrinsic max, abs
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
*   ..set n2
      n2 = (nx+2)*(ny+2)
*   ..test nx, ny
      call xertst (0.le.nx.and.0.le.ny, 'faulty argument nx, ny')
      
      fchdia = 0.0e0_R8
      facdriftm = maxval(facdrift)
      if (facdriftm.ne.0.0_R8 .and. dia_cur.ne.0.0e0_R8) then            !srv 16.10.17
*       X-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT

!         do iy=-1,ny
        do iy=0,ny-1
          do ix=0,nx
            if(leftix(ix,iy).ne.-2) then                                 !srv 26.11.02
             fchdia(ix,iy,0)=fchdia(ix,iy,0)+
     &                       dia_cur*facdrift(ix,iy)* ! signmf*          !srv 25.04.00 05.02.05 11.01.13 16.10.17
     &      (pz(ix,iy)*bb(ix,iy,2)+                                      !srv 31.03.00 18.01.02
     &        pz(leftix(ix,iy),leftiy(ix,iy))                            !srv 31.03.00
     &       *bb(leftix(ix,iy),leftiy(ix,iy),2))*0.5e0_R8                !srv 31.03.00
     &      *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))
     &      /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))
     &      *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))
!     &       (ne(ix,iy)+ne(leftix(ix,iy),leftiy(ix,iy)))/2.0e0_R8       !srv
!     &      *(te(ix,iy)+te(leftix(ix,iy),leftiy(ix,iy))                 !srv
!     &       +ti(ix,iy)+ti(leftix(ix,iy),leftiy(ix,iy)))/2.0e0_R8       !srv
!     &      *(bb(ix,iy,2)+bb(leftix(ix,iy),leftiy(ix,iy),2))/2.0e0_R8   !srv
!     &      *(vol(ix,iy)+vol(leftix(ix,iy),leftiy(ix,iy)))              !srv
!     &      /(gs(ix,iy,2)+gs(leftix(ix,iy),leftiy(ix,iy),2))            !srv
!     &      *(OnedBsqc(topix(ix,iy),topiy(ix,iy))-OnedBsqc(ix,iy))      !srv
            endif                                                        !srv 26.11.02
          enddo
        enddo

*       Y-COMPONENT OF DIVERGENT PART OF DIAMAGNETIC CURRENT

        do iy=0,ny
          do ix=0,nx-1
            if(leftix(ix,iy).ne.-2 .and. rightix(ix,iy).ne.nx+1)then     !srv 26.11.02
!           do ix=-1,nx                                                  !srv
!             fchdia(ix,iy,1)=fchdia(ix,iy,1) - dia_cur*signmf*          !srv 25.04.00 21.04.00 21.03.02
             fchdia(ix,iy,1)= - dia_cur*facdrift(ix,iy)* ! signmf*       !srv 25.04.00 21.04.00 21.03.02 11.01.13
!     &      (vol(bottomix(ix,iy),bottomiy(ix,iy))*pz(ix,iy)*bb(ix,iy,2)+ !srv 31.03.00 18.01.02 { interpolaion for bottom
!     &       vol(ix,iy)*pz(bottomix(ix,iy),bottomiy(ix,iy))             !srv 31.03.00
!     &       *bb(bottomix(ix,iy),bottomiy(ix,iy),2))/                   !srv 31.03.00
!     &      (vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))           !srv 18.01.02 }
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))          !srv
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))        !srv
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))  !srv
     &       (pz(ix,iy)*bb(ix,iy,2)+                                     !srv 31.03.00 18.01.02 {
     &        pz(bottomix(ix,iy),bottomiy(ix,iy))                        !srv 31.03.00
     &       *bb(bottomix(ix,iy),bottomiy(ix,iy),2))*0.5e0_R8            !srv 31.03.00 18.01.02 }
     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))           !srv
     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))         !srv
     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))   !srv

!     &      -(ne(ix,iy)+ne(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8   !srv
!     &      *(te(ix,iy)+te(bottomix(ix,iy),bottomiy(ix,iy))             !srv
!     &       +ti(ix,iy)+ti(bottomix(ix,iy),bottomiy(ix,iy)))/2.0e0_R8   !srv
!     &      *(bb(ix,iy,2)+bb(bottomix(ix,iy),bottomiy(ix,iy),2))        !srv
!     ^       /2.0e0_R8                                                  !srv
!     &      *(vol(ix,iy)+vol(bottomix(ix,iy),bottomiy(ix,iy)))          !srv
!     &      /(gs(ix,iy,2)+gs(bottomix(ix,iy),bottomiy(ix,iy),2))        !srv
!     &      *(OnedBsqc(rightix(ix,iy),rightiy(ix,iy))-OnedBsqc(ix,iy))  !srv
            endif                                                        !srv 26.11.02
          enddo                                                          !srv
        enddo                                                            !srv
      endif
* ..return
      ncall = ncall+1
      call subend ()
      return
*-----------------------------------------------------------------------
*.end b2tdia

      end subroutine b2tdia                                              !srv 16.10.17 }
