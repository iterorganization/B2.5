*.specification

      subroutine b2trfl (ix,iy,Ainv,zeffne,amti,zetap,cF_lim) 
      use b2mod_types
      use b2mod_b2cmpa
      use b2mod_b2cmpt
      use b2mod_indirect
      use b2mod_constants     
      use b2mod_zhfrtf
      use b2mod_b2zhco , only: zhflcorr
      use b2mod_subsys
      implicit none
*   ..input arguments (unchanged on exit)
      integer ix,iy
      real (kind=R8) :: zeffne,zetap,amti(0:nspecies-1),
     & Ainv(0:2*nspecies-1,0:2*nspecies-1)
*   ..output arguments (unspecified on entry)
      real (kind=R8) :: 
     & cF_lim(0:nspecies-1,0:nspecies-1)
*   ..common blocks
*-----------------------------------------------------------------------
*.documentation
*-----------------------------------------------------------------------
*.declarations
*   ..local variables
      integer i,j
      real (kind=R8) :: coef,t0,t1
*-----------------------------------------------------------------------
*.computation

*   ..subprogram start-up calls
      call subini ('b2trfl')
*   ..main
      if (cflim(4).ne.0.0_R8.and..not.on_closed_surface(ix,iy)) then
        coef=2.5_R8*zetap/mp/(zeffne**2)
!       print *,coef
        do i=0, nspecies-1
          t0=abs(sum(Ainv(i,0:nspecies-1)*nal(ix,iy,0:nspecies-1)*
     & gtalc(ix,iy,0:nspecies-1)))
!       print *,t0*coef
!       print *,cflim(4)
          t1=zhflcorr*cflim(4)/sqrt(amti(i))
        
!       print *,t1
!       print *,cflim(4)
!       print *,sqrt(amti(i))

          cF_lim(i,0:nspecies-1)=1/(1+coef*t0/t1)
        enddo
      else
        cF_lim(0:nspecies-1,0:nspecies-1)=1
      endif 

!       print *, 'checkpoint_cF_lim'
!       print *, 'ix,iy:',ix,iy
!       print *, 'cF_lim:'
!       do, i=0,nspecies-1
!         print *, (cF_lim(i,j), j=0,nspecies-1 )
!       enddo
!       print *

      call subend ()
      return


      end subroutine b2trfl
